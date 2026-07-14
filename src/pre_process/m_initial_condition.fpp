!>
!! @file
!! @brief Contains module m_initial_condition

!> @brief Assembles initial conditions by layering prioritized patches via constructive solid geometry
module m_initial_condition

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_variables_conversion
    use m_icpp_patches
    use m_assign_variables
    use m_perturbation
    use m_chemistry
    use m_boundary_conditions

    implicit none

    type(ic_context) :: ic  !< Initial-condition state (fields, bc types, patch ids)

contains

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_initial_condition_module

        integer :: i, j, k, l

        allocate (ic%q_prim_vf(1:sys_size))
        allocate (ic%q_cons_vf(1:sys_size))

        do i = 1, sys_size
            allocate (ic%q_prim_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
            allocate (ic%q_cons_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end do

        if (chemistry) then
            allocate (ic%q_T_sf%sf(idwbuff(1)%beg:idwbuff(1)%end,idwbuff(2)%beg:idwbuff(2)%end,idwbuff(3)%beg:idwbuff(3)%end))
        end if

        allocate (ic%patch_id_fp(0:m,0:n,0:p))

        if (qbmm .and. .not. polytropic) then
            allocate (pb%sf(0:m,0:n,0:p,1:nnode,1:nb))
            allocate (mv%sf(0:m,0:n,0:p,1:nnode,1:nb))
        end if

        do i = 1, sys_size
            ic%q_cons_vf(i)%sf = -1.e-6_stp  ! real(dflt_real, kind=stp) ! TODO :: remove this magic number
            ic%q_prim_vf(i)%sf = -1.e-6_stp  ! real(dflt_real, kind=stp)
        end do

        allocate (ic%bc_type(1:num_dims,1:2))

        allocate (ic%bc_type(1, 1)%sf(0:0,0:n,0:p))
        allocate (ic%bc_type(1, 2)%sf(0:0,0:n,0:p))

        do l = 0, p
            do k = 0, n
                ic%bc_type(1, 1)%sf(0, k, l) = int(min(bc_x%beg, 0), kind=1)
                ic%bc_type(1, 2)%sf(0, k, l) = int(min(bc_x%end, 0), kind=1)
            end do
        end do

        if (n > 0) then
            allocate (ic%bc_type(2, 1)%sf(-buff_size:m + buff_size,0:0,0:p))
            allocate (ic%bc_type(2, 2)%sf(-buff_size:m + buff_size,0:0,0:p))

            do l = 0, p
                do j = -buff_size, m + buff_size
                    ic%bc_type(2, 1)%sf(j, 0, l) = int(min(bc_y%beg, 0), kind=1)
                    ic%bc_type(2, 2)%sf(j, 0, l) = int(min(bc_y%end, 0), kind=1)
                end do
            end do

            if (p > 0) then
                allocate (ic%bc_type(3, 1)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (ic%bc_type(3, 2)%sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))

                do k = -buff_size, n + buff_size
                    do j = -buff_size, m + buff_size
                        ic%bc_type(3, 1)%sf(j, k, 0) = int(min(bc_z%beg, 0), kind=1)
                        ic%bc_type(3, 2)%sf(j, k, 0) = int(min(bc_z%end, 0), kind=1)
                    end do
                end do
            end if
        end if

        ! Initial damage state is always zero
        if (cont_damage) then
            ic%q_cons_vf(eqn_idx%damage)%sf = 0._wp
            ic%q_prim_vf(eqn_idx%damage)%sf = 0._wp
        end if

        ! Initial hyper_cleaning state is always zero TODO more general
        if (hyper_cleaning) then
            ic%q_cons_vf(eqn_idx%psi)%sf = 0._wp
            ic%q_prim_vf(eqn_idx%psi)%sf = 0._wp
        end if

        ! Setting default values for patch identities bookkeeping variable. This is necessary to avoid any confusion in the
        ! assessment of the extent of application that the overwrite permissions give a patch when it is being applied in the
        ! domain.
        ic%patch_id_fp = 0

    end subroutine s_initialize_initial_condition_module

    !> Iterate over patches and, depending on the geometry type, call the related subroutine to setup the said geometry on the grid
    !! using the primitive variables included with the patch parameters. The subroutine is complete once the primitive variables are
    !! converted to conservative ones.
    impure subroutine s_generate_initial_condition

        integer :: i

        if (old_ic) then
            call s_convert_conservative_to_primitive_variables(ic%q_cons_vf, ic%q_T_sf, ic%q_prim_vf, idwbuff)
        end if

        call s_apply_icpp_patches(ic%patch_id_fp, ic%q_prim_vf)

        if (num_bc_patches > 0) call s_apply_boundary_patches(ic%q_prim_vf, ic%bc_type)

        if (perturb_flow) call s_perturb_surrounding_flow(ic%q_prim_vf)
        if (perturb_sph) call s_perturb_sphere(ic%q_prim_vf)
        if (mixlayer_perturb) call s_perturb_mixlayer(ic%q_prim_vf)
        if (simplex_perturb) call s_perturb_simplex(ic%q_prim_vf)
        if (chemistry) call s_compute_T_from_primitives(ic%q_T_sf, ic%q_prim_vf, idwint)

        if (elliptic_smoothing .and. chemistry) then
            call s_elliptic_smoothing(ic%q_prim_vf, ic%bc_type, ic%q_T_sf)
            call s_compute_T_from_primitives(ic%q_T_sf, ic%q_prim_vf, idwint)
        else if (elliptic_smoothing) then
            call s_elliptic_smoothing(ic%q_prim_vf, ic%bc_type)
        end if

        call s_convert_primitive_to_conservative_variables(ic%q_prim_vf, ic%q_cons_vf)

        if (qbmm .and. .not. polytropic) then
            call s_initialize_mv(ic%q_cons_vf, mv%sf)
            call s_initialize_pb(ic%q_cons_vf, mv%sf, pb%sf)
        end if

    end subroutine s_generate_initial_condition

    !> AMR: generate the fine block's initial condition at FINE resolution and write it as lustre_amr_0, so the simulation
    !! initializes the fine level from a sharp interface instead of prolonging the coarse-width IC (which pins the fine interface to
    !! the coarse cell width). Swaps the global grid to the fine block, re-stamps the icpp patches at fine dx, converts to
    !! conservative, writes the AMR-restart record, and restores the base grid. P0 scope: parallel_io, single rank, one static
    !! un-tiled cartesian block; patches only.
    impure subroutine s_generate_amr_fine_ic

        integer                         :: bg(3), en(3), rr, fm, fn, fp, i, un, ierr, sm, sn, sp
        real(wp)                        :: sdx, sdy, sdz, xL, yL, zL
        real(wp), allocatable           :: sx_cc(:), sy_cc(:), sz_cc(:), sx_cb(:), sy_cb(:), sz_cb(:)
        type(int_bounds_info)           :: sidw(3)
        type(scalar_field), allocatable :: fq_prim(:), fq_cons(:)

#ifdef MFC_MIXED_PRECISION
        integer(kind=1), allocatable :: fpid(:,:,:)
#else
        integer, allocatable :: fpid(:,:,:)
#endif
        character(LEN=path_len + name_len) :: file_loc

        if (.not. amr) return
        if (.not. parallel_io .or. num_procs > 1 .or. cyl_coord) return  ! P0: parallel_io, single rank, cartesian
        bg = 0; en = 0
        bg(1) = amr_block_beg(1); en(1) = amr_block_end(1)
        if (n > 0) then
            bg(2) = amr_block_beg(2); en(2) = amr_block_end(2)
        end if
        if (p > 0) then
            bg(3) = amr_block_beg(3); en(3) = amr_block_end(3)
        end if
        if (en(1) < bg(1)) return  ! no static block defined (dynamic-only regrid)

        rr = amr_ref_ratio
        fm = rr*(en(1) - bg(1) + 1) - 1
        fn = 0; if (n > 0) fn = rr*(en(2) - bg(2) + 1) - 1
        fp = 0; if (p > 0) fp = rr*(en(3) - bg(3) + 1) - 1

        ! save the base grid
        sm = m; sn = n; sp = p; sdx = dx; sdy = dy; sdz = dz; sidw = idwbuff
        allocate (sx_cc, source=x_cc); allocate (sx_cb, source=x_cb)
        if (n > 0) then
            allocate (sy_cc, source=y_cc); allocate (sy_cb, source=y_cb)
        end if
        if (p > 0) then
            allocate (sz_cc, source=z_cc); allocate (sz_cb, source=z_cb)
        end if

        ! swap the global grid to a uniform subdivision of the block's level-0 cells
        m = fm; n = fn; p = fp
        dx = sdx/rr; if (n > 0) dy = sdy/rr; if (p > 0) dz = sdz/rr
        idwbuff(1)%beg = 0; idwbuff(1)%end = fm
        idwbuff(2)%beg = 0; idwbuff(2)%end = fn
        idwbuff(3)%beg = 0; idwbuff(3)%end = fp
        deallocate (x_cc, x_cb); allocate (x_cc(0:fm), x_cb(-1:fm))
        xL = sx_cb(bg(1) - 1)
        do i = -1, fm
            x_cb(i) = xL + real(i + 1, wp)*dx
        end do
        do i = 0, fm
            x_cc(i) = 0.5_wp*(x_cb(i - 1) + x_cb(i))
        end do
        if (n > 0) then
            deallocate (y_cc, y_cb); allocate (y_cc(0:fn), y_cb(-1:fn))
            yL = sy_cb(bg(2) - 1)
            do i = -1, fn
                y_cb(i) = yL + real(i + 1, wp)*dy
            end do
            do i = 0, fn
                y_cc(i) = 0.5_wp*(y_cb(i - 1) + y_cb(i))
            end do
        end if
        if (p > 0) then
            deallocate (z_cc, z_cb); allocate (z_cc(0:fp), z_cb(-1:fp))
            zL = sz_cb(bg(3) - 1)
            do i = -1, fp
                z_cb(i) = zL + real(i + 1, wp)*dz
            end do
            do i = 0, fp
                z_cc(i) = 0.5_wp*(z_cb(i - 1) + z_cb(i))
            end do
        end if

        ! stamp the icpp patches on the fine grid, then convert to conservative variables
        allocate (fq_prim(sys_size), fq_cons(sys_size))
        do i = 1, sys_size
            allocate (fq_prim(i)%sf(0:fm,0:fn,0:fp)); allocate (fq_cons(i)%sf(0:fm,0:fn,0:fp))
            fq_prim(i)%sf = -1.e-6_stp  ! background sentinel, exactly as s_initialize_initial_condition_module seeds the base
        end do  ! grid before patching - without it, cells no patch fully overwrites keep allocation garbage
        allocate (fpid(0:fm,0:fn,0:fp)); fpid = 0
        call s_apply_icpp_patches(fpid, fq_prim)
        call s_convert_primitive_to_conservative_variables(fq_prim, fq_cons)

        ! write lustre_amr_0 - same record layout as s_write_amr_restart (np=1 => contiguous stream bytes)
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // 'amr_0.dat'
        open (newunit=un, FILE=trim(file_loc), form='unformatted', access='stream', status='replace', IOSTAT=ierr)
        if (ierr == 0) then
            write (un) 1, 1, sys_size  ! num_procs, num_blocks, sys_size
            write (un) bg(1), bg(2), bg(3), en(1), en(2), en(3)  ! block region lo, hi (level-0 cell indices)
            write (un) fm, fn, fp  ! this rank's fine extents
            do i = 1, sys_size
                write (un) fq_cons(i)%sf(0:fm,0:fn,0:fp)
            end do
            close (un)
            if (proc_rank == 0) print '(A,I0,A,I0,A,I0,A)', ' [amr] wrote fine-resolution IC: level-0 ', en(1) - bg(1) + 1, &
                & '-cell block -> ', fm + 1, 'x', fn + 1, ' fine cells (lustre_amr_0)'
        end if

        ! restore the base grid
        m = sm; n = sn; p = sp; dx = sdx; dy = sdy; dz = sdz; idwbuff = sidw
        deallocate (x_cc, x_cb); call move_alloc(sx_cc, x_cc); call move_alloc(sx_cb, x_cb)
        if (n > 0) then
            deallocate (y_cc, y_cb); call move_alloc(sy_cc, y_cc); call move_alloc(sy_cb, y_cb)
        end if
        if (p > 0) then
            deallocate (z_cc, z_cb); call move_alloc(sz_cc, z_cc); call move_alloc(sz_cb, z_cb)
        end if
        do i = 1, sys_size
            deallocate (fq_prim(i)%sf, fq_cons(i)%sf)
        end do
        deallocate (fq_prim, fq_cons, fpid)

    end subroutine s_generate_amr_fine_ic

    !> Deallocation procedures for the module
    impure subroutine s_finalize_initial_condition_module

        integer :: i

        do i = 1, sys_size
            deallocate (ic%q_prim_vf(i)%sf)
            deallocate (ic%q_cons_vf(i)%sf)
        end do

        deallocate (ic%q_prim_vf)
        deallocate (ic%q_cons_vf)

        if (chemistry) then
            deallocate (ic%q_T_sf%sf)
        end if

        deallocate (ic%patch_id_fp)

        deallocate (ic%bc_type(1, 1)%sf)
        deallocate (ic%bc_type(1, 2)%sf)

        if (n > 0) then
            deallocate (ic%bc_type(2, 1)%sf)
            deallocate (ic%bc_type(2, 2)%sf)
        end if

        if (p > 0) then
            deallocate (ic%bc_type(3, 1)%sf)
            deallocate (ic%bc_type(3, 2)%sf)
        end if

        deallocate (ic%bc_type)

    end subroutine s_finalize_initial_condition_module

end module m_initial_condition
