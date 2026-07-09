!>
!! @file
!! @brief Contains module m_surface_tension

#:include 'case.fpp'
#:include 'macros.fpp'
#:include 'inline_capillary.fpp'

!> @brief Computes capillary source fluxes and color-function gradients for the diffuse-interface surface tension model
module m_surface_tension

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_weno
    use m_muscl
    use m_helper
    use m_boundary_common
    use m_sim_helpers

    implicit none

    private; public :: s_initialize_surface_tension_module, s_compute_capillary_source_flux, s_get_capillary, &
        & s_compute_surfactant_diffusion_flux, s_finalize_surface_tension_module

    !> @name color function gradient components and magnitude
    !> @{
    type(scalar_field), allocatable, dimension(:) :: c_divs
    !> @}
    $:GPU_DECLARE(create='[c_divs]')

    !> @name cell boundary reconstructed gradient components and magnitude
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: gL_x, gR_x
    !> @}
    $:GPU_DECLARE(create='[gL_x, gR_x]')

    !> @name cell-centered surface tension for the linear thermal closure sigma(T)
    !> (always allocated when surface tension is active; only written/read for sigma_model == 1,
    !> so the constant-sigma path never dereferences it but the device mapping is always valid)
    !> @{
    real(wp), allocatable, dimension(:,:,:) :: c_sigma
    !> @}
    $:GPU_DECLARE(create='[c_sigma]')

    !> @name cell-centered gradient of the surfactant density Gamma*|grad c|, used to build the
    !> tangential surface-diffusion flux (allocated unconditionally when surface tension is active
    !> so the device descriptor is valid; only written/read when surfactant .and. surf_diff > 0)
    !> @{
    real(wp), allocatable, dimension(:,:,:,:) :: dGamma_tilde
    !> @}
    $:GPU_DECLARE(create='[dGamma_tilde]')

    type(int_bounds_info) :: is1, is2, is3, iv
    $:GPU_DECLARE(create='[is1, is2, is3, iv]')

contains

    !> Allocate and initialize surface tension module arrays
    impure subroutine s_initialize_surface_tension_module

        integer :: j

        @:ALLOCATE(c_divs(1:num_dims + 1))

        do j = 1, num_dims + 1
            @:ALLOCATE(c_divs(j)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))
            @:ACC_SETUP_SFs(c_divs(j))
        end do

        @:ALLOCATE(gL_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))
        @:ALLOCATE(gR_x(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, num_dims + 1))

        ! Allocated unconditionally so the device descriptor is always valid in the capillary
        ! source-flux kernel; it is only written (s_get_capillary) and read when sigma_model == 1.
        @:ALLOCATE(c_sigma(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

        ! Surfactant-density gradient work array (tangential surface diffusion); allocated
        ! unconditionally for a valid device descriptor, written/read only for surf_diff > 0.
        @:ALLOCATE(dGamma_tilde(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end, &
                   & 1:num_dims))

    end subroutine s_initialize_surface_tension_module

    !> Compute the capillary source flux from reconstructed color-gradient fields
    subroutine s_compute_capillary_source_flux(vSrc_rsx_vf, flux_src_vf, id, isx, isy, isz)

        real(wp), dimension(-1:,-1:,-1:,1:), intent(in)        :: vSrc_rsx_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        integer, intent(in)                                    :: id
        type(int_bounds_info), intent(in)                      :: isx, isy, isz

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3, 3) :: Omega
        #:else
            real(wp), dimension(num_dims, num_dims) :: Omega
        #:endif
        real(wp) :: w1L, w1R, w2L, w2R, w3L, w3R, w1, w2, w3
        real(wp) :: normWL, normWR, normW
        real(wp) :: sigma_face
        integer  :: j, k, l, i

        if (id == 1) then
            $:GPU_PARALLEL_LOOP(collapse=3, &
                                & private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, sigma_face]')
            do l = isz%beg, isz%end
                do k = isy%beg, isy%end
                    do j = isx%beg, isx%end
                        w1L = gL_x(j, k, l, 1)
                        w2L = gL_x(j, k, l, 2)
                        w3L = 0._wp
                        if (p > 0) w3L = gL_x(j, k, l, 3)

                        w1R = gR_x(j + 1, k, l, 1)
                        w2R = gR_x(j + 1, k, l, 2)
                        w3R = 0._wp
                        if (p > 0) w3R = gR_x(j + 1, k, l, 3)

                        normWL = gL_x(j, k, l, num_dims + 1)
                        normWR = gR_x(j + 1, k, l, num_dims + 1)

                        w1 = (w1L + w1R)/2._wp
                        w2 = (w2L + w2R)/2._wp
                        w3 = (w3L + w3R)/2._wp
                        normW = (normWL + normWR)/2._wp

                        if (normW > capillary_cutoff) then
                            sigma_face = sigma
                            if (sigma_model /= 0) sigma_face = (c_sigma(j, k, l) + c_sigma(j + 1, k, l))/2._wp

                            @:compute_capillary_stress_tensor(sigma_face)

                            do i = 1, num_dims
                                flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, &
                                            & l) + Omega(1, i)

                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) + Omega(1, &
                                            & i)*vSrc_rsx_vf(j, k, l, i)
                            end do

                            ! Continuum surface force capillary stress, Schmidmayer et al. JCP (2017)
                            flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                        & l) + sigma_face*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 1)
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        else if (id == 2) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 1
                $:GPU_PARALLEL_LOOP(collapse=3, &
                                    & private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, sigma_face]')
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end
                            w1L = gL_x(j, k, l, 1)
                            w2L = gL_x(j, k, l, 2)
                            w3L = 0._wp
                            if (p > 0) w3L = gL_x(j, k, l, 3)

                            w1R = gR_x(j, k + 1, l, 1)
                            w2R = gR_x(j, k + 1, l, 2)
                            w3R = 0._wp
                            if (p > 0) w3R = gR_x(j, k + 1, l, 3)

                            normWL = gL_x(j, k, l, num_dims + 1)
                            normWR = gR_x(j, k + 1, l, num_dims + 1)

                            w1 = (w1L + w1R)/2._wp
                            w2 = (w2L + w2R)/2._wp
                            w3 = (w3L + w3R)/2._wp
                            normW = (normWL + normWR)/2._wp

                            if (normW > capillary_cutoff) then
                                sigma_face = sigma
                                if (sigma_model /= 0) sigma_face = (c_sigma(j, k, l) + c_sigma(j, k + 1, l))/2._wp

                                @:compute_capillary_stress_tensor(sigma_face)

                                do i = 1, num_dims
                                    flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, &
                                                & k, l) + Omega(2, i)

                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) + Omega(2, &
                                                & i)*vSrc_rsx_vf(j, k, l, i)
                                end do

                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                            & l) + sigma_face*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 2)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endif
        else if (id == 3) then
            #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
                $:GPU_PARALLEL_LOOP(collapse=3, &
                                    & private='[Omega, w1L, w2L, w3L, w1R, w2R, w3R, w1, w2, w3, normWL, normWR, normW, sigma_face]')
                do l = isz%beg, isz%end
                    do k = isy%beg, isy%end
                        do j = isx%beg, isx%end
                            w1L = gL_x(j, k, l, 1)
                            w2L = gL_x(j, k, l, 2)
                            w3L = 0._wp
                            if (p > 0) w3L = gL_x(j, k, l, 3)

                            w1R = gR_x(j, k, l + 1, 1)
                            w2R = gR_x(j, k, l + 1, 2)
                            w3R = 0._wp
                            if (p > 0) w3R = gR_x(j, k, l + 1, 3)

                            normWL = gL_x(j, k, l, num_dims + 1)
                            normWR = gR_x(j, k, l + 1, num_dims + 1)

                            w1 = (w1L + w1R)/2._wp
                            w2 = (w2L + w2R)/2._wp
                            w3 = (w3L + w3R)/2._wp
                            normW = (normWL + normWR)/2._wp

                            if (normW > capillary_cutoff) then
                                sigma_face = sigma
                                if (sigma_model /= 0) sigma_face = (c_sigma(j, k, l) + c_sigma(j, k, l + 1))/2._wp

                                @:compute_capillary_stress_tensor(sigma_face)

                                do i = 1, num_dims
                                    flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l) = flux_src_vf(eqn_idx%mom%beg + i - 1)%sf(j, &
                                                & k, l) + Omega(3, i)

                                    flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, l) + Omega(3, &
                                                & i)*vSrc_rsx_vf(j, k, l, i)
                                end do

                                flux_src_vf(eqn_idx%E)%sf(j, k, l) = flux_src_vf(eqn_idx%E)%sf(j, k, &
                                            & l) + sigma_face*c_divs(num_dims + 1)%sf(j, k, l)*vSrc_rsx_vf(j, k, l, 3)
                            end if
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            #:endif
        end if

    end subroutine s_compute_capillary_source_flux

    !> Compute color-function gradients and reconstruct them at cell boundaries
    impure subroutine s_get_capillary(q_prim_vf, bc_type)

        type(scalar_field), dimension(sys_size), intent(in)        :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(int_bounds_info)                                      :: isx, isy, isz
        integer                                                    :: j, k, l, i
        real(wp)                                                   :: T_cell, normc, Gamma_surf

        isx%beg = -1; isy%beg = 0; isz%beg = 0

        if (m > 0) isy%beg = -1; if (p > 0) isz%beg = -1

        isx%end = m; isy%end = n; isz%end = p

        ! compute gradient components
        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(1)%sf(j, k, l) = 1._wp/(x_cc(j + 1) - x_cc(j - 1))*(q_prim_vf(eqn_idx%c)%sf(j + 1, k, &
                           & l) - q_prim_vf(eqn_idx%c)%sf(j - 1, k, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(2)%sf(j, k, l) = 1._wp/(y_cc(k + 1) - y_cc(k - 1))*(q_prim_vf(eqn_idx%c)%sf(j, k + 1, &
                           & l) - q_prim_vf(eqn_idx%c)%sf(j, k - 1, l))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        if (p > 0) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        c_divs(3)%sf(j, k, l) = 1._wp/(z_cc(l + 1) - z_cc(l - 1))*(q_prim_vf(eqn_idx%c)%sf(j, k, &
                               & l + 1) - q_prim_vf(eqn_idx%c)%sf(j, k, l - 1))
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    c_divs(num_dims + 1)%sf(j, k, l) = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_dims
                        c_divs(num_dims + 1)%sf(j, k, l) = c_divs(num_dims + 1)%sf(j, k, l) + c_divs(i)%sf(j, k, l)**2._wp
                    end do

                    c_divs(num_dims + 1)%sf(j, k, l) = sqrt(real(c_divs(num_dims + 1)%sf(j, k, l), kind=wp))
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_populate_capillary_buffers(c_divs, bc_type, bc_xyz_info(bc_x, bc_y, bc_z))

        iv%beg = 1; iv%end = num_dims + 1

        ! reconstruct gradient components at cell boundaries
        call s_reconstruct_cell_boundary_values_capillary(c_divs, gL_x, gR_x, i)

        ! Variable surface-tension closures fill the cell-centered field c_sigma over the full
        ! buffer range; contributions are additive so closures compose. The face value used by
        ! the capillary force is later formed by averaging adjacent cells. q_prim_vf ghost cells
        ! are already populated, so no extra halo exchange is needed.
        !   sigma_model == 1: linear thermal closure sigma(T) = sigma + dsigma/dT*(T - T_ref),
        !     T recovered from the mixture stiffened-gas EOS
        !     (T = ((gamma_mix + 1)*p + pi_inf_mix)/mCP, mCP = sum(alpha*rho*cv*gamma)).
        !   sigma_model == 2: linear solutocapillary closure sigma(Gamma) = sigma + dsigma/dGamma*Gamma,
        !     interfacial concentration Gamma = (Gamma*|grad c|)/|grad c| recovered on the interface
        !     band (|grad c| > capillary_cutoff); off-band cells keep the clean value sigma.
        if (sigma_model /= 0) then
            $:GPU_PARALLEL_LOOP(collapse=3, private='[T_cell, normc, Gamma_surf]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        c_sigma(j, k, l) = sigma
                        if (sigma_model == 1) then
                            T_cell = f_compute_mixture_temperature(q_prim_vf, j, k, l)
                            c_sigma(j, k, l) = c_sigma(j, k, l) + sigma_dTdT*(T_cell - sigma_T_ref)
                        end if
                        if (sigma_model == 2) then
                            normc = c_divs(num_dims + 1)%sf(j, k, l)
                            if (normc > capillary_cutoff) then
                                Gamma_surf = q_prim_vf(eqn_idx%surf)%sf(j, k, l)/normc
                                c_sigma(j, k, l) = c_sigma(j, k, l) + sigma_dGamma*Gamma_surf
                            end if
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        ! Cell-centered gradient of the surfactant density Gamma*|grad c| (= q_prim(surf)) over
        ! the interior plus one buffer layer; the diffusion flux producer averages these to faces
        ! for the tangential (non-idir) components. q_prim_vf ghost cells are already populated.
        if (surfactant .and. surf_diff > 0._wp) then
            $:GPU_PARALLEL_LOOP(collapse=3)
            do l = idwbuff(3)%beg + 1, idwbuff(3)%end - 1
                do k = idwbuff(2)%beg + 1, idwbuff(2)%end - 1
                    do j = idwbuff(1)%beg + 1, idwbuff(1)%end - 1
                        dGamma_tilde(j, k, l, 1) = (q_prim_vf(eqn_idx%surf)%sf(j + 1, k, l) - q_prim_vf(eqn_idx%surf)%sf(j - 1, &
                                     & k, l))/(x_cc(j + 1) - x_cc(j - 1))
                        dGamma_tilde(j, k, l, 2) = (q_prim_vf(eqn_idx%surf)%sf(j, k + 1, l) - q_prim_vf(eqn_idx%surf)%sf(j, &
                                     & k - 1, l))/(y_cc(k + 1) - y_cc(k - 1))
                        if (p > 0) then
                            dGamma_tilde(j, k, l, 3) = (q_prim_vf(eqn_idx%surf)%sf(j, k, l + 1) - q_prim_vf(eqn_idx%surf)%sf(j, &
                                         & k, l - 1))/(z_cc(l + 1) - z_cc(l - 1))
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_get_capillary

    !> Accumulate the tangential surfactant-diffusion face flux -D_s*(I - n(x)n)*grad(Gamma_tilde) into the surfactant slot of
    !! flux_src_vf for sweep direction idir. The projection (I - n(x)n) removes the across-interface (normal) gradient so surfactant
    !! diffuses only along the interface, and the divergence of the resulting flux conserves total surfactant. The idir gradient
    !! component uses the compact two-point difference; the tangential components are the cell-centered gradients averaged to the
    !! face. The flux at index x is the face between cells x and x+1 (thermal-conduction convention).
    subroutine s_compute_surfactant_diffusion_flux(idir, q_prim_vf, flux_src_vf, irx, iry, irz)

        integer, intent(in) :: idir
        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in) :: irx, iry, irz
        real(wp) :: normc_f, nrm1, nrm2, nrm3, gt1, gt2, gt3, ndotg, grid_spacing, flux_idir
        integer :: x, y, z
        integer, dimension(3) :: off

        is1 = irx; is2 = iry; is3 = irz
        $:GPU_UPDATE(device='[is1, is2, is3]')

        off = 0
        off(idir) = 1

        $:GPU_PARALLEL_LOOP(collapse=3, private='[x, y, z, normc_f, nrm1, nrm2, nrm3, gt1, gt2, gt3, ndotg, grid_spacing, &
                            & flux_idir]', copyin='[off]')
        do z = is3%beg, is3%end
            do y = is2%beg, is2%end
                do x = is1%beg, is1%end
                    normc_f = 0.5_wp*(c_divs(num_dims + 1)%sf(x, y, z) + c_divs(num_dims + 1)%sf(x + off(1), y + off(2), &
                                      & z + off(3)))

                    if (normc_f > capillary_cutoff) then
                        ! Face-averaged interface unit normal n = grad(c)/|grad c|
                        nrm1 = 0.5_wp*(c_divs(1)%sf(x, y, z) + c_divs(1)%sf(x + off(1), y + off(2), z + off(3)))/normc_f
                        nrm2 = 0.5_wp*(c_divs(2)%sf(x, y, z) + c_divs(2)%sf(x + off(1), y + off(2), z + off(3)))/normc_f
                        nrm3 = 0._wp
                        if (p > 0) nrm3 = 0.5_wp*(c_divs(3)%sf(x, y, z) + c_divs(3)%sf(x + off(1), y + off(2), z + off(3)))/normc_f

                        ! Face gradient of Gamma_tilde: tangential components averaged from cell gradients
                        gt1 = 0.5_wp*(dGamma_tilde(x, y, z, 1) + dGamma_tilde(x + off(1), y + off(2), z + off(3), 1))
                        gt2 = 0.5_wp*(dGamma_tilde(x, y, z, 2) + dGamma_tilde(x + off(1), y + off(2), z + off(3), 2))
                        gt3 = 0._wp
                        if (p > 0) gt3 = 0.5_wp*(dGamma_tilde(x, y, z, 3) + dGamma_tilde(x + off(1), y + off(2), z + off(3), 3))

                        ! Sharper compact two-point difference for the idir component
                        select case (idir)
                        case (1)
                            grid_spacing = x_cc(x + 1) - x_cc(x)
                            gt1 = (q_prim_vf(eqn_idx%surf)%sf(x + 1, y, z) - q_prim_vf(eqn_idx%surf)%sf(x, y, z))/grid_spacing
                        case (2)
                            grid_spacing = y_cc(y + 1) - y_cc(y)
                            gt2 = (q_prim_vf(eqn_idx%surf)%sf(x, y + 1, z) - q_prim_vf(eqn_idx%surf)%sf(x, y, z))/grid_spacing
                        case (3)
                            grid_spacing = z_cc(z + 1) - z_cc(z)
                            gt3 = (q_prim_vf(eqn_idx%surf)%sf(x, y, z + 1) - q_prim_vf(eqn_idx%surf)%sf(x, y, z))/grid_spacing
                        end select

                        ! Tangential projection: subtract the normal component of grad(Gamma_tilde)
                        ndotg = nrm1*gt1 + nrm2*gt2 + nrm3*gt3
                        select case (idir)
                        case (1)
                            flux_idir = surf_diff*(gt1 - nrm1*ndotg)
                        case (2)
                            flux_idir = surf_diff*(gt2 - nrm2*ndotg)
                        case (3)
                            flux_idir = surf_diff*(gt3 - nrm3*ndotg)
                        end select

                        flux_src_vf(eqn_idx%surf)%sf(x, y, z) = flux_src_vf(eqn_idx%surf)%sf(x, y, z) - flux_idir
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_surfactant_diffusion_flux

    !> Reconstruct left and right cell-boundary values of capillary variables
    subroutine s_reconstruct_cell_boundary_values_capillary(v_vf, vL_x, vR_x, norm_dir)

        type(scalar_field), dimension(iv%beg:iv%end), intent(in)                                  :: v_vf
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,iv%beg:), intent(out) :: vL_x
        real(wp), dimension(idwbuff(1)%beg:,idwbuff(2)%beg:,idwbuff(3)%beg:,iv%beg:), intent(out) :: vR_x
        integer, intent(in)                                                                       :: norm_dir
        integer                                                                                   :: i, j, k, l

        $:GPU_UPDATE(device='[iv]')

        $:GPU_PARALLEL_LOOP(collapse=4, private='[i, j, k, l]')
        do i = iv%beg, iv%end
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        vL_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                        vR_x(j, k, l, i) = v_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_reconstruct_cell_boundary_values_capillary

    !> Finalize the surface tension module
    impure subroutine s_finalize_surface_tension_module

        integer :: j

        do j = 1, num_dims + 1
            @:DEALLOCATE(c_divs(j)%sf)
        end do

        @:DEALLOCATE(c_divs)

        @:DEALLOCATE(gL_x, gR_x)

        @:DEALLOCATE(c_sigma)

        @:DEALLOCATE(dGamma_tilde)

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension
