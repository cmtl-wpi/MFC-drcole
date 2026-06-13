!>
!! @file
!! @brief Contains module m_thermal_conduction

#:include 'macros.fpp'

!> @brief Bulk Fourier heat conduction for non-reacting flows: an explicit -k*grad(T) face flux accumulated into the energy slot of
!! the source flux array. Temperature is recovered from the mixture stiffened-gas EOS and the cell conductivity follows the harmonic
!! mixture closure 1/k = sum(alpha_i/k_i) (Samareh et al. 2014, Eq. 8) over the per-fluid constants fluid_pp(i)%k_therm. The face
!! flux stored at index x is the face between cells x and x+1 (chemistry diffusion convention) and is differenced in
!! s_compute_additional_physics_rhs. Independent of the chemistry module's chem_params%diffusion path.
module m_thermal_conduction

    use m_derived_types      !< Definitions of the derived types
    use m_global_parameters  !< Definitions of the global parameters
    use m_sim_helpers        !< Mixture stiffened-gas temperature helper

    implicit none

    private; public :: s_initialize_thermal_conduction_module, s_get_thermal_conduction, s_compute_conductive_flux, &
        & s_finalize_thermal_conduction_module

    real(wp), allocatable, dimension(:,:,:) :: T_tc  !< Cell-centered mixture temperature
    $:GPU_DECLARE(create='[T_tc]')

    type(int_bounds_info) :: isc1_tc, isc2_tc, isc3_tc
    $:GPU_DECLARE(create='[isc1_tc, isc2_tc, isc3_tc]')

contains

    !> Allocate the temperature work array over the full buffer extent
    impure subroutine s_initialize_thermal_conduction_module

        @:ALLOCATE(T_tc(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, idwbuff(3)%beg:idwbuff(3)%end))

    end subroutine s_initialize_thermal_conduction_module

    !> Refresh the cell-centered temperature from primitive variables over the full buffer extent. q_prim_vf ghost cells are already
    !! populated when this is called, so no extra halo exchange is needed. Isothermal boundaries then overwrite the ghost
    !! temperature so that an imposed wall/far-field temperature can drive the conductive flux.
    subroutine s_get_thermal_conduction(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer                                             :: j, k, l

        $:GPU_PARALLEL_LOOP(collapse=3)
        do l = idwbuff(3)%beg, idwbuff(3)%end
            do k = idwbuff(2)%beg, idwbuff(2)%end
                do j = idwbuff(1)%beg, idwbuff(1)%end
                    T_tc(j, k, l) = f_compute_mixture_temperature(q_prim_vf, j, k, l)
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

        call s_apply_thermal_conduction_bc()

    end subroutine s_get_thermal_conduction

    !> Overwrite the ghost-cell temperature at isothermal boundaries with a Dirichlet reflection T_ghost = 2*Twall - T_interior, so
    !! the face temperature equals Twall and the conductive flux through the boundary is set by the prescribed wall/far-field
    !! temperature. Mirrors the chemistry q_T_sf isothermal handling in m_boundary_common, but is standalone: it applies for any
    !! boundary type (e.g. an open boundary holding a far-field temperature), not just no-slip/slip walls. Boundaries left
    !! non-isothermal keep the extrapolated (zero-gradient, adiabatic) ghost temperature.
    subroutine s_apply_thermal_conduction_bc()

        integer :: j, k, l

        if (bc_x%isothermal_in) then
            $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do j = 1, buff_size
                        T_tc(-j, k, l) = 2._wp*bc_x%Twall_in - T_tc(j - 1, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (bc_x%isothermal_out) then
            $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    $:GPU_LOOP(parallelism='[seq]')
                    do j = 1, buff_size
                        T_tc(m + j, k, l) = 2._wp*bc_x%Twall_out - T_tc(m - (j - 1), k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

        if (n > 0) then
            if (bc_y%isothermal_in) then
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do j = 1, buff_size
                            T_tc(k, -j, l) = 2._wp*bc_y%Twall_in - T_tc(k, j - 1, l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bc_y%isothermal_out) then
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                do l = idwbuff(3)%beg, idwbuff(3)%end
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do j = 1, buff_size
                            T_tc(k, n + j, l) = 2._wp*bc_y%Twall_out - T_tc(k, n - (j - 1), l)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

        if (p > 0) then
            if (bc_z%isothermal_in) then
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                do l = idwbuff(2)%beg, idwbuff(2)%end
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do j = 1, buff_size
                            T_tc(k, l, -j) = 2._wp*bc_z%Twall_in - T_tc(k, l, j - 1)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            if (bc_z%isothermal_out) then
                $:GPU_PARALLEL_LOOP(collapse=2, private='[j, k, l]')
                do l = idwbuff(2)%beg, idwbuff(2)%end
                    do k = idwbuff(1)%beg, idwbuff(1)%end
                        $:GPU_LOOP(parallelism='[seq]')
                        do j = 1, buff_size
                            T_tc(k, l, p + j) = 2._wp*bc_z%Twall_out - T_tc(k, l, p - (j - 1))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if
        end if

    end subroutine s_apply_thermal_conduction_bc

    !> Accumulate the conductive face flux -k_face*dT/dxi into the energy slot of flux_src_vf for sweep direction idir. The 2-point
    !! difference of cell-center temperatures is divided by the cell-center spacing, so stretched grids are handled correctly.
    subroutine s_compute_conductive_flux(idir, q_prim_vf, flux_src_vf, irx, iry, irz)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in)                      :: irx, iry, irz
        real(wp)                                               :: k_L, k_R, k_face, dT_dxi, grid_spacing
        integer                                                :: x, y, z, i
        integer, dimension(3)                                  :: offsets

        isc1_tc = irx; isc2_tc = iry; isc3_tc = irz

        $:GPU_UPDATE(device='[isc1_tc, isc2_tc, isc3_tc]')

        offsets = 0
        offsets(idir) = 1

        $:GPU_PARALLEL_LOOP(collapse=3, private='[x, y, z, i, k_L, k_R, k_face, dT_dxi, grid_spacing]', copyin='[offsets]')
        do z = isc3_tc%beg, isc3_tc%end
            do y = isc2_tc%beg, isc2_tc%end
                do x = isc1_tc%beg, isc1_tc%end
                    select case (idir)
                    case (1)
                        grid_spacing = x_cc(x + 1) - x_cc(x)
                    case (2)
                        grid_spacing = y_cc(y + 1) - y_cc(y)
                    case (3)
                        grid_spacing = z_cc(z + 1) - z_cc(z)
                    end select

                    ! Samareh et al. (2014) Eq. 8 harmonic mixture closure 1/k = sum(alpha_i/k_i)
                    ! (series resistance across the diffuse interface). Volume fractions can
                    ! slightly under/overshoot near interfaces, so clamp them to [0,1] and guard
                    ! the reciprocals with sgm_eps; the checker enforces k_therm > 0 per fluid.
                    k_L = 0._wp; k_R = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_fluids
                        k_L = k_L + min(max(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(x, y, z), 0._wp), 1._wp)/max(kappas(i), sgm_eps)
                        k_R = k_R + min(max(q_prim_vf(eqn_idx%adv%beg + i - 1)%sf(x + offsets(1), y + offsets(2), &
                                        & z + offsets(3)), 0._wp), 1._wp)/max(kappas(i), sgm_eps)
                    end do
                    k_L = 1._wp/max(k_L, sgm_eps)
                    k_R = 1._wp/max(k_R, sgm_eps)

                    ! Face value from the two cell values, matching the chemistry template
                    k_face = 0.5_wp*(k_L + k_R)

                    dT_dxi = (T_tc(x + offsets(1), y + offsets(2), z + offsets(3)) - T_tc(x, y, z))/grid_spacing

                    flux_src_vf(eqn_idx%E)%sf(x, y, z) = flux_src_vf(eqn_idx%E)%sf(x, y, z) - k_face*dT_dxi
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_compute_conductive_flux

    impure subroutine s_finalize_thermal_conduction_module

        @:DEALLOCATE(T_tc)

    end subroutine s_finalize_thermal_conduction_module

end module m_thermal_conduction
