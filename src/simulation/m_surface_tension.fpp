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
                        ! Floor sigma at a small positive value: a strongly negative sigma(Gamma) or
                        ! sigma(T) closure (e.g. a saturating surfactant) must not drive sigma <= 0,
                        ! which crashes the capillary force. Inert in normal operation (c_sigma >> floor).
                        c_sigma(j, k, l) = max(c_sigma(j, k, l), 1.e-3_wp*sigma)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end if

    end subroutine s_get_capillary

    !> Accumulate the Jain (2024) interface-confined surfactant-diffusion face flux into the surfactant slot of flux_src_vf for
    !! sweep direction idir. Flux = D_s*( d(Gamma_tilde)/dx_idir - 2*(0.5 - c)*n_idir*Gamma_tilde/eps ): isotropic diffusion of the
    !! smeared density Gamma_tilde plus an interfacial SHARPENING flux (second term, coefficient a = 2) that re-confines it to the
    !! interface from both sides, reproducing surface diffusion at the exact Laplace-Beltrami rate without the fragile
    !! Gamma_tilde/|grad c| division. c = q_prim(eqn_idx%c) is the color function, n = grad(c)/|grad c| the interface unit normal,
    !! eps ~= dx the interface thickness. The divergence of this flux conserves total surfactant. Flux at index x is the face
    !! between cells x and x+1 (thermal-conduction convention). Ref: Jain, JCP 515 (2024) 113277, Eq. (6).
    subroutine s_compute_surfactant_diffusion_flux(idir, q_prim_vf, flux_src_vf, irx, iry, irz)

        integer, intent(in)                                    :: idir
        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: flux_src_vf
        type(int_bounds_info), intent(in)                      :: irx, iry, irz
        real(wp)                                               :: normc_f, n_idir, c_face, gtil_face, dgtil, eps, flux_idir
        integer                                                :: x, y, z
        integer, dimension(3)                                  :: off

        is1 = irx; is2 = iry; is3 = irz
        $:GPU_UPDATE(device='[is1, is2, is3]')

        off = 0
        off(idir) = 1

        $:GPU_PARALLEL_LOOP(collapse=3, private='[x, y, z, normc_f, n_idir, c_face, gtil_face, dgtil, eps, flux_idir]', &
                            & copyin='[off]')
        do z = is3%beg, is3%end
            do y = is2%beg, is2%end
                do x = is1%beg, is1%end
                    normc_f = 0.5_wp*(c_divs(num_dims + 1)%sf(x, y, z) + c_divs(num_dims + 1)%sf(x + off(1), y + off(2), &
                                      & z + off(3)))

                    if (normc_f > capillary_cutoff) then
                        ! idir component of the interface unit normal n = grad(c)/|grad c|, face-averaged
                        n_idir = 0.5_wp*(c_divs(idir)%sf(x, y, z) + c_divs(idir)%sf(x + off(1), y + off(2), z + off(3)))/normc_f
                        ! face-averaged color function and surfactant density; compact idir gradient of Gamma_tilde
                        c_face = 0.5_wp*(q_prim_vf(eqn_idx%c)%sf(x, y, z) + q_prim_vf(eqn_idx%c)%sf(x + off(1), y + off(2), &
                                         & z + off(3)))
                        gtil_face = 0.5_wp*(q_prim_vf(eqn_idx%surf)%sf(x, y, z) + q_prim_vf(eqn_idx%surf)%sf(x + off(1), &
                                            & y + off(2), z + off(3)))
                        select case (idir)
                        case (1)
                            eps = x_cc(x + 1) - x_cc(x)
                            dgtil = (q_prim_vf(eqn_idx%surf)%sf(x + 1, y, z) - q_prim_vf(eqn_idx%surf)%sf(x, y, z))/eps
                        case (2)
                            eps = y_cc(y + 1) - y_cc(y)
                            dgtil = (q_prim_vf(eqn_idx%surf)%sf(x, y + 1, z) - q_prim_vf(eqn_idx%surf)%sf(x, y, z))/eps
                        case (3)
                            eps = z_cc(z + 1) - z_cc(z)
                            dgtil = (q_prim_vf(eqn_idx%surf)%sf(x, y, z + 1) - q_prim_vf(eqn_idx%surf)%sf(x, y, z))/eps
                        end select

                        ! Jain Eq. (6): isotropic diffusion minus interfacial sharpening flux (sharpening coefficient a = 2)
                        flux_idir = surf_diff*(dgtil - 2._wp*(0.5_wp - c_face)*n_idir*gtil_face/eps)

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

    end subroutine s_finalize_surface_tension_module

end module m_surface_tension
