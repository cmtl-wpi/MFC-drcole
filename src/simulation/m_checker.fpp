!>
!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Validates simulation input parameters for consistency and supported configurations
module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_helper_basic

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers

        if (igr) then
            call s_check_inputs_nvidia_uvm
        else
            if (recon_type == WENO_TYPE) then
                call s_check_inputs_weno
            else if (recon_type == MUSCL_TYPE) then
                call s_check_inputs_muscl
            end if
        end if

        call s_check_inputs_time_stepping

        call s_check_inputs_thermal_conduction

        call s_check_inputs_thermal_scalar

        @:PROHIBIT(ib_state_wrt .and. .not. ib, "ib_state_wrt requires ib to be enabled")

    end subroutine s_check_inputs

    !> Checks constraints on bulk thermal conduction parameters
    impure subroutine s_check_inputs_thermal_conduction

        integer :: i

        if (thermal_conduction) then
            @:PROHIBIT(chemistry, "thermal_conduction is not supported with chemistry; use chem_params%diffusion instead")
            @:PROHIBIT(igr, "thermal_conduction is not supported with igr")
            @:PROHIBIT(cyl_coord, "thermal_conduction is not supported with cyl_coord")
            @:PROHIBIT(bubbles_euler .or. bubbles_lagrange, "thermal_conduction is not supported with bubble models")
            @:PROHIBIT(hypoelasticity .or. hyperelasticity, "thermal_conduction is not supported with elasticity")
            @:PROHIBIT(mhd, "thermal_conduction is not supported with mhd")
            @:PROHIBIT(relax, "thermal_conduction is not supported with relax (phase change)")
            @:PROHIBIT(ib, "thermal_conduction is not supported with immersed boundaries")
            @:PROHIBIT(model_eqns == 1 .or. model_eqns == 4, &
                       & "thermal_conduction requires model_eqns = 2 or 3 (mixture stiffened-gas temperature)")

            do i = 1, num_fluids
                @:PROHIBIT(fluid_pp(i)%cv <= 0._wp, &
                           & "thermal_conduction requires fluid_pp(i)%cv > 0 for every fluid to evaluate temperature")
                @:PROHIBIT(fluid_pp(i)%k_therm <= 0._wp, &
                           & "thermal_conduction requires fluid_pp(i)%k_therm > 0 for every fluid (harmonic mixture closure)")
            end do
        end if

    end subroutine s_check_inputs_thermal_conduction

    !> Checks constraints on the independent temperature scalar (thermal_scalar)
    impure subroutine s_check_inputs_thermal_scalar

        integer :: i

        if (thermal_scalar) then
            @:PROHIBIT(igr, "thermal_scalar is not supported with igr")
            @:PROHIBIT(cyl_coord, "thermal_scalar is not supported with cyl_coord")
            @:PROHIBIT(chemistry, "thermal_scalar is not supported with chemistry")
            @:PROHIBIT(bubbles_euler .or. bubbles_lagrange, "thermal_scalar is not supported with bubble models")
            @:PROHIBIT(hypoelasticity .or. hyperelasticity, "thermal_scalar is not supported with elasticity")
            @:PROHIBIT(mhd, "thermal_scalar is not supported with mhd")
            @:PROHIBIT(relax, "thermal_scalar is not supported with relax (phase change)")
            @:PROHIBIT(ib, "thermal_scalar is not supported with immersed boundaries")
            @:PROHIBIT(model_eqns == 1 .or. model_eqns == 4, "thermal_scalar requires model_eqns = 2 or 3")
            @:PROHIBIT(riemann_solver /= 1 .and. riemann_solver /= 2, &
                       & "thermal_scalar requires riemann_solver = 1 (HLL) or 2 (HLLC) for passive-scalar advection")

            if (thermal_conduction) then
                do i = 1, num_fluids
                    @:PROHIBIT(fluid_pp(i)%cv <= 0._wp, &
                               & "thermal_scalar with thermal_conduction requires fluid_pp(i)%cv > 0 for every fluid (rho*cp diffusivity)")
                end do
            end if
        end if

    end subroutine s_check_inputs_thermal_scalar

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

#if !defined(MFC_OpenACC) && !(defined(__PGI) || defined(_CRAYFTN))
        @:PROHIBIT(rdma_mpi, "Unsupported value of rdma_mpi for the current compiler")
#endif

    end subroutine s_check_inputs_compilers

    !> Checks constraints on WENO scheme parameters
    impure subroutine s_check_inputs_weno

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*weno_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*weno_order, &
                   & "m must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*weno_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*weno_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*weno_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_weno

    !> Validate that the grid resolution is sufficient for the MUSCL reconstruction order
    impure subroutine s_check_inputs_muscl

        character(len=5) :: numStr  !< for int to string conversion

        call s_int_to_str(num_stcls_min*muscl_order, numStr)
        @:PROHIBIT(m + 1 < num_stcls_min*muscl_order, &
                   & "m must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " // trim(numStr))
        @:PROHIBIT(n + 1 < min(1, n)*num_stcls_min*muscl_order, &
                   & "For 2D simulation, n must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))
        @:PROHIBIT(p + 1 < min(1, p)*num_stcls_min*muscl_order, &
                   & "For 3D simulation, p must be greater than or equal to (num_stcls_min*muscl_order - 1), whose value is " &
                   & // trim(numStr))

    end subroutine s_check_inputs_muscl

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

    !> Validate NVIDIA unified virtual memory configuration parameters
    impure subroutine s_check_inputs_nvidia_uvm

#ifdef __NVCOMPILER_GPU_UNIFIED_MEM
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu > 3 .or. nv_uvm_igr_temps_on_gpu < 0, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 3]")
        @:PROHIBIT(nv_uvm_igr_temps_on_gpu == 3 .and. igr_iter_solver == 2, &
                   & "nv_uvm_igr_temps_on_gpu must be in the range [0, 2] for igr_iter_solver == 2")
#endif

    end subroutine s_check_inputs_nvidia_uvm

end module m_checker
