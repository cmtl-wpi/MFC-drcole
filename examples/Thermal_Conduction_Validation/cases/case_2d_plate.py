#!/usr/bin/env python3
# 2D heat equation, steady-state plate. Three edges held at T=0, top edge at
# T=100; marched to steady state (Laplace). thermal_scalar -> u stays 0.
import json
import math

Lx = Ly = 1.0
gam, p_inf, cv = 2.0, 100.0, 12.5
p0, Tcold, Thot = 25.0, 10.0, 110.0  # cold edges 10, top 110 -> dT=100 (MFC requires Twall>0)
alpha = 0.05
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * Tcold)
k_therm = alpha * rho0 * cv * gam

Nx = Ny = 127
dx = Lx / (Nx + 1)
c0 = (gam * (p0 + p_inf)) ** 0.5
dt = min(0.4 * dx / c0, 0.35 * dx**2 / (4.0 * alpha))
# slowest transient decays at alpha*((pi/Lx)^2 + (pi/Ly)^2); run 8 e-folds to steady
t_end = 8.0 / (alpha * ((math.pi / Lx) ** 2 + (math.pi / Ly) ** 2))
Nt = int(round(t_end / dt))

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": Lx,
            "y_domain%beg": 0.0,
            "y_domain%end": Ly,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(1, Nt // 30),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            "bc_y%beg": -3,
            "bc_y%end": -3,
            "bc_x%isothermal_in": "T",
            "bc_x%isothermal_out": "T",
            "bc_y%isothermal_in": "T",
            "bc_y%isothermal_out": "T",
            "bc_x%Twall_in": Tcold,
            "bc_x%Twall_out": Tcold,
            "bc_y%Twall_in": Tcold,
            "bc_y%Twall_out": Thot,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "T",
            # Thermal conduction
            "thermal_conduction": "T",
            "thermal_scalar": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
            "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
            "fluid_pp(1)%cv": cv,
            "fluid_pp(1)%k_therm": k_therm,
            # Patch 1: full domain, cold start (T_s = 0)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": Lx / 2,
            "patch_icpp(1)%y_centroid": Ly / 2,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%T_temp_val": Tcold,  # cold start at the cold-edge temperature
        }
    )
)
