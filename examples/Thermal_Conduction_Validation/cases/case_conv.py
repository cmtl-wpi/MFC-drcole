#!/usr/bin/env python3
# Convergence driver: 1D single Fourier mode on a periodic domain (no BC error).
# T_s = T0 + A*sin(2*pi*x/L). The harness sets grid and time step via the
# environment so it can hold the diffusion number fixed (spatial sweep) or vary
# dt on a fixed grid (temporal sweep). thermal_scalar -> u stays 0.
import json
import math
import os

L = 1.0
gam, p_inf, cv = 2.0, 100.0, 12.5
p0, T0, A = 25.0, 10.0, 3.0
alpha = 0.05
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * T0)
k_therm = alpha * rho0 * cv * gam
kw = 2.0 * math.pi / L

Nx = int(os.environ.get("CONV_N", "128")) - 1
dx = L / (Nx + 1)
dt = float(os.environ.get("CONV_DT", str(0.2 * dx**2 / alpha)))
Nt = int(os.environ.get("CONV_NSTEPS", str(int(round(0.3 / (alpha * kw**2) / dt)))))

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "F",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": Nt,  # write only initial + final
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
            "bc_x%beg": -1,
            "bc_x%end": -1,
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
            # Patch 1: full domain, periodic sine mode
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": L / 2,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            "patch_icpp(1)%T_temp_val": f"{T0} + {A}*sin({kw:.12f}*x)",
        }
    )
)
