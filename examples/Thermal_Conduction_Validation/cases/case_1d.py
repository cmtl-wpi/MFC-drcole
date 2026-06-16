#!/usr/bin/env python3
# 1D heat equation, Dirichlet sine decay. Temperature is recovered from the
# stiffened-gas EOS: the sine profile is imposed through the density at uniform
# pressure, rho(x) = (p0+p_inf)/((gam-1)*cv*T(x)), T(x) = Twall + A*sin(pi*x/L);
# walls held at Twall.
import json
import math

L = 1.0
gam, p_inf, cv = 2.0, 100.0, 12.5
p0, Twall, A = 25.0, 10.0, 3.0
alpha = 0.05
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * Twall)  # = 1.0
k_therm = alpha * rho0 * cv * gam
kappa = math.pi / L

Nx = 255
dx = L / (Nx + 1)
c0 = (gam * (p0 + p_inf)) ** 0.5
dt = min(0.4 * dx / c0, 0.35 * dx**2 / (2.0 * alpha))
t_end = 1.5 / (alpha * kappa**2)
Nt = int(round(t_end / dt))

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "m": Nx,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(1, Nt // 40),
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
            "bc_x%isothermal_in": "T",
            "bc_x%isothermal_out": "T",
            "bc_x%Twall_in": Twall,
            "bc_x%Twall_out": Twall,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "T_wrt": "T",  # write temperature recovered from the EOS so post_process/viz can plot T
            "parallel_io": "T",
            # Thermal conduction (EOS temperature)
            "thermal_conduction": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
            "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
            "fluid_pp(1)%cv": cv,
            "fluid_pp(1)%k_therm": k_therm,
            # Patch 1: full domain, sine T imposed through density
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": L / 2,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": f"{rho0 * Twall:.12f}/({Twall} + {A}*sin({kappa:.12f}*x))",
            "patch_icpp(1)%alpha(1)": 1.0,
        }
    )
)
