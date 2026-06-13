#!/usr/bin/env python3
# Isolation test for the bulk-conduction flux + Dirichlet temperature BC, with NO droplet,
# NO surface tension, NO viscosity. Single uniform-composition fluid carrying the same
# density-gradient-imposed linear temperature T(x) = T0 + gradT*x at uniform pressure that the
# thermocapillary example uses. For a LINEAR T and uniform k, div(k grad T) = 0, so this state is
# an exact static equilibrium: a correct conduction implementation must keep velocity ~ 0 and the
# temperature profile linear. Any growing velocity isolates a bug in the flux or the BC, separate
# from the Marangoni / droplet coupling.
import json

m = 128
L = 7.5
gam, p_inf, cv, p0 = 2.0, 100.0, 12.5, 25.0
T0, gradT = 10.0, 2.0 / 15.0
rho_coeff = (p0 + p_inf) / ((gam - 1.0) * cv)  # = 10
# alpha_T from Ma = 0.3 (same as the example's 2D case); k = alpha_T * rho_ref * cv * gam
Ur = 0.1 * gradT * 0.5 / 0.1
alpha_T = Ur * 0.5 / 0.3
k_therm = alpha_T * (rho_coeff / T0) * cv * gam
dx = L / (m + 1)
c0 = (gam * (p0 + p_inf)) ** 0.5
dt = min(0.4 * dx / c0, 0.35 * dx**2 / (2.0 * 1 * alpha_T))
n_tau = 3.0
tau = (rho_coeff / T0) * 0.5**2 / 0.1
t_step_stop = int(round(n_tau * tau / dt))

rho_expr = f"{rho_coeff:.9f}/({T0} + {gradT:.9f}*x)"

data = {
    "run_time_info": "T",
    "x_domain%beg": -L / 2,
    "x_domain%end": L / 2,
    "m": m,
    "n": 0,
    "p": 0,
    "cyl_coord": "F",
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": max(1, t_step_stop // 30),
    "model_eqns": 2,
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "num_patches": 1,
    # bulk conduction with the Dirichlet far-field temperatures
    "thermal_conduction": "T",
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%k_therm": k_therm,
    "bc_x%isothermal_in": "T",
    "bc_x%isothermal_out": "T",
    "bc_x%Twall_in": T0 + gradT * (-L / 2),
    "bc_x%Twall_out": T0 + gradT * (L / 2),
    # single uniform-composition patch spanning the domain, u = 0, density carries the linear T
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha(1)": 1.0,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "parallel_io": "T",
}
print(json.dumps(data))
