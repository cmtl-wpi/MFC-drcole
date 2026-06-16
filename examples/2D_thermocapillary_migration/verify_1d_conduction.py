#!/usr/bin/env python3
# Isolation test for the bulk-conduction flux + Dirichlet temperature BC, with NO droplet,
# NO surface tension, NO viscosity. Single uniform-composition fluid carrying the same
# density-gradient-imposed linear temperature T(x) = T0 + gradT*x at uniform pressure that the
# thermocapillary example uses. For a LINEAR T and uniform k, div(k grad T) = 0, so this state is
# an exact static equilibrium: a correct conduction implementation must keep velocity ~ 0 and the
# temperature profile linear. Any growing velocity isolates a bug in the flux or the BC, separate
# from the Marangoni / droplet coupling.
import json

# -- Samareh (2014) parameters (same values as the 2D example) --
D = 1.0       # droplet diameter (reference length)
r = D / 2.0   # droplet radius = 0.5
mu = 0.1      # dynamic viscosity of both phases
dsigma_dT = -0.1  # surface-tension slope dsigma/dT < 0
sigma0 = 0.1      # surface tension at T_ref
Ma_th = 0.3       # thermal Marangoni number (matching the 2D case)

# -- EOS (single fluid, no droplet) --
# Chosen so rho_coeff / T0 = 1.0 at the reference density, keeping the
# compressible acoustic CFL manageable. The values differ from the 2D case
# (which tunes cv for rho_drop=0.2) because this test has a single fluid.
# T = (p + p_inf) / (gamma - 1) * rho * cv
gam = 2.0
p_inf, cv = 100.0, 12.5
p0 = 25.0

# -- Imposed linear temperature field T(x) = T0 + gradT*x --
T0 = 10.0
gradT = 2.0 / 15.0  # |dT/dx| = 0.13333 (same as Samareh)

rho_coeff = (p0 + p_inf) / ((gam - 1.0) * cv)  # = 10
rho_ref = rho_coeff / T0  # reference density (= 1.0 for this test)

# -- Thermal diffusivity from the Marangoni number --
# U_r = |sigma_T| * gradT * r / mu  (interfacial velocity scale)
# alpha_T = U_r * r / Ma            (thermal diffusivity, k* = 1 limit)
U_r = (-dsigma_dT) * gradT * r / mu
alpha_T = U_r * r / Ma_th
k_therm = alpha_T * rho_ref * cv * gam

# -- Grid and time stepping --
m = 128               # cells in x
Lx = 7.5              # domain length (same as the 2D box height)
dx = Lx / (m + 1)
c0 = (gam * (p0 + p_inf)) ** 0.5  # sound speed

cfl_acoustic = 0.4          # acoustic CFL safety factor
cfl_diffusion = 0.35        # explicit-diffusion CFL factor
n_dims = 1                  # 1D conduction
dt = min(cfl_acoustic * dx / c0,
         cfl_diffusion * dx**2 / (2.0 * n_dims * alpha_T))

# -- Run for n_tau viscous diffusion times --
n_tau = 3.0
tau = rho_ref * r**2 / mu  # viscous time rho * r^2 / mu
t_step_stop = int(round(n_tau * tau / dt))

rho_expr = f"{rho_coeff:.9f}/({T0} + {gradT:.9f}*x)"

data = {
    "run_time_info": "T",
    "x_domain%beg": -Lx / 2,
    "x_domain%end": Lx / 2,
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
    "bc_x%Twall_in": T0 + gradT * (-Lx / 2),
    "bc_x%Twall_out": T0 + gradT * (Lx / 2),
    # single uniform-composition patch spanning the domain, u = 0, density carries the linear T
    "patch_icpp(1)%geometry": 1,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%length_x": Lx,
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
