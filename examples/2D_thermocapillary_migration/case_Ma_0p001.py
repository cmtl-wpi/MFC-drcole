#!/usr/bin/env python3
# Thermocapillary migration of a 2D drop in the near-zero-Marangoni limit, realized by BULK
# CONDUCTION -- Samareh, Mostaghimi & Moreau, Int. J. Heat Mass Transfer 73 (2014) 616-626,
# Sec. 4.1.1 / Fig. 5. This is the conduction companion to the frozen-T case_Ma_0.py: instead of
# pinning the temperature by faking it through density, we evolve an independent temperature scalar
# T_s with a large thermal conductivity (small Marangoni number Ma), so conduction actively holds T
# near the imposed linear gradient. As Ma -> 0 this approaches Samareh's invariant-T limit
# (v_t / v_YGB ~ 0.80). Ma is very small here (deep limit), so dt is conduction-limited (dt ~ Ma)
# and the run is long -- coarse grids only are practical.
#
# Geometry, gradient, and surface tension are TC1's, identical to case_Ma_0.py (5D x 7.5D slip-wall
# box, drop 1.5D above the cold floor, gradT = 2/15, sigma0 = 0.1, sigma_T = -0.1, identical fluids
# rho = 0.2, mu = 0.1). The conduction machinery (independent T_s scalar + isothermal gradient walls
# + k_therm from Ma) mirrors case_Ma_20.py. Keep `e`-notation out of analytic strings: MFC's IC
# parser reads `e` as Euler's number.

import json

# Marangoni number: the small-Ma realization of Samareh's Ma = 0 limit (smaller -> closer to invariant T)
Ma = 0.001

# Geometry: D = 1 drop, 5D wide x 7.5D tall box; drop 1.5D above the cold floor, rise axis = y
D = 1.0
r = D / 2.0
W = 5.0 * D
Ly = 7.5 * D
y_bottom = -Ly / 2.0
y_drop = y_bottom + 1.5 * D

Nx = 128  # cells across the box width (Samareh used 64, 128, 256)
dx = W / Nx
Ny = round(Ly / dx)

# Stiffened-gas EOS: two identical fluids (gamma = 2); low-Mach, c ~ 20 at rho = 0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
rho_b = 0.2  # Samareh: rho_d = rho_b = 0.2 (identical fluids)
mu_b = 0.1   # dynamic viscosity of both phases (MFC takes Re = 1/mu)
cv_b = 1.0   # EOS heat capacity (arbitrary; T_s is the thermal field, sigma reads T_s)

# Imposed linear field and surface tension (TC1 / Samareh Sec. 4.1.1)
gradT = 2.0 / 15.0  # |dT/dy| = 0.1333 (Samareh)
sigma0 = 0.1
sigma_T = -0.1  # dsigma/dT
G = abs(sigma_T * gradT)  # Marangoni stress scale = 0.013333

# Conduction properties from Ma: smaller Ma = larger diffusivity = T held closer to the gradient
alpha_b = G * r**2 / (mu_b * Ma)  # bulk thermal diffusivity
cp_b = gam * cv_b
k_b = alpha_b * rho_b * cp_b  # bulk conductivity
rho_d, mu_d, cv_d, k_d = rho_b, mu_b, cv_b, k_b  # identical fluids

# Imposed field T(y) = T0 + gradT*y, carried by the scalar T_s. T0 = 10 keeps T well above 0 for the
# isothermal-wall validator; only gradT and sigma_T (which set v_YGB) are physical.
T0 = 10.0


def T_of_y(y):
    return T0 + gradT * y


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the drop's initial position
T_expr = f"{T0} + {gradT:.9f}*y"  # plain decimals only (no `e` notation)

eps = 1.0e-9  # trace volume fraction of the "other" fluid in each patch

# Time stepping: min(acoustic CFL, explicit-conduction limit). Small Ma => conduction-limited dt.
c_max = (gam * (p0 + p_inf) / rho_b) ** 0.5  # uniform density (T carried by T_s, not by density)
t_r = mu_b / G  # capillary-thermal time = 7.5
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_b))  # 2D explicit diffusion number
t_step_stop = round(2.0 * t_r / mydt)  # 2 capillary-thermal times (the conduction-sweep window)
t_step_save = max(1, t_step_stop // 80)

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation Algorithm
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Slip walls on all sides (Samareh's box); isothermal gradient walls on y set below
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + sigma(T); T carried by an independent scalar T_s
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    "thermal_scalar": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "T_s_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 (continuous phase)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (identical properties)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    "fluid_pp(2)%k_therm": k_d,
    # Patch 1: background medium, uniform density, color c = 0, linear T_s(y)
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_b,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_d,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%T_temp_val": T_expr,
    # Patch 2: droplet (circle), color c = 1, same T_s(y) as the bulk (T is continuous across the
    # interface); identical fluid to patch 1. Pressure carries the Laplace jump p0 + sigma/r.
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 + sigma0 / r,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_d,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_expr,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Ly / 2),
    "bc_y%Twall_out": T_of_y(Ly / 2),
}

print(json.dumps(data))
