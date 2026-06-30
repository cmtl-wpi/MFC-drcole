#!/usr/bin/env python3
# CONTROLLED EXPERIMENT: conventional MFC "patch-per-region" version of case_Ma_20.py (TC2).
# Identical physics/EOS/grid/dt to case_Ma_20.py; the ONLY change is the initial condition layout:
#   - Patch 1 (bulk): full-box rectangle, analytic per-fluid density proxy rho_b(y) for fluid 1
#                     (so the bulk carries the imposed linear T(y), exactly as the single patch does).
#   - Patch 2 (drop): circle, smoothen=T against patch 1, with a CONSTANT (isothermal) drop density
#                     rho_d -- the textbook two-fluid drop. This is the MFC patch-per-region convention.
# Because the drop is isothermal while the bulk is stratified, T is discontinuous across the interface
# (a jump that varies around the drop) -- the prediction we are testing. Surface tension, conduction,
# sigma(T), pressure (Laplace via smoothing), color, and the bulk gradient are all unchanged, so any
# difference in the measured migration vs the validated single-patch run isolates this one choice.

import json

# Non-dimensional targets (Samareh Fig. 7)
Re = 5.0
Ma = 20.0
Ca = 0.01666
prop_ratio = 0.5  # droplet / bulk material-property ratio (rho, mu, cv, k all at 0.5)

# Geometry: D = 1 drop in a 2D wide x 4D tall box; drop 1D above the cold floor, gradient axis = y
D = 1.0
r = D / 2.0
Wx = 2.0 * D
Hy = 4.0 * D
y_bottom = -Hy / 2.0
y_drop = y_bottom + 1.0 * D

Nx = 64  # cells across the box width (Samareh used 64, 128)
dx = Wx / Nx
Ny = round(Hy / dx)

# Bulk reference scales (arbitrary; only Re, Ma, Ca are physical)
rho_b = 1.0
mu_b = 0.02
gam = 2.0

# Invert the non-dimensional numbers for the physical surface-tension and conduction properties
G = Re * mu_b**2 / (rho_b * r**2)  # |sigma_T*gradT| = 0.008
gradT = 1.0 / Hy  # T runs 1 (cold) .. 2 (hot) across the box height = 0.25
sigma_T = -G / gradT  # dsigma/dT < 0 = -0.032
sigma0 = G * r / Ca  # surface tension at T_ref ~ 0.240
alpha_b = G * r**2 / (mu_b * Ma)  # bulk thermal diffusivity = 0.005

# Imposed field T(y) = T_base + gradT*(y - y_bottom).
T_base = 1.0


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the drop's initial position; T_ref = 1.25

# Stiffened-gas EOS, bulk fluid.
c_ref = 15.0
p0 = 5.0
p_inf_b = c_ref**2 * rho_b / gam - p0  # bulk stiffening
cv_b = (p0 + p_inf_b) / ((gam - 1.0) * rho_b * T_ref)
cp_b = gam * cv_b
k_b = alpha_b * rho_b * cp_b  # bulk conductivity

# Droplet = prop_ratio * bulk for every material property.
rho_d = prop_ratio * rho_b
mu_d = prop_ratio * mu_b
cv_d = prop_ratio * cv_b
k_d = prop_ratio * k_b
p_inf_d = (gam - 1.0) * rho_d * cv_d * T_ref - p0

U_r = G * r / mu_b  # Marangoni velocity scale = 0.2
t_r = mu_b / G  # capillary-thermal time = 2.5

# Time stepping: min(acoustic CFL, explicit-conduction limit) -- identical to case_Ma_20.py.
T_hot = T_of_y(Hy / 2.0)
rho_hot_min = prop_ratio * rho_b * T_ref / T_hot
c_max = (gam * (p0 + p_inf_d) / rho_hot_min) ** 0.5
alpha_d = k_d / (rho_d * cp_b * prop_ratio)
alpha_max = max(alpha_b, alpha_d)
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_max))
t_step_stop = round(15.0 * t_r / mydt)
t_step_save = max(1, t_step_stop // 100)

# Drop center/radius (decimals -- a bare r/e/eps token would be substituted by the IC parser)
xc_d, yc_d, r_d = 0.0, y_drop, r
laplace = sigma0 / r  # Laplace pressure jump sigma/r
cb = (gam - 1.0) * cv_b
cd = (gam - 1.0) * cv_d

# Patch 1 (bulk) carries the imposed gradient analytically: rho_b(y) = (p0 + p_inf_b)/((gam-1) cv_b T(y)),
# evaluated at the bulk pressure p0 (away from the drop). alpha_1 = 1 over the whole patch.
Texpr = f"({T_base:.9f} + {gradT:.9f}*(y - ({y_bottom:.9f})))"
bulk_arho1_expr = f"({p0:.9f} + {p_inf_b:.9f})/({cb:.9f}*{Texpr})"
vac = 1.0e-8  # absent-phase partial-density floor

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Hy / 2,
    "y_domain%end": Hy / 2,
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
    "mpp_lim": "T",
    "time_stepper": 3,
    "weno_order": 7,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",  # mp_weno is WENO5-only; off for weno_order=7
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Closed slip-wall box; isothermal gradient walls set below
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + temperature-dependent surface tension sigma(T)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 (bulk)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf_b / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (droplet)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    "fluid_pp(2)%k_therm": k_d,
    # Patch 1 -- bulk background (full box). Analytic density carries the linear T(y); alpha_1 = 1.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Hy,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": bulk_arho1_expr,
    "patch_icpp(1)%alpha_rho(2)": vac,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(1)%cf_val": 0.0,
    # Patch 2 -- drop (circle), smoothed against patch 1. CONSTANT isothermal drop density rho_d:
    # the conventional two-fluid drop. This is the patch-per-region setup under test.
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": xc_d,
    "patch_icpp(2)%y_centroid": yc_d,
    "patch_icpp(2)%radius": r_d,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 + laplace,
    "patch_icpp(2)%alpha_rho(1)": vac,
    "patch_icpp(2)%alpha_rho(2)": rho_d,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%cf_val": 1.0,
    # Isothermal Dirichlet gradient walls
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Hy / 2),
    "bc_y%Twall_out": T_of_y(Hy / 2),
}

print(json.dumps(data))
