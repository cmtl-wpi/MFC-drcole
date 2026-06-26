#!/usr/bin/env python3
# CONTROLLED EXPERIMENT: conventional MFC "patch-per-region" version of case_Ma_0p1.py (TC1,
# conduction at Ma=0.1). Identical physics/EOS/grid/dt/algorithm to case_Ma_0p1.py; the ONLY
# change is the initial-condition layout:
#   - Patch 1 (bulk): full-box rectangle, analytic density proxy rho_b(y) = rho_coeff/T(y) for
#                     fluid 1, so the bulk carries the imposed linear T(y) exactly as the single
#                     patch does.
#   - Patch 2 (drop): circle, smoothen=T against patch 1, with a CONSTANT (isothermal) drop
#                     density rho_d = rho_coeff/T_ref -- the textbook two-fluid drop, fluid 2.
# In TC1 both fluids are identical (rho_d = rho_b, equal mu/cv/gamma), so at t=0 the constant drop
# density equals the bulk density at the drop centroid -- there is no initial jump. But the drop is
# ISOTHERMAL while the bulk is STRATIFIED, so as the drop rises into warmer (lower-density) bulk a
# temperature jump opens across the interface -- exactly where the Marangoni force lives. This is the
# failure mode the single-patch layout was built to avoid; this case measures its cost on Fig 5.

import json

# Marangoni number: the small-Ma realization of Samareh's Ma = 0 limit
Ma = 0.1

# Geometry: D = 1 drop, 5D wide x 7.5D tall box; drop 1.5D above the cold floor, rise axis = y
D = 1.0
r = D / 2.0
W = 5.0 * D
Ly = 7.5 * D
y_drop = -Ly / 2.0 + 1.5 * D

Nx = 64  # cells across the box width
dx = W / Nx
Ny = round(Ly / dx)

# Stiffened-gas EOS: two identical fluids (gamma = 2); low-Mach, c ~ 20 at rho = 0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1  # dynamic viscosity of both phases (MFC takes Re = 1/mu)

# Imposed linear field T(y) = T0 + gradT*y, encoded as rho(y) = rho_coeff/T(y) at uniform p0
T0 = 10.0
gradT = 2.0 / 15.0  # |dT/dy| = 0.1333 (Samareh)
sigma0 = 0.1
sigma_T = -0.1  # dsigma/dT
G = abs(sigma_T * gradT)  # Marangoni stress scale = 0.013333
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2 (density at the drop's reference temperature T0)
rho_coeff = rho_drop * T0
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # closes the EOS so rho(T0) = rho_drop
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop's initial position

# Conduction properties from Ma: smaller Ma = larger diffusivity = T held closer to the gradient.
alpha_b = G * r**2 / (mu * Ma)  # bulk thermal diffusivity
cp = gam * cv
k_b = alpha_b * rho_drop * cp  # bulk conductivity

# Drop center/radius (decimals -- a bare r/e token would be substituted by the IC parser)
xc_d, yc_d, r_d = 0.0, y_drop, r
laplace = sigma0 / r  # Laplace pressure jump sigma/r
cb = (gam - 1.0) * cv  # so rho_b(y) = (p0 + p_inf)/(cb*T(y))

# Patch 1 (bulk) carries the imposed gradient analytically: rho_b(y) = (p0 + p_inf)/((gam-1) cv T(y)),
# evaluated at the bulk pressure p0. alpha_1 = 1 over the whole patch.
Texpr = f"({T0:.9f} + {gradT:.9f}*y)"
bulk_arho1_expr = f"({p0 + p_inf:.9f})/({cb:.9f}*{Texpr})"
rho_d = (p0 + p_inf) / (cb * T_ref)  # constant isothermal drop density (= bulk density at y_drop)
vac = 1.0e-8  # absent-phase partial-density floor


def T_of_y(y):
    return T0 + gradT * y


# Time stepping: min(acoustic CFL, explicit-conduction limit) -- identical to case_Ma_0p1.py.
rho_min = rho_coeff / (T0 + gradT * Ly / 2.0)  # hot wall: lowest density, max sound speed
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
t_r = mu / G  # capillary-thermal time = 7.5
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_b))  # 2D explicit diffusion number
t_step_stop = round(10.0 * t_r / mydt)
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
    # Simulation Algorithm (identical to the single-patch parent: mpp_lim=F isolates the IC change)
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
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
    # Open (ghost-cell extrapolation) side walls in x; isothermal gradient walls on y (set below).
    # The y walls MUST stay walls: they pin the cold-floor/hot-ceiling Dirichlet T that drives the
    # gradient. Opening them with -3 + isothermal is documented to reverse the migration here.
    "bc_x%beg": -3,
    "bc_x%end": -3,
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
    # Fluid 1 (bulk / continuous phase)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (droplet -- identical properties to fluid 1)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%k_therm": k_b,
    # Patch 1 -- bulk background (full box). Analytic density carries the linear T(y); alpha_1 = 1.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
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
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Ly / 2),
    "bc_y%Twall_out": T_of_y(Ly / 2),
}

print(json.dumps(data))
