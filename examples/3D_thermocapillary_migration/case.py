#!/usr/bin/env python3
# 3D thermocapillary droplet migration -- Samareh, Mostaghimi & Moreau (2014), Fig 6:
# a fully 3D SPHERE in an imposed linear temperature gradient, target v_t/v_YGB ~ 0.95.
#
# This is the 3D sibling of ../2D_thermocapillary_migration (which reproduces Samareh's Fig 5, the
# planar 2D *cylinder* -> ~0.80). Same physics and EOS realization; the only changes are the third
# dimension (a sphere patch, p > 0) and -- crucially -- BULK CONDUCTION IS MANDATORY here.
#
# WHY CONDUCTION IS REQUIRED IN 3D (and optional in 2D)
# In the no-conduction (frozen-T) limit the 2D rise reaches a quasi-steady plateau, but the 3D rise
# does NOT: the toroidal internal circulation continuously steepens the frozen interfacial gradient,
# so the velocity drifts past v_YGB without saturating (finer grid -> faster). There is then no
# validatable 3D number. Bulk Fourier conduction (the thermal_conduction feature) diffuses the
# temperature, tames that runaway, and restores a steady plateau that can be compared to Samareh's
# 0.95. So unlike the 2D example, this one runs with thermal_conduction = T and isothermal walls by
# construction (SAMAREH3D_MA > 0 is enforced).
#
# PARAMETERS (Samareh Sec. 4.1.1, identical to the 2D case): D = 1 sphere, box 5D x 5D x 7.5D,
# rho_d = rho_b = 0.2, mu_d = mu_b = 0.1, sigma0 = 0.1, sigma_T = -0.1, |gradT| = 2/15 -> v_YGB = 8.889e-3.
# The drop sits 1.5D above the cold (y%beg) floor; slip walls on all six faces. The one deviation is
# the absolute temperature baseline (shifted up by T0 = 10 so the density proxy rho = rho_coeff/T stays
# positive); the gradient and the slope sigma_T -- all v_YGB and the Marangoni stress depend on -- are exact.
#
# GOTCHA (analytic-IC parser): MFC expands a bare `e` to Euler's number even inside `1e-9`. Keep
# `e`-notation out of every analytic patch string (see rho_expr); eps is folded into a plain decimal.

import json
import os

# -- Variant selection (env vars) --
width_cells = int(os.environ.get("SAMAREH3D_NX", "64"))  # cells per box WIDTH (Samareh 3D used 64, 128)
Ma_th = float(os.environ.get("SAMAREH3D_MA", "1.0"))  # thermal Marangoni number; conduction REQUIRED in 3D
n_tr = float(os.environ.get("SAMAREH3D_TR", "2"))  # run length in capillary-thermal times t_r
assert Ma_th > 0, "3D needs conduction (SAMAREH3D_MA > 0): the frozen-T 3D rise has no plateau to validate."

# -- Geometry (Samareh Fig 6: D=1 sphere, 5D x 5D x 7.5D box; rise axis = y) --
D = 1.0  # droplet diameter
r = D / 2.0  # droplet radius = 0.5
W = 5.0 * D  # short-axis extent (x and z), 2.5D clearance each side
Ly = 7.5 * D  # long-axis extent (the gradient / rise axis, +y)
y_drop = -Ly / 2 + 1.5 * D  # Samareh: drop 1.5D above the cold floor = -2.25

dx = W / width_cells  # isotropic cell size set by the box-width resolution
long_cells = round(Ly / dx)  # cells along the 7.5D rise axis (= 1.5*width_cells)
m = width_cells - 1  # x (short)
n = long_cells - 1  # y (long, rise axis)
p = width_cells - 1  # z (short) -- 3D

# Diffuse color interface ~2 cells wide at every resolution (smooth_coeff = dx/w with w = 2*dx -> 0.5).
cf_smooth_coeff = 0.5

# -- Equation of state (two IDENTICAL stiffened-gas fluids, gamma = 2; mu* = 1, k* = 1) --
gam = 2.0
p_inf, p0 = 32.0, 8.0  # with rho ~ 0.2: c = sqrt(gam*(p0+p_inf)/rho) ~ 20
mu = 0.1  # dynamic viscosity of both phases; MFC takes Re = 1/mu

# -- Imposed linear temperature field T(y) = T0 + gradT*y (centered on the box at y=0) --
T0 = 10.0  # temperature at the box center; shifted up from Samareh's T(cold wall)=0 so rho stays positive
gradT = 2.0 / 15.0  # |dT/dy| = 0.13333, Samareh's imposed gradient
sigma0 = 0.1  # surface tension at T_ref
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2
rho_coeff = rho_drop * T0  # = 2.0 ; rho(center) = rho_coeff/T0 = 0.2
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # = 62.5 ; closes the EOS so rho(center) = 0.2
T_ref = T0 + gradT * y_drop  # sigma = sigma0 AT THE DROP
dsigma_dT = -0.1  # dsigma/dT (Samareh)

eps = 1.0e-9  # trace volume fraction of the (identical) second phase
rho_num = (1.0 - eps) * rho_coeff  # = 1.999999998 (folded to a plain decimal; never embed eps=1e-9 in a string)
rho_expr = f"{rho_num:.9f}/({T0} + {gradT:.9f}*y)"  # ~ 2.0/(10.0 + 0.133333333*y)

# -- Young-Goldstein-Block terminal speed and the diagnostic time scales --
v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu  # = 8.889e-3
rho_min = rho_coeff / (T0 + gradT * (Ly / 2.0))  # density at the hot top wall (max sound speed)
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5  # ~ 20.5
tau = rho_drop * r**2 / mu  # viscous time = 0.5 (Samareh's)
t_r = mu / abs(dsigma_dT * gradT)  # capillary-thermal time = 7.5 (Samareh scale)

# -- Bulk thermal conduction (k* = 1): alpha_T from the requested thermal Marangoni number Ma --
U_r = (-dsigma_dT) * gradT * r / mu  # Marangoni interfacial velocity scale
alpha_T = U_r * r / Ma_th  # thermal diffusivity
k_therm = alpha_T * (rho_coeff / T0) * cv * gam  # k = alpha*rho*cp with cp = gam*cv at rho = rho_coeff/T0

# -- Time stepping: acoustic-CFL limited, also capped by the 3D explicit-diffusion number (d = 3) --
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (6.0 * alpha_T))  # dt <= 0.35*dx^2/(2*d*alpha), d = 3 (3D)
t_end = n_tr * t_r
t_step_stop = int(round(t_end / mydt))
t_step_save = max(1, t_step_stop // 60)  # ~60 snapshots

data = {
    "run_time_info": "T",
    # Computational domain: rise (gradient) axis y in [-Ly/2, Ly/2]; short axes x,z in [-W/2, W/2]
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "z_domain%beg": -W / 2,
    "z_domain%end": W / 2,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation algorithm (6-equation model; proven WENO5/HLLC settings)
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
    # Boundaries: slip walls (-2) on all six faces (Samareh's box)
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + temperature-dependent surface tension sigma(T) = sigma0 + sigma_dTdT*(T-T_ref)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": dsigma_dT,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Continuous phase (fluid 1)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    # Second phase (fluid 2) -- identical properties (mu* = 1, k* = 1)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    # Patch 1 -- background medium (3D cuboid spanning the domain): analytic linear-T density, color c=0.
    "patch_icpp(1)%geometry": 9,  # 3D cuboid
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": W,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    # Patch 2 -- droplet (3D sphere): marks c=1, smeared over ~2 cells. IDENTICAL density/composition/
    # pressure to patch 1, so the capillary stress acts purely on the c interface (the mu*=1, k*=1 YGB limit).
    "patch_icpp(2)%geometry": 8,  # 3D sphere
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": cf_smooth_coeff,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": rho_expr,
    "patch_icpp(2)%alpha_rho(2)": eps,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1.0,
    # Bulk Fourier conduction (REQUIRED in 3D) + isothermal Dirichlet y-walls pinned to the imposed
    # gradient (cold floor / hot ceiling). This is what tames the frozen-T 3D runaway into a plateau.
    "thermal_conduction": "T",
    "fluid_pp(1)%k_therm": k_therm,
    "fluid_pp(2)%k_therm": k_therm,
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T0 + gradT * (-Ly / 2.0),  # cold floor (y%beg)
    "bc_y%Twall_out": T0 + gradT * (Ly / 2.0),  # hot ceiling (y%end)
}

print(json.dumps(data))
