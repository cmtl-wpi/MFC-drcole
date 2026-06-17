#!/usr/bin/env python3
# Thermocapillary droplet migration, 2D -- validation against Samareh, Mostaghimi & Moreau,
# Int. J. Heat Mass Transfer 73 (2014) 616-626 (Sec. 4.1.1 / Fig 5, "zero Marangoni number").
#
# A neutrally-buoyant drop sits in an imposed linear vertical temperature field. Surface tension
# falls with temperature (sigma(T) = sigma0 + sigma_T*(T - T_ref), sigma_T < 0), so Marangoni
# stress drags the interface hot->cold and the drop rises toward the hot top. In the zero-Ma
# creeping-flow limit Young-Goldstein-Block give the terminal speed
#     v_YGB = (2/15)*|sigma_T|*gradT*(D/2)/mu;   Samareh's converged 2D ratio v_t/v_YGB ~ 0.80.
#
# The imposed T is encoded through density at uniform pressure (rho = rho_coeff/T) and recovered
# from the stiffened-gas EOS; both fluids are identical, so capillary stress acts only on the color
# field. T is shifted up by T0=10 to keep rho positive -- only gradT and sigma_T (which set v_YGB)
# are physical. 2D only: the 3D rise has no quasi-steady plateau on this branch (drifts with grid).
#
# Env vars (one build serves the whole sweep; the analytic IC depends only on the fixed T0/gradT):
#   SAMAREH_NX    cells across the 5D box width (default 128)
#   SAMAREH_DSDT  dsigma/dT slope (default -0.1)
#   SAMAREH_TR    run length in capillary-thermal times t_r (default 4)
#   SAMAREH_WALL  1 = slip-wall box, drop 1.5D above cold floor (default); 0 = open/centered
#   SAMAREH_MA    thermal Marangoni number; >0 enables bulk Fourier conduction (default 0 = frozen T)
#   SAMAREH_TS    1 = carry T as an independent advected scalar at uniform density (default 0)
#   SAMAREH_UNBALANCED  1 = bare uniform-pressure IC (rings); default balances the t=0 Laplace jump
#
# Note: keep `e`-notation out of analytic patch strings -- MFC's IC parser reads `e` as Euler's number.

import json
import os

# Variant selection (defaults = 2D headline case at Samareh's medium grid)
width_cells = int(os.environ.get("SAMAREH_NX", "128"))
dsigma_dT = float(os.environ.get("SAMAREH_DSDT", "-0.1"))
n_tr = float(os.environ.get("SAMAREH_TR", "4"))
wall = os.environ.get("SAMAREH_WALL", "1") == "1"
Ma_th = float(os.environ.get("SAMAREH_MA", "0"))  # >0 enables bulk conduction
ts_mode = os.environ.get("SAMAREH_TS", "0") == "1"  # independent temperature scalar
unbalanced_ic = os.environ.get("SAMAREH_UNBALANCED", "0") == "1"
assert Ma_th >= 0, "SAMAREH_MA must be >= 0 (0 disables conduction)"

# Geometry: D=1 drop, 5D wide x 7.5D tall box, rise axis = y
D = 1.0
r = D / 2.0
W = 5.0 * D
Ly = 7.5 * D
y_drop = (-Ly / 2 + 1.5 * D) if wall else 0.0  # 1.5D above cold floor (wall) or centered (open)
bc = -2 if wall else -3  # slip walls, else open

dx = W / width_cells
long_cells = round(Ly / dx)
m = width_cells - 1
n = long_cells - 1
p = 0

cf_smooth_coeff = 0.5  # ~2-cell diffuse color interface

# Equation of state: two identical stiffened-gas fluids (gamma=2); low-Mach, c ~ 20 at rho=0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1  # dynamic viscosity of both phases (MFC takes Re = 1/mu)

# Imposed linear field T(y) = T0 + gradT*y, encoded as rho(y) = rho_coeff/T(y) at uniform p0
T0 = 10.0  # temperature at box center (shifted up from Samareh's T=0 to keep rho positive)
gradT = 2.0 / 15.0  # |dT/dy| = 0.1333 (Samareh)
sigma0 = 0.1
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2
rho_coeff = rho_drop * T0
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # closes the EOS so rho(center) = 0.2
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop

eps = 1.0e-9  # trace volume fraction of the second phase
rho_num = (1.0 - eps) * rho_coeff  # plain decimal: eps in the string would render "1e-09" (parser footgun)
if ts_mode:
    rho_expr = (1.0 - eps) * rho_coeff / T0  # uniform density; T imposed on the scalar instead
    T_expr = f"{T0} + {gradT:.9f}*y"
else:
    rho_expr = f"{rho_num:.9f}/({T0} + {gradT:.9f}*y)"

# Stability / time scales
rho_min = rho_coeff / (T0 + gradT * (Ly / 2.0))  # hot top: lowest density -> max sound speed
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
t_r = mu / abs(dsigma_dT * gradT)  # capillary-thermal time = 7.5

# Bulk conduction: k from the requested Marangoni number, alpha = U_r*r/Ma
if Ma_th > 0:
    U_r = (-dsigma_dT) * gradT * r / mu
    alpha_T = U_r * r / Ma_th
    k_therm = alpha_T * (rho_coeff / T0) * cv * gam

# Time stepping: acoustic-CFL limited (Mach ~ 4e-4); add the diffusion cap when conducting
mydt = 0.35 * dx / c_max
if Ma_th > 0:
    mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_T))
t_end = n_tr * t_r
t_step_stop = int(round(t_end / mydt))
t_step_save = max(1, t_step_stop // 80)  # ~80 snapshots

data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation Algorithm Parameters
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
    # Boundaries: slip walls (Samareh's box) or open (approximates an unbounded domain)
    "bc_x%beg": bc,
    "bc_x%end": bc,
    "bc_y%beg": bc,
    "bc_y%end": bc,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + temperature-dependent surface tension sigma(T)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": dsigma_dT,
    # Formatted Database Files Structure Parameters
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
    # Second phase (fluid 2), identical properties
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    # Patch 1: background medium, analytic linear-T density, color c=0
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    # Patch 2: droplet (circle), color c=1, smeared over ~2 cells. Same density/composition as
    # patch 1 (no real fluid jump). Pressure is the Laplace overpressure p0 + sigma/r unless
    # unbalanced, which removes the t=0 acoustic ring.
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": cf_smooth_coeff,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 if unbalanced_ic else p0 + sigma0 / r,
    "patch_icpp(2)%alpha_rho(1)": rho_expr,
    "patch_icpp(2)%alpha_rho(2)": eps,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1.0,
}

if Ma_th > 0:
    # Bulk Fourier conduction (equal k -> the YGB k* = 1 limit)
    data.update({"thermal_conduction": "T", "fluid_pp(1)%k_therm": k_therm, "fluid_pp(2)%k_therm": k_therm})
    if wall:
        # Isothermal Dirichlet walls pin the cold floor / hot ceiling to the imposed gradient
        data.update(
            {
                "bc_y%isothermal_in": "T",
                "bc_y%isothermal_out": "T",
                "bc_y%Twall_in": T0 + gradT * (-Ly / 2.0),
                "bc_y%Twall_out": T0 + gradT * (Ly / 2.0),
            }
        )

if ts_mode:
    # Independent temperature scalar: T(y) imposed on T_s for both patches; sigma(T) reads T_s
    data.update(
        {
            "thermal_scalar": "T",
            "T_s_wrt": "T",
            "patch_icpp(1)%T_temp_val": T_expr,
            "patch_icpp(2)%T_temp_val": T_expr,
        }
    )

print(json.dumps(data))
