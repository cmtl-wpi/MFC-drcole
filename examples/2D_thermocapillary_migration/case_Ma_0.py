#!/usr/bin/env python3
# Thermocapillary migration of a 2D drop in the zero-Marangoni limit -- Samareh, Mostaghimi &
# Moreau, Int. J. Heat Mass Transfer 73 (2014) 616-626, Sec. 4.1.1 / Fig. 5. A neutrally buoyant
# drop of diameter D = 1 sits in an imposed linear temperature field in a 5D x 7.5D slip-wall box,
# its center 1.5D above the cold floor. Surface tension falls with temperature
# (sigma(T) = sigma0 + sigma_T*(T - T_ref), sigma_T < 0), so Marangoni stress drags the interface
# hot->cold and the drop rises toward the hot wall. With Ma = 0 the temperature is frozen at the
# imposed profile and Young-Goldstein-Block give the terminal speed
#     v_YGB = |sigma_T| |gradT| D / (6 mu_b + 9 mu_d) = 8.889e-3.
# Samareh's converged 2D ratio is v_t / v_YGB ~ 0.80.
#
# The imposed field is recovered from the stiffened-gas EOS as T = (p+p_inf)/((gam-1)*rho*cv). A
# single analytic patch sets the color blob, the drop's Laplace pressure jump, and the density from
# one shared smooth circle profile eta(x,y), so rho tracks p and the recovered T is exactly the
# imposed linear field everywhere -- drop included. Both fluids are identical, so capillary stress
# acts only on the color field. T is shifted up by T0 = 10 (Samareh's T0 = 0) to keep rho positive -- only gradT and
# sigma_T, which set v_YGB, are physical. Keep `e`-notation out of analytic strings: MFC's IC parser
# reads `e` as Euler's number.

import json

# Geometry: D = 1 drop, 5D wide x 7.5D tall box; drop 1.5D above the cold floor, rise axis = y
D = 1.0
r = D / 2.0
W = 5.0 * D
Ly = 7.5 * D
y_drop = -Ly / 2 + 1.5 * D

Nx = 64  # cells across the box width (Samareh used 64, 128, 256)
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
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2
rho_coeff = rho_drop * T0
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # closes the EOS so rho(center) = rho_drop
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop

eps = 1.0e-9  # trace volume fraction of the second phase

# One smooth circle profile eta(x,y) -- ~1 in the drop, ~0 outside -- drives the color field, the
# Laplace pressure jump, AND the density together. Then the EOS-recovered T is the imposed linear
# field BY CONSTRUCTION everywhere (drop included): rho tracks p so the (p+p_inf) factors cancel.
# The jump (sigma/r)*eta stays balanced against the surface-tension force (~sigma*kappa*grad(cf),
# kappa~1/r), so there is no t=0 acoustic ring. A single patch => no analytic-density leak across
# patches. Hardcode the drop center/radius as decimals -- a bare `r`/`e` token would be substituted.
xc_d, yc_d, r_d = 0.0, y_drop, r  # drop center and radius, as literals in the string
w_if = 0.75 * dx  # interface half-width (~3-cell transition)
laplace = sigma0 / r  # Laplace pressure jump sigma/r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
cf_expr = f"({eta})"
pres_expr = f"{p0:.9f} + {laplace:.9f}*({eta})"
# rho = rho_coeff/T * (p+p_inf)/(p0+p_inf), with p = p0 + (sigma/r)*eta, so the recovered T = T(y):
rho_num = (1.0 - eps) * rho_coeff
rho_expr = f"{rho_num:.9f}*({p0 + p_inf:.9f} + {laplace:.9f}*({eta}))/({p0 + p_inf:.9f}*({T0} + {gradT:.9f}*y))"

# Time stepping: acoustic-CFL limited (migration Mach ~ 4e-4), run to a clear terminal plateau
rho_min = rho_coeff / (T0 + gradT * Ly / 2.0)  # hot wall: lowest density, max sound speed
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
t_r = mu / abs(sigma_T * gradT)  # capillary-thermal time = 7.5
mydt = 0.35 * dx / c_max
t_step_stop = round(10.0 * t_r / mydt)  # 4 capillary-thermal times
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
    # Slip walls on all sides (Samareh's box)
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 1,
    "num_fluids": 2,
    # Physics: viscosity + temperature-dependent surface tension sigma(T)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 (continuous phase)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    # Fluid 2 (identical properties)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    # Single patch over the whole box. The drop lives entirely in the analytic profiles: cf, pres,
    # and alpha_rho all share the one eta(x,y) above, so rho tracks the (jump-carrying) pressure and
    # the EOS-recovered T is exactly the imposed linear field everywhere -- drop included, no bump.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": cf_expr,
}

print(json.dumps(data))
