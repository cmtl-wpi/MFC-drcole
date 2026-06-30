#!/usr/bin/env python3
# Thermocapillary migration of a single 2D drop in the CREEPING limit -- Nas & Tryggvason,
# "Thermocapillary interaction of two bubbles or drops", Int. J. Multiphase Flow 29 (2003)
# 1117-1135, Sec. 3.1 / Fig. 2 (resolution test and validation). A drop of diameter D = 1 sits
# centered in a SQUARE box "extending four drop radii in each direction" (2D x 2D), in an imposed
# linear temperature field (cold floor, hot ceiling). Surface tension falls with temperature, so
# Marangoni stress drags the interface hot->cold and the drop migrates toward the hot wall. The
# migration velocity V* = V/U_r should rise monotonically to a steady-state plateau ~ 0.13.
#   Non-dimensional targets: Re = Ma = 2.5e-3, Ca = 1e-3; ALL property ratios (rho, mu, cv, k) = 0.5.
#
# STATUS: NOT yet reproduced. The committed Nx=32 run completes cleanly but is acoustically
# contaminated and ~3x fast (mean V* ~ 0.4): the very low c_ref below (chosen to keep dt large in this
# conduction-step-limited creeping run, ~1.2M steps) lets the IC's acoustic ringing swamp the tiny
# (Ma~1e-5) migration. A faithful run needs c_ref ~ 15 AND a finer grid (~5x the steps, ~day-scale).
# (Nas-Tryggvason Fig. 3, the finite-Re/Ma test Re=5/Ma=20/Ca=0.01666, is already case_Ma_20.py.)
#
# GEOMETRY matches Nas & Tryggvason Fig. 1: periodic in x, no-slip isothermal walls on the
# gradient axis (here +y). This differs from the Samareh slip-wall box used by case_Ma_0/20; in the
# tight 4-radius square the wall type sets the confined terminal velocity, and 0.13 is the no-slip
# value. The conduction isothermal-wall BC (m_thermal_conduction.fpp) fires for any physical face
# (bc%beg/end < 0), so it works with no-slip (-16) exactly as with reflective (-2).
#
# TEMPERATURE WITH DISTINCT FLUIDS (compressible code, incompressible reference) -- same construction
# as case_Ma_20.py. MFC is compressible: at this vanishing migration Mach pressure is ~uniform and
# the EOS locks T = (p+p_inf)/((gam-1)*rho*cv), so a temperature gradient is a density gradient. With
# two distinct fluids we impose T(y) by a PER-FLUID density proxy,
#   rho_i(x,y) = (p + p_inf_i)/((gam-1)*cv_i*T(y)),
# so the fluid-to-fluid density RATIO is height-independent (T(y) cancels) and the Nas-Tryggvason
# ratio 0.5 is preserved by the per-fluid stiffening p_inf_i; the MIXTURE EOS recovers T_mix(y)=T(y)
# exactly across the interface. Bulk conduction (here FAST: Ma=2.5e-3 => large alpha) plus isothermal
# walls hold the field at the imposed profile against the drop's (tiny) advection. Honest caveat:
# each fluid's ABSOLUTE density stratifies as ~1/T (a compressibility artifact absent in the
# front-tracking reference), so the plateau matches Fig 2 qualitatively; the monotonic-rise-to-plateau
# shape and the ~0.13 level are the validation target. Keep `e`-notation out of analytic strings.
#
# COST: creeping flow (Re=Ma=2.5e-3) => t_r/tau_diff = 1/Re = 400, so reaching the plateau takes
# O(1/Re) diffusion times and the run is conduction-dt-limited: ~1.2M steps at Nx=32 (16 cells/D, the
# coarsest paper grid), ~4.7M at Nx=64, ~19M at Nx=128. Start coarse.

import json

# Non-dimensional targets (Nas & Tryggvason Fig. 2)
Re = 2.5e-3
Ma = 2.5e-3
Ca = 1.0e-3
prop_ratio = 0.5  # drop / bulk material-property ratio (rho, mu, cv, k all at 0.5)

# Geometry: D = 1 drop centered in a square box 4 drop radii on a side (= 2D x 2D); gradient axis = y
D = 1.0
r = D / 2.0
Wx = 2.0 * D
Hy = 2.0 * D
y_bottom = -Hy / 2.0
y_drop = 0.0  # drop centered (Fig 2 is a symmetric square resolution test)

Nx = 32  # cells across the box width (Nas & Tryggvason: 32, 64, 128 = 16, 32, 64 cells/D)
dx = Wx / Nx
Ny = round(Hy / dx)

# Bulk reference scales (arbitrary; only Re, Ma, Ca are physical)
rho_b = 1.0
mu_b = 0.02
gam = 2.0

# Invert the non-dimensional numbers for the physical surface-tension and conduction properties
G = Re * mu_b**2 / (rho_b * r**2)  # |sigma_T*gradT| = 4e-6
gradT = 1.0 / Hy  # T runs 1 (cold) .. 2 (hot) across the box height = 0.5
sigma_T = -G / gradT  # dsigma/dT < 0 = -8e-6
sigma0 = G * r / Ca  # surface tension at T_ref = 2e-3
alpha_b = G * r**2 / (mu_b * Ma)  # bulk thermal diffusivity = 0.02

# Imposed field T(y) = T_base + gradT*(y - y_bottom). T_base > 0 keeps rho positive and the
# isothermal-wall validator happy; a larger offset reduces the ~1/T compressibility stratification.
T_base = 1.0


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the drop's initial position; T_ref = 1.5

# Stiffened-gas EOS, bulk fluid. The migration Mach (U_r/c ~ 5e-5) leaves enormous headroom, so a low
# reference sound speed is chosen to keep the acoustic step from binding below the conduction step
# (creeping flow is conduction-dt-limited). p0 sits well above the Laplace jump sigma0/r = 4e-3.
c_ref = 2.5
p0 = 0.5
p_inf_b = c_ref**2 * rho_b / gam - p0  # bulk stiffening = 2.625
cv_b = (p0 + p_inf_b) / ((gam - 1.0) * rho_b * T_ref)
cp_b = gam * cv_b
k_b = alpha_b * rho_b * cp_b  # bulk conductivity

# Droplet = prop_ratio * bulk for every material property. The EOS stiffening p_inf_d follows so the
# drop density is rho_d at (p0, T_ref): ratio rho_d/rho_b = prop_ratio.
rho_d = prop_ratio * rho_b
mu_d = prop_ratio * mu_b
cv_d = prop_ratio * cv_b
k_d = prop_ratio * k_b
p_inf_d = (gam - 1.0) * rho_d * cv_d * T_ref - p0

U_r = G * r / mu_b  # Marangoni velocity scale = 1e-4
t_r = mu_b / G  # capillary-thermal time = 5000

# Time stepping: min(acoustic CFL, explicit-conduction limit). Density is lowest (sound speed
# highest) at the hot wall, where rho ~ rho_ref*T_ref/T_hot.
T_hot = T_of_y(Hy / 2.0)
rho_hot_min = prop_ratio * rho_b * T_ref / T_hot  # drop fluid at the hot wall: global density min
c_max = (gam * (p0 + p_inf_d) / rho_hot_min) ** 0.5
alpha_d = k_d / (rho_d * cp_b * prop_ratio)  # = k_d/(rho_d*cp_d), the fastest-diffusing phase
alpha_max = max(alpha_b, alpha_d)
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_max))
t_step_stop = round(2.0 * t_r / mydt)  # 2 capillary-thermal times (plateau is reached by t* ~ 1.2)
t_step_save = max(1, t_step_stop // 100)

# One analytic patch over the whole box (a second patch's analytic density would leak globally). The
# smooth circle eta(x,y) ~ 1 in the drop / 0 outside drives the color, the volume-fraction split, the
# Laplace pressure jump, and BOTH per-fluid densities together. Hardcode the drop center/radius as
# decimals -- a bare r/e/eps token would be substituted by the IC parser.
xc_d, yc_d, r_d = 0.0, y_drop, r
w_if = 0.75 * dx  # interface half-width (~3-cell transition)
laplace = sigma0 / r  # Laplace pressure jump sigma/r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
Texpr = f"({T_base:.9f} + {gradT:.9f}*(y - ({y_bottom:.9f})))"
pres_expr = f"({p0:.9f} + {laplace:.9f}*({eta}))"
# volume fractions: drop = eta, bulk = 1 - eta
alpha2_expr = f"({eta})"
alpha1_expr = f"(1.0 - ({eta}))"
# partial densities: alpha_i * rho_i, with rho_i = (p + p_inf_i)/((gam-1)*cv_i*T(y))
cb = (gam - 1.0) * cv_b
cd = (gam - 1.0) * cv_d
arho1_expr = f"(1.0 - ({eta}))*({pres_expr} + {p_inf_b:.9f})/({cb:.9f}*{Texpr})"
arho2_expr = f"({eta})*({pres_expr} + {p_inf_d:.9f})/({cd:.9f}*{Texpr})"
cf_expr = f"({eta})"

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
    "weno_order": 5,  # weno5 needs >=25 cells/dir; the 32-cell square (16/D) cannot fit weno7 (>=35)
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "F",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Nas & Tryggvason Fig 1: periodic in x, no-slip isothermal walls on the gradient axis (+y)
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "num_patches": 1,
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
    # Fluid 2 (droplet), all material properties at prop_ratio * bulk; stiffening keeps rho_d/rho_b = 0.5
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    "fluid_pp(2)%k_therm": k_d,
    # Single analytic patch: the drop lives entirely in eta(x,y). Color, volume fractions, pressure,
    # and both per-fluid densities all share it, so the mixture EOS recovers the linear T(y) exactly
    # (drop included) while the distinct fluid-to-fluid density ratio (0.5) is preserved.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Hy,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": arho1_expr,
    "patch_icpp(1)%alpha_rho(2)": arho2_expr,
    "patch_icpp(1)%alpha(1)": alpha1_expr,
    "patch_icpp(1)%alpha(2)": alpha2_expr,
    "patch_icpp(1)%cf_val": cf_expr,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Hy / 2),
    "bc_y%Twall_out": T_of_y(Hy / 2),
}

print(json.dumps(data))
