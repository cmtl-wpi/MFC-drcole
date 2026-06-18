#!/usr/bin/env python3
# Thermocapillary migration of a 3D drop at large Marangoni number -- Samareh, Mostaghimi & Moreau,
# Int. J. Heat Mass Transfer 73 (2014) 616-626, Sec. 4.2 / Figs. 8, 13, matched to the LMS Space
# Shuttle microgravity experiment (Hadland et al.). A Fluorinert FC-75 drop (D = 10.7 mm) rises
# through DC-200 silicon oil in a 45 x 60 x 45 mm cell with a 1000 K/m vertical gradient (cold floor
# 283 K, hot ceiling 343 K). At the T0 = 313 K cell center Re = 17.79 and Ma = 1723.
#
# Three coupled features are active: bulk Fourier conduction (thermal_conduction), the sigma(T)
# closure (sigma_model = 1), and temperature-dependent viscosity mu(T) = exp(C + D/T) (visc_model =
# 1) -- the silicon oil's viscosity varies substantially across the 60 K cell, and that drag produces
# the experiment's non-monotonic rise-velocity loop. Drop and bulk have distinct near-constant
# densities, so temperature is carried by an independent scalar T_s (in Kelvin); sigma(T) and mu(T)
# read it. The stiffened-gas EOS only sets the (low-Mach, insensitive) acoustics. Sec. 4.2.1: the
# drop starts on the local linear field (Fig. 8). A converged comparison is heavy (Samareh used
# 120 x 320 x 120 up to 240 x 640 x 240); this case runs a coarse demo grid by default.

import json
import math

# Geometry (LMS test cell; gradient / rise axis = y)
D = 10.7e-3
r = D / 2.0
Wx = 45.0e-3  # square cross-section (x, z)
Ly = 60.0e-3  # gradient / rise axis (y)
y_drop = -Ly / 2 + 15.0e-3  # released 15 mm above the cold wall

Nx = 30  # cells across the 45 mm cross-section (Samareh used 120, 180, 240)
dx = Wx / Nx
Ny = round(Ly / dx)

# Imposed temperature field (Kelvin)
T_c, T_h = 283.0, 343.0  # cold floor / hot ceiling wall temperatures
gradT = (T_h - T_c) / Ly  # = 1000 K/m
T0 = 313.0  # reference temperature (cell center)
sigma0 = 0.007  # surface tension at T0 (N/m)
sigma_T = -3.6e-5  # dsigma/dT (N/m/K)

# Per-fluid properties (fluid 1 = silicon oil bulk, fluid 2 = Fluorinert drop); mu(T) = exp(C + D/T)
rho_b, rho_d = 918.3, 1727.7  # densities (kg/m^3)
k_b, k_d = 0.13389, 0.063  # thermal conductivities (W/mK)
cp_b, cp_d = 1778.2, 1047.0  # specific heats (J/kgK)
C_b, D_b = -10.17, 1643.0  # silicon oil Arrhenius coefficients
C_d, D_d = -11.76, 1540.0  # Fluorinert Arrhenius coefficients

# Stiffened-gas EOS, softened for a stable low-Mach state (acoustics decoupled from T via T_s)
gam = 2.0
cv_b, cv_d = cp_b / gam, cp_d / gam  # cp = gam*cv
c_snd = 30.0  # softened sound speed (m/s); migration Mach ~ 1e-3, so c only sets the acoustic CFL
p0 = 1.0e5  # background pressure (~1 atm)
p_inf_b = rho_b * c_snd**2 / gam - p0  # stiffening so c = sqrt(gam*(p0+p_inf)/rho) = c_snd
p_inf_d = rho_d * c_snd**2 / gam - p0

# Reference viscosity at the drop-center temperature at t = 0 (visc_model = 1 evaluates mu(T) live)
T_visc_ref = T0 + gradT * y_drop  # = 298 K
mu_b_ref = math.exp(C_b + D_b / T_visc_ref)
mu_d_ref = math.exp(C_d + D_d / T_visc_ref)

# Migration scales (Samareh): U_r = |sigma_T*gradT|*r/mu_b(T0), t_r = mu_b(T0)/|sigma_T*gradT|
mu_b0 = math.exp(C_b + D_b / T0)  # silicon oil viscosity at T0 (~7.3e-3 Pa.s)
G = abs(sigma_T * gradT)  # Marangoni stress scale
t_r = mu_b0 / G  # capillary-thermal time ~ 0.20 s

# Time stepping: acoustic CFL + 3D explicit-diffusion cap on the smallest cell
alpha_b, alpha_d = k_b / (rho_b * cp_b), k_d / (rho_d * cp_d)
mydt = 0.35 * dx / c_snd
mydt = min(mydt, 0.35 * dx**2 / (6.0 * max(alpha_b, alpha_d)))
t_step_stop = round(1.0 * t_r / mydt)  # 1 capillary-thermal time
t_step_save = max(1, t_step_stop // 60)

# Temperature-scalar IC: both phases start on the imposed linear field (Sec. 4.2.1)
T_field = f"{T0} + {gradT:.6f}*y"  # plain decimals only (no `e` notation)

eps = 1.0e-9  # trace volume fraction of the other phase

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "z_domain%beg": -Wx / 2,
    "z_domain%end": Wx / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": Nx - 1,
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
    # Slip cross-section walls; isothermal gradient (y) walls set below
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: mu(T) viscosity + sigma(T) + bulk conduction + independent temperature scalar
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T0,
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
    # Fluid 1 -- silicon oil (bulk): Arrhenius mu(T), real k / cp
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf_b / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b_ref,
    "fluid_pp(1)%k_therm": k_b,
    "fluid_pp(1)%visc_model": 1,
    "fluid_pp(1)%visc_c": C_b,
    "fluid_pp(1)%visc_d": D_b,
    # Fluid 2 -- Fluorinert (drop): Arrhenius mu(T), real k / cp
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d_ref,
    "fluid_pp(2)%k_therm": k_d,
    "fluid_pp(2)%visc_model": 1,
    "fluid_pp(2)%visc_c": C_d,
    "fluid_pp(2)%visc_d": D_d,
    # Patch 1 -- silicon oil filling the cell (3D cuboid)
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": Wx,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_b,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_d,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%T_temp_val": T_field,
    # Patch 2 -- Fluorinert drop (3D sphere), distinct density / properties from the bulk
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_d,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_field,
    # Isothermal Dirichlet gradient walls pin the cold floor / hot ceiling
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_c,
    "bc_y%Twall_out": T_h,
}

print(json.dumps(data))
