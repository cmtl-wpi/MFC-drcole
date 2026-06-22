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
# the experiment's non-monotonic rise-velocity loop.
#
# TEMPERATURE WITH DISTINCT FLUIDS (compressible code, dimensional / SI).
# MFC is compressible: T = (p+p_inf)/((gam-1)*rho*cv), so the imposed field T(y) (in Kelvin) is
# encoded by a PER-FLUID density proxy, rho_i(x,y,z) = (p + p_inf_i)/((gam-1)*cv_i*T(y)). The
# fluid-to-fluid density ratio is height-independent (T cancels), so the real Fluorinert/oil ratio
# (~1.88) is preserved; the mixture EOS recovers T_mix = T(y) exactly across the interface; and
# sigma(T)/mu(T) read it. Bulk conduction + isothermal gradient walls sustain the field.
#
# SOUND-SPEED SOFTENING. The real liquid sound speed (~750 m/s with the true cp) makes dt tiny. We
# soften it to c_snd for tractable dt; that softens cv = cp/gam, so we RESCALE each fluid's k to keep
# the physical thermal diffusivity alpha = k/(rho*cp) -- and Ma = U_r*r/alpha -- unchanged. The
# absolute cp and k are then soft-EOS artifacts; the migration numbers (Re, Ma, Ca, the property
# ratios, mu(T)) are physical. Honest caveat: each fluid's density stratifies as ~1/T (a
# compressibility artifact absent in the incompressible reference; ~ +-10% over the 60 K cell).
# A converged comparison is heavy (Samareh used 120 x 320 x 120 up to 240 x 640 x 240); this case
# runs a coarse demo grid by default. Keep `e`/`r`/`eps` tokens out of analytic strings.

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

# Imposed temperature field (Kelvin): T(y) = T0 + gradT*y, T0 at the cell center y = 0
T_c, T_h = 283.0, 343.0  # cold floor / hot ceiling wall temperatures
gradT = (T_h - T_c) / Ly  # = 1000 K/m
T0 = 313.0  # reference temperature (cell center); densities are quoted here
sigma0 = 0.007  # surface tension at T0 (N/m)
sigma_T = -3.6e-5  # dsigma/dT (N/m/K)

# Per-fluid material properties (fluid 1 = silicon oil bulk, fluid 2 = Fluorinert drop)
rho_b, rho_d = 918.3, 1727.7  # densities at T0 (kg/m^3)
k_b_real, k_d_real = 0.13389, 0.063  # thermal conductivities (W/mK)
cp_b_real, cp_d_real = 1778.2, 1047.0  # specific heats (J/kgK)
C_b, D_b = -10.17, 1643.0  # silicon oil Arrhenius coefficients mu = exp(C + D/T)
C_d, D_d = -11.76, 1540.0  # Fluorinert Arrhenius coefficients

# Stiffened-gas EOS with a softened sound speed so the density proxy stays cheap. p_inf_i and cv_i
# follow so c = c_snd and rho = rho_i at (p0, T0); cv comes out equal for both (= c_snd^2/(gam(gam-1)T0)).
gam = 2.0
c_snd = 30.0  # softened sound speed (m/s); migration Mach ~ 1e-3
p0 = 1.0e5  # background pressure (~1 atm)
p_inf_b = rho_b * c_snd**2 / gam - p0
p_inf_d = rho_d * c_snd**2 / gam - p0
cv_b = (p0 + p_inf_b) / ((gam - 1.0) * rho_b * T0)
cv_d = (p0 + p_inf_d) / ((gam - 1.0) * rho_d * T0)
cp_b, cp_d = gam * cv_b, gam * cv_d  # soft cp

# Rescale k to keep the PHYSICAL thermal diffusivity (hence Ma) with the soft cp
alpha_b = k_b_real / (rho_b * cp_b_real)  # physical diffusivity
alpha_d = k_d_real / (rho_d * cp_d_real)
k_b = alpha_b * rho_b * cp_b  # rescaled conductivity (soft EOS)
k_d = alpha_d * rho_d * cp_d

# Reference viscosity at the drop-center temperature at t = 0 (visc_model = 1 evaluates mu(T) live)
T_visc_ref = T0 + gradT * y_drop  # = 298 K
mu_b_ref = math.exp(C_b + D_b / T_visc_ref)
mu_d_ref = math.exp(C_d + D_d / T_visc_ref)

# Migration scales (Samareh): U_r = |sigma_T*gradT|*r/mu_b(T0), t_r = mu_b(T0)/|sigma_T*gradT|
mu_b0 = math.exp(C_b + D_b / T0)  # silicon oil viscosity at T0 (~7.3e-3 Pa.s)
G = abs(sigma_T * gradT)  # Marangoni stress scale
t_r = mu_b0 / G  # capillary-thermal time ~ 0.20 s

# Time stepping: acoustic CFL + 3D explicit-diffusion cap on the smallest cell
mydt = 0.35 * dx / c_snd
mydt = min(mydt, 0.35 * dx**2 / (6.0 * max(alpha_b, alpha_d)))
t_step_stop = round(1.0 * t_r / mydt)  # 1 capillary-thermal time
t_step_save = max(1, t_step_stop // 60)

# One analytic patch (a second patch's analytic density would leak globally). The smooth sphere
# eta(x,y,z) ~ 1 in the drop / 0 outside drives color, the volume-fraction split, the Laplace
# pressure jump, and BOTH per-fluid densities. Hardcode center/radius as decimals.
xc_d, yc_d, zc_d, r_d = 0.0, y_drop, 0.0, r
w_if = 0.75 * dx
laplace = sigma0 / r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2 + (z - ({zc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
Texpr = f"({T0:.6f} + {gradT:.6f}*y)"
pres_expr = f"({p0:.6f} + {laplace:.6f}*({eta}))"
alpha1_expr = f"(1.0 - ({eta}))"
alpha2_expr = f"({eta})"
cb = (gam - 1.0) * cv_b
cd = (gam - 1.0) * cv_d
arho1_expr = f"(1.0 - ({eta}))*({pres_expr} + {p_inf_b:.6f})/({cb:.9f}*{Texpr})"
arho2_expr = f"({eta})*({pres_expr} + {p_inf_d:.6f})/({cd:.9f}*{Texpr})"
cf_expr = f"({eta})"

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
    "mpp_lim": "T",
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
    "num_patches": 1,
    "num_fluids": 2,
    # Physics: mu(T) viscosity + sigma(T) + bulk conduction (T from the per-fluid density proxy)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T0,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "T_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 -- silicon oil (bulk): Arrhenius mu(T), rescaled k, soft cv
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf_b / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b_ref,
    "fluid_pp(1)%k_therm": k_b,
    "fluid_pp(1)%visc_model": 1,
    "fluid_pp(1)%visc_c": C_b,
    "fluid_pp(1)%visc_d": D_b,
    # Fluid 2 -- Fluorinert (drop): Arrhenius mu(T), rescaled k, soft cv; stiffening keeps rho_d/rho_b
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d_ref,
    "fluid_pp(2)%k_therm": k_d,
    "fluid_pp(2)%visc_model": 1,
    "fluid_pp(2)%visc_c": C_d,
    "fluid_pp(2)%visc_d": D_d,
    # Single analytic patch: drop lives entirely in eta(x,y,z). Color, volume fractions, pressure,
    # and both per-fluid densities share it, so the mixture EOS recovers the linear T(y) exactly
    # while the distinct Fluorinert/oil density ratio is preserved.
    "patch_icpp(1)%geometry": 9,  # 3D cuboid spanning the cell
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": Wx,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": arho1_expr,
    "patch_icpp(1)%alpha_rho(2)": arho2_expr,
    "patch_icpp(1)%alpha(1)": alpha1_expr,
    "patch_icpp(1)%alpha(2)": alpha2_expr,
    "patch_icpp(1)%cf_val": cf_expr,
    # Isothermal Dirichlet gradient walls pin the cold floor / hot ceiling
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_c,
    "bc_y%Twall_out": T_h,
}

print(json.dumps(data))
