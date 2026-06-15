#!/usr/bin/env python3
# Thermocapillary droplet migration -- MFC validation against Samareh, Mostaghimi & Moreau,
# "Thermocapillary migration of a deformable droplet", Int. J. Heat Mass Transfer 73 (2014) 616-626.
#
# WHAT THIS REPRODUCES
# Samareh's first validation case in 2D (their Section 4.1.1 / Fig 5, "Motion of a drop in the
# limit of zero Marangoni number"). A neutrally-buoyant drop of diameter D sits 1.5D above the
# COLD bottom wall of a slip-walled box, in an imposed LINEAR vertical temperature field
# (T = 0 on the cold floor, increasing upward). The surface tension falls with temperature,
#       sigma(T) = sigma0 + (dsigma/dT) * (T - T_ref),   dsigma/dT < 0,
# so the hot (upper) side of the interface has lower sigma. The resulting tangential (Marangoni)
# stress drags interfacial fluid hot->cold and, by reaction, the drop RISES toward the hot top.
# In the zero-gravity, creeping-flow, ZERO-Marangoni-number limit Young, Goldstein & Block (1959)
# [Samareh Eq. 29] give the terminal rise speed
#       v_YGB = |sigma_T * gradT| * D / (6*mu_b + 9*mu_d).
# For equal viscosities (mu_d = mu_b = mu) this is v_YGB = (2/15) * |sigma_T| * gradT * (D/2) / mu.
# Samareh report the converged 2D ratio v_t/v_YGB ~ 0.80 (Fig 5).
#
# THIS EXAMPLE IS 2D ONLY. Samareh's 3D companion (Fig 6, ~0.95) is deliberately NOT here: it
# belongs in a separate 3D_thermocapillary_migration example, and -- more to the point -- on this
# no-conduction branch the 3D thermocapillary rise has NO validatable result. The frozen-T field is
# advected in all three directions and the 3D rise velocity drifts UNBOUNDEDLY (finer grid -> faster;
# see mfc-sigmaT-3d-drift-no-conduction), so there is no quasi-steady plateau to compare against the
# 0.95. Adding 3D would only manufacture a number off a curve that never plateaus.
#
# HOW MFC REALIZES THE "Ma = 0" PREMISE (and where that breaks down)
# Samareh's Test Case 1 holds the temperature field INVARIANT (their word): both fluids get an
# effectively INFINITE thermal diffusivity, so the imposed linear T never responds to the flow
# (Ma = 0). In the DEFAULT mode here (SAMAREH_MA = 0) MFC uses the opposite diffusivity limit
# (alpha -> 0, thermal Peclet -> infinity): the linear T is a frozen initial condition that the flow
# slowly advects/distorts. The two limits agree at EARLY times, before interfacial parcels (speed
# ~ U_r = |sigma_T|*gradT*r/mu) have traveled far enough to reshape the frozen gradient. All
# quoted ratios therefore come from a stated measurement window; the long-time diagnostics measure
# where the frozen-IC approximation holds (2D: a true plateau) and where it fails (3D: an unbounded
# velocity drift). Bulk Fourier conduction is now available (SAMAREH_MA > 0, the thermal_conduction
# feature) for finite-Ma physics, and an independent temperature scalar (SAMAREH_TS) decouples T from
# density; both are documented in the env-var list below. The finite-Ma cases (Samareh Sec. 4.1.2 /
# 4.2, Figs 7/8/12/13/16) need bulk conduction (now available) and temperature-dependent viscosity
# mu(T) (which MFC still does not provide), so a fully faithful finite-Ma reproduction is not yet
# possible here; SAMAREH_MA enables the conduction half of that physics.
#
# SAMAREH PARAMETERS (their Sec. 4.1.1) AND THE ONE DEVIATION
#   D = 1, box 5D wide x 7.5D tall, rho_d = rho_b = 0.2, mu_d = mu_b = 0.1,
#   sigma0 = 0.1, sigma_T = dsigma/dT = -0.1, |gradT| = 0.1333 (= 2/15)  ->  v_YGB = 8.889e-3.
# These are matched, INCLUDING the density rho = 0.2 (set by the EOS realization below), so the
# viscous time tau = rho*r^2/mu = 0.5 is Samareh's and the time axes are comparable. The ONE
# deviation is the absolute temperature baseline: Samareh use T = 0 on the cold wall, but the
# density-proxy IC below (rho = rho_coeff/T) would diverge as T -> 0, so T is shifted up by T0 = 10
# (T spans ~9.5..10.5 across the box). Only the absolute T level changes; the gradient |gradT| and
# the slope sigma_T -- the only things v_YGB and the Marangoni stress depend on -- are exact. The
# compressible-EOS knobs (pi_inf, cv, p0) are MFC-specific, chosen for a stable low-Mach state.
#
# GEOMETRY MODES
#   SAMAREH_WALL=1 (DEFAULT): Samareh's actual Fig 5/6 geometry -- SLIP (-2) walls on ALL sides,
#                  drop 1.5D above the cold (y.beg) wall. The closed box defines the rest frame
#                  (lab-frame velocity, no drift correction) -> compares against THEIR 0.80 / 0.95.
#   SAMAREH_WALL=0: drop CENTERED, OPEN (-3) boundaries -> approximates the UNBOUNDED domain; the
#                  2D anchor is then the unbounded-cylinder analytic 15/16*v_YGB; measure.py
#                  subtracts the small open-box return drift.
#
# PARAMETERIZATION (env vars, so one build serves the whole validation sweep)
#   SAMAREH_NX    : cells per box WIDTH (the 5D short axis; Samareh used 64,128,256)  (default 128)
#   SAMAREH_DSDT  : dsigma/dT slope (default -0.1 = Samareh)  -- the Marangoni-strength sweep knob
#   SAMAREH_TR    : run length in capillary-thermal times t_r = mu/|sigma_T*gradT| = 7.5 (default 4)
#   SAMAREH_WALL  : 1 = Samareh slip-wall geometry (default), 0 = open/centered
#   SAMAREH_MA    : thermal Marangoni number; > 0 enables bulk Fourier conduction (default 0 = frozen-T).
#                   k follows from alpha = U_r*r/Ma. The isothermal Dirichlet wall BC (cold floor / hot
#                   ceiling) is imposed with the slip-wall geometry; the open box stays adiabatic (a
#                   Dirichlet wall BC belongs at a closed wall, not an outflow boundary). Wall +
#                   conduction gives a stable forward plateau v_t/v_YGB ~ 0.85 (Samareh ~0.80).
#   SAMAREH_TS    : 1 = carry T as an INDEPENDENT advected scalar (thermal_scalar) at UNIFORM density,
#                   decoupled from the EOS; T(y) imposed via patch_icpp%T_temp_val, output via T_s_wrt.
# The analytic density IC depends only on (T0, gradT), which are FIXED across the sweep, so its
# compiled-in Fortran never changes: build once, then run every variant with --no-build.
#
# GOTCHA (analytic-IC parser): MFC expands a bare `e` to Euler's number -- even inside a literal
# like 1e-9. Keep `e`-notation OUT of every analytic patch string (see rho_expr below).

import json
import os

# -- Variant selection (env vars; defaults = the 2D headline case at Samareh's medium grid) --
width_cells = int(os.environ.get("SAMAREH_NX", "128"))  # cells across the 5D box width
dsigma_dT = float(os.environ.get("SAMAREH_DSDT", "-0.1"))  # dsigma/dT (Samareh: -0.1)
n_tr = float(os.environ.get("SAMAREH_TR", "4"))  # run length in capillary-thermal times t_r
wall = os.environ.get("SAMAREH_WALL", "1") == "1"  # Samareh slip-wall geometry (default) vs open
# Finite-Ma conduction / independent-temperature modes (ported from the thermal_scalar feature):
Ma_th = float(os.environ.get("SAMAREH_MA", "0"))  # thermal Marangoni number; > 0 enables bulk Fourier conduction
ts_mode = os.environ.get("SAMAREH_TS", "0") == "1"  # carry T as an independent advected scalar (decoupled from rho)
assert Ma_th >= 0, "SAMAREH_MA must be >= 0 (0 disables conduction)"

# -- Geometry (Samareh Sec. 4.1.1: D=1 drop, 5D wide x 7.5D tall box; rise axis = y) --
D = 1.0  # droplet diameter
r = D / 2.0  # droplet radius = 0.5
W = 5.0 * D  # short-axis extent (box width x), 2.5D clearance each side (= Samareh)
Ly = 7.5 * D  # long-axis extent (the gradient / RISE axis, +y)
y_drop = (-Ly / 2 + 1.5 * D) if wall else 0.0  # Samareh: drop 1.5D above the cold floor; open: centered
bc = -2 if wall else -3  # slip (reflective) walls for the Samareh geometry, else open (extrapolation)

dx = W / width_cells  # isotropic cell size set by the box-width resolution
long_cells = round(Ly / dx)  # cells along the 7.5D rise axis (= 1.5*width_cells)
m = width_cells - 1  # x (short) -- MFC index = cells - 1
n = long_cells - 1  # y (long, rise axis)
p = 0  # 2D

# Diffuse color interface ~2 cells wide at every resolution (sharp-interface limit as dx->0).
# smooth_coeff = dx/w with w = 2*dx  ->  smooth_coeff = 0.5 (constant), per the smoothing kernel
# eta = 0.5*(1 - tanh((smooth_coeff/min(dx..))*(|r_vec| - radius))).
cf_smooth_coeff = 0.5

# -- Equation of state (two IDENTICAL stiffened-gas fluids, gamma = 2; mu* = 1, k* = 1) --
# Background pressure kept well above the Laplace jump sigma0/r = 0.2 so the uniform-pressure IC is
# only mildly out of capillary equilibrium (low-Mach, stable). The flow is deeply incompressible
# (Mach ~ v_YGB/c ~ 4e-4), so the sound speed -- set by (p_inf, cv, p0) -- affects only the acoustic
# CFL, NOT v_YGB or the migration ratio. cv is chosen so rho(center) = rho_coeff/T0 = 0.2 (Samareh).
# c = sqrt(gam*(p+p_inf)/rho) ~ 20 at rho = 0.2: a moderate stiffness that keeps the capillary
# perturbation small (rho'/rho = (sigma0/r)/(rho*c^2) ~ 0.25%) without an over-stiff (expensive) dt.
gam = 2.0
p_inf, p0 = 32.0, 8.0  # with rho ~ 0.2: c = sqrt(gam*(p0+p_inf)/rho) ~ 20
mu = 0.1  # dynamic viscosity of both phases (Samareh: 0.1); MFC takes Re = 1/mu

# -- Imposed linear temperature field T(y) = T0 + gradT*y (centered on the box at y=0) --
# Stiffened-gas EOS at uniform pressure p0:  T = (p0 + p_inf)/((gam-1)*rho*cv)
#   => rho(y) = (p0 + p_inf)/((gam-1)*cv*T(y)) = rho_coeff / (T0 + gradT*y).
# T0 = 10 keeps T (hence rho and sigma(T)) comfortably positive across the domain.
T0 = 10.0  # temperature at the box center (y = 0); shifted up from Samareh's T(cold wall)=0
gradT = 2.0 / 15.0  # |dT/dy| = 0.13333, Samareh's imposed gradient
sigma0 = 0.1  # surface tension at T_ref (Samareh)
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2
rho_coeff = rho_drop * T0  # = 2.0 ; rho(center) = rho_coeff/T0 = 0.2
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # = 62.5 ; closes the EOS so rho(center) = 0.2
T_ref = T0 + gradT * y_drop  # sigma = sigma0 AT THE DROP (drop sits at y_drop; = T0 in open mode)

eps = 1.0e-9  # trace volume fraction of the (identical) second phase
# Precompute (1-eps)*rho_coeff as a plain decimal: embedding eps=1e-9 in the analytic string would
# render "1e-09" and MFC would expand the `e` to Euler's number, corrupting rho to a NEGATIVE value.
rho_num = (1.0 - eps) * rho_coeff  # = 1.999999998
GRAD = "y"  # gradient / rise axis
if ts_mode:
    # Decoupled mode: hold density UNIFORM (= rho_coeff/T0 = rho_drop = 0.2) and impose T(y) directly
    # on the independent scalar T_s. A numeric value (not a string) so it is written unquoted to the
    # namelist -- a quoted numeric string would be read into the real alpha_rho as a datatype mismatch.
    rho_expr = (1.0 - eps) * rho_coeff / T0  # uniform background density (= 0.2)
    T_expr = f"{T0} + {gradT:.9f}*{GRAD}"  # T_s(y) = T0 + gradT*y imposed on the advected scalar
else:
    rho_expr = f"{rho_num:.9f}/({T0} + {gradT:.9f}*{GRAD})"  # ~ 2.0/(10.0 + 0.133333333*y)

# -- Young-Goldstein-Block terminal speed (mu* = 1, k* = 1) and the diagnostic time scales --
v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu  # = 8.88e-3 at dsigma_dT = -0.1
# Max sound speed over the domain (stiffened gas: c = sqrt(gam*(p+p_inf)/rho)); the lowest density
# is at the hot top (highest T), so that sets the acoustic CFL. NOTE: c scales as 1/sqrt(rho), so
# the rho=0.2 here makes c ~ sqrt(5) larger than the naive sqrt(gam*(p0+p_inf)) -- the dt must use
# this, or ICFL exceeds 1 within ~20 steps.
rho_min = rho_coeff / (T0 + gradT * (Ly / 2.0))  # density at the hot top wall (max sound speed)
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5  # ~ 20.5
tau = rho_drop * r**2 / mu  # viscous time rho*r^2/mu = 0.5 (Samareh's)
t_r = mu / abs(dsigma_dT * gradT)  # capillary-thermal time mu/|sigma_T*gradT| = 7.5 (Samareh scale)

# -- Bulk thermal conduction (k* = 1): alpha_T from the requested thermal Marangoni number Ma --
# U_r = |sigma_T|*gradT*r/mu is the Marangoni interfacial velocity scale; alpha_T = U_r*r/Ma. The
# conductivity follows from alpha_T = k/(rho*cv*gam) (stiffened-gas cp = gam*cv) at the drop density
# rho = rho_coeff/T0 = 0.2; both fluids get the same k (the YGB k* = 1 limit).
if Ma_th > 0:
    U_r = (-dsigma_dT) * gradT * r / mu
    alpha_T = U_r * r / Ma_th
    k_therm = alpha_T * (rho_coeff / T0) * cv * gam

# -- Time stepping: acoustic-CFL limited (Ma ~ U/c ~ 4e-4, so cost is set by step count) --
# dt = ICFL*dx/c_max with ICFL = 0.35 (RK3 + WENO5 stable, with margin for the brief acoustic
# capillary-relaxation transient that briefly raises |u|+c). Run n_tr capillary-thermal times t_r
# to reach -- and hold -- the quasi-steady migration plateau (Samareh's Fig 5/6 plateau by t/t_r ~ 1-2).
mydt = 0.35 * dx / c_max
# With conduction on, also cap by the explicit-diffusion number: dt <= 0.35*dx^2/(2*d*alpha), d=2 (2D).
if Ma_th > 0:
    mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_T))
t_end = n_tr * t_r  # default 4*t_r = 30; covers the plateau with margin
t_step_stop = int(round(t_end / mydt))
t_step_save = max(1, t_step_stop // 80)  # ~80 snapshots, enough to resolve the rise-velocity curve

data = {
    # Logistics
    "run_time_info": "T",
    # Computational domain: rise (gradient) axis y in [-Ly/2, Ly/2]; short axes in [-W/2, W/2]
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
    # Boundaries: slip walls (-2, Samareh's box, default) or open (-3, approximates unbounded domain)
    "bc_x%beg": bc,
    "bc_x%end": bc,
    "bc_y%beg": bc,
    "bc_y%end": bc,
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
    # Patch 1 -- background medium spanning the domain: analytic linear-T density, ambient color c=0.
    "patch_icpp(1)%geometry": 3,  # 2D rectangle
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
    # Patch 2 -- droplet (2D circle): marks c=1, smeared over ~2 cells. IDENTICAL density/composition/
    # pressure to patch 1, so the capillary stress acts purely on the c interface with no real
    # fluid/density jump (the mu*=1, k*=1 / undistorted-T limit of YGB).
    "patch_icpp(2)%geometry": 2,  # 2D circle
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": cf_smooth_coeff,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": rho_expr,
    "patch_icpp(2)%alpha_rho(2)": eps,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1.0,
}

if Ma_th > 0:
    # Bulk Fourier conduction (harmonic mixture closure; equal k -> the YGB k* = 1 limit).
    data.update({"thermal_conduction": "T", "fluid_pp(1)%k_therm": k_therm, "fluid_pp(2)%k_therm": k_therm})
    # Isothermal Dirichlet wall BC pins the cold floor / hot ceiling to the imposed gradient -- the
    # faithful Samareh closed-box setup. Applied with the slip-wall geometry only (a Dirichlet wall BC
    # does not belong at the open/outflow box, which stays adiabatic). Verified: wall + conduction
    # gives a stable forward plateau v_t/v_YGB ~ 0.85, matching Samareh's ~0.80. (This case used to
    # REVERSE -- that was an MPI halo-corruption bug, interior ranks overwriting their halo cells in
    # the isothermal conduction BC, since fixed in s_apply_thermal_conduction_bc; NOT the density
    # proxy or advective throughflow as earlier believed. See CONDUCTION_REVERSAL_SAGA.md.)
    if wall:
        data.update(
            {
                "bc_y%isothermal_in": "T",
                "bc_y%isothermal_out": "T",
                "bc_y%Twall_in": T0 + gradT * (-Ly / 2.0),  # cold floor (y%beg)
                "bc_y%Twall_out": T0 + gradT * (Ly / 2.0),  # hot ceiling (y%end)
            }
        )

if ts_mode:
    # Independent temperature scalar: T(y) imposed directly on T_s for both patches (T is continuous
    # across the color interface; only the color function jumps). sigma(T) reads T_s; with Ma > 0 the
    # conduction block above diffuses T_s at alpha = k/(rho cp). T_s_wrt exposes it for visualization.
    data.update(
        {
            "thermal_scalar": "T",
            "T_s_wrt": "T",
            "patch_icpp(1)%T_temp_val": T_expr,
            "patch_icpp(2)%T_temp_val": T_expr,
        }
    )

print(json.dumps(data))
