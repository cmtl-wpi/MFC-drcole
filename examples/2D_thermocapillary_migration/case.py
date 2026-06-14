#!/usr/bin/env python3
# Thermocapillary droplet migration -- MFC validation against Samareh, Mostaghimi & Moreau,
# "Thermocapillary migration of a deformable droplet", Int. J. Heat Mass Transfer 73 (2014) 616-626.
#
# WHAT THIS REPRODUCES
# Samareh's first validation case (their Section 4.1.1, "Motion of a drop in the limit of zero
# Marangoni number"). A neutrally-buoyant drop of diameter D sits in an imposed LINEAR temperature
# field. The surface tension falls with temperature,
#       sigma(T) = sigma0 + (dsigma/dT) * (T - T_ref),   dsigma/dT < 0,
# so the hot side of the interface has lower sigma. The resulting tangential (Marangoni) stress
# drags interfacial fluid hot->cold and, by reaction, the drop migrates toward the HOT side.
# In the zero-gravity, creeping-flow, ZERO-Marangoni-number limit Young, Goldstein & Block (1959)
# [Samareh Eq. 29] give the terminal speed
#       v_YGB = |sigma_T * gradT| * D / (6*mu_b + 9*mu_d).
# For equal viscosities (mu_d = mu_b = mu) this is v_YGB = (2/15) * |sigma_T| * gradT * (D/2) / mu.
#
# HOW MFC REALIZES THE "Ma = 0" PREMISE (updated: bulk conduction now available)
# Samareh's Test Case 1 holds the temperature field INVARIANT (their word): both fluids get an
# effectively INFINITE thermal diffusivity, so the imposed linear T never responds to the flow
# (Ma = 0). MFC now provides bulk Fourier conduction (the thermal_conduction feature: an explicit
# -k grad(T) energy flux with the harmonic mixture closure 1/k = sum(alpha_i/k_i), Samareh Eq. 8),
# selected here through the thermal Marangoni number knob SAMAREH_MA = U_r*r/alpha_T:
#   SAMAREH_MA = 0    : legacy frozen-T limit (alpha = 0, thermal Peclet -> infinity). The linear
#                       T is an initial condition that the flow slowly advects/distorts; ratios
#                       are only meaningful in the early quasi-steady window, and the 3D runs
#                       drift unboundedly past v_YGB at long times (README Sec. "Long-time
#                       behaviour", CONDUCTION_SCOPE.md).
#   SAMAREH_MA > 0    : finite-Ma physics. Conduction restores the imposed gradient on the
#                       timescale r^2/alpha = Ma*tau*(mu/(rho*U_r*r)); at the default Ma = 0.3
#                       that is roughly the viscous time, so the migration velocity has a true
#                       steady state to converge to (the frozen-T runaway disappears). Emulating
#                       the invariant-T premise outright (Ma -> 0, alpha -> inf) is the EXPENSIVE
#                       direction (explicit diffusion dt ~ Ma) -- see CONDUCTION_SCOPE.md Part I.
# The Ma = 0 anchors below remain exact only in the invariant-T limit; at Ma = 0.3 the measured
# ratio sits slightly below them (the gradient is mildly distorted near the interface).
#
# SAMAREH PARAMETERS (their Sec. 4.1.1) AND THE ONE DEVIATION
#   D = 1, box 5D wide x 7.5D long (3D: 5D x 5D x 7.5D), mu_d = mu_b = 0.1, sigma0 = 0.1,
#   sigma_T = dsigma/dT = -0.1, |gradT| = 0.1333 (= 2/15)  ->  v_YGB = 8.889e-3.
# These are matched exactly EXCEPT the density: Samareh use rho_d = rho_b = 0.2, while here
# rho(drop) = 1.0 (set by the EOS realization below). v_YGB and the terminal ratio are
# rho-independent in creeping flow, but the viscous time tau = rho*r^2/mu is 5x longer, so time
# axes are NOT comparable to the paper's. The compressible-EOS knobs (pi_inf, cv, p0) are
# MFC-specific, chosen for a stable low-Mach state. The migration axis is +x (the long, 7.5D
# axis); with no gravity the orientation is arbitrary, and an x-gradient keeps measure.py simple.
#
# GEOMETRY MODES (the validation needs both)
#   default          : drop CENTERED, OPEN (-3) boundaries -> approximates the UNBOUNDED domain;
#                      the 2D anchor is the unbounded-cylinder analytic 15/16*v_YGB (README).
#                      measure.py subtracts the small open-box return drift.
#   SAMAREH_WALL=1   : Samareh's actual Test Case 1 geometry -- SLIP (-2) walls on ALL boundaries,
#                      drop 1.5D from the cold (x.beg) wall -> compares against THEIR ~0.8 (2D).
#                      The closed box defines the rest frame (lab-frame velocity, no drift corr.).
#
# PARAMETERIZATION (env vars, so one build serves the whole validation sweep)
#   SAMAREH_DIM   : 2 or 3        (default 2)   -- 2D circle / 3D sphere
#   SAMAREH_NX    : cells per box WIDTH (the 5D short axis; Samareh used 64,128,256)  (default 128)
#   SAMAREH_DSDT  : dsigma/dT slope (default -0.1 = Samareh)  -- the Marangoni-strength sweep knob
#   SAMAREH_TAU   : run length in viscous times tau = rho*r^2/mu (default 3)
#   SAMAREH_WALL  : 1 = Samareh slip-wall geometry, 0 = open/centered (default 0)
#   SAMAREH_TS    : 1 = carry T as an independent scalar at UNIFORM density (decoupled; fixes the
#                   2D conduction reversal), 0 = legacy density-proxy temperature (default 0)
# The analytic density IC depends only on (T0, gradT), which are FIXED across the sweep, so its
# compiled-in Fortran never changes: build once, then run every variant with --no-build.
#
# GOTCHA (analytic-IC parser): MFC expands a bare `e` to Euler's number -- even inside a literal
# like 1e-9. Keep `e`-notation OUT of every analytic patch string (see rho_expr below).

import json
import os

# -- Variant selection (env vars; defaults = the 2D headline case at Samareh's medium grid) --
dim = int(os.environ.get("SAMAREH_DIM", "2"))  # 2 or 3
width_cells = int(os.environ.get("SAMAREH_NX", "128"))  # cells across the 5D box width
dsigma_dT = float(os.environ.get("SAMAREH_DSDT", "-0.1"))  # dsigma/dT (Samareh: -0.1)
n_tau = float(os.environ.get("SAMAREH_TAU", "3"))  # run length in viscous times tau (default 3)
wall = os.environ.get("SAMAREH_WALL", "0") == "1"  # Samareh slip-wall geometry vs open/centered
# Thermal Marangoni number Ma = U_r*r/alpha_T selecting the bulk-conduction strength
# (thermal_conduction feature, k* = 1). Ma = 0 (DEFAULT) disables conduction -> the validated
# frozen-T sigma(T) mode. Ma > 0 enables bulk conduction with a Dirichlet far-field temperature
# BC. CAVEAT: the conduction flux + BC are verified correct in isolation (see
# verify_1d_conduction.py: a static linear-T state stays static to ~3e-8), but on THIS example's
# density-gradient-imposed temperature the conduction couples with the moving droplet to reverse
# the apparent migration in 2D (correct toward hot for ~0.5 tau, then back). The density-proxy
# temperature makes this a poor conduction-validation vehicle; conduction mode is exploratory.
# The 3D OPEN box without the BC does tame the frozen-T runaway (run with SAMAREH_MA>0, SAMAREH_DIM=3
# and remove the isothermal block below to reproduce). See CONDUCTION_SCOPE.md.
Ma_th = float(os.environ.get("SAMAREH_MA", "0"))
# SAMAREH_TS = 1 carries temperature as an INDEPENDENT advected scalar (the thermal_scalar feature)
# instead of as a density proxy. The density is then held UNIFORM and T(x) is imposed directly on
# the scalar via patch_icpp%T_temp_val, so temperature is no longer slaved to density (1D check:
# verify_1d_thermal_scalar.py -- the scalar diffuses at exactly alpha*kappa^2 with velocity at
# machine zero). With SAMAREH_MA > 0 the scalar additionally diffuses at alpha = k/(rho cp); with
# SAMAREH_MA = 0 it is purely advected (genuinely frozen, but decoupled).
#
# IMPORTANT, measured result: decoupling does NOT remove the 2D conduction reversal. SAMAREH_TS=1
# SAMAREH_MA=0 (advection only) migrates correctly toward hot at +0.78 v_YGB, confirming the scalar
# advection and sigma(T_s) coupling are right and that the density proxy was NOT the cause. But
# SAMAREH_TS=1 SAMAREH_MA=0.3 (conduction + isothermal far-field BC) still reverses to ~-2.7 v_YGB,
# essentially identical to the density-proxy mode and to its incompressible-limit extrapolation. The
# reversal is therefore an artifact of the conduction + imposed-temperature-BC coupling with the
# migrating drop, not of the density proxy (revising the earlier diagnosis). Use SAMAREH_TS=1 with
# SAMAREH_MA=0 for a clean, density-decoupled frozen-T run; the finite-Ma conduction reversal is an
# open issue (see README Sec. 5 and THERMAL_SCALAR_SCOPE.md).
ts_mode = os.environ.get("SAMAREH_TS", "0") == "1"
# SAMAREH_NOBC = 1 (diagnostic) runs conduction with ADIABATIC x-boundaries (no isothermal far-field
# Dirichlet BC). Isolates whether the 2D conduction reversal comes from the imposed-temperature BC or
# from the bulk conduction operator itself. Adiabatic ends let the imposed gradient slowly relax.
no_bc = os.environ.get("SAMAREH_NOBC", "0") == "1"
assert dim in (2, 3), "SAMAREH_DIM must be 2 or 3"
assert Ma_th >= 0, "SAMAREH_MA must be >= 0 (0 disables conduction)"

# -- Geometry (Samareh Sec. 4.1.1: D=1 drop, 5D wide x 7.5D long box) --
D = 1.0  # droplet diameter
r = D / 2.0  # droplet radius = 0.5
W = 5.0 * D  # short-axis extent (box width), 2.5D clearance each side (= Samareh)
Lx = 7.5 * D  # long-axis extent (the gradient / migration axis, +x)
x_drop = (-Lx / 2 + 1.5 * D) if wall else 0.0  # Samareh: drop 1.5D from the cold wall; open: centered
bc = -2 if wall else -3  # slip (reflective) walls for the Samareh geometry, else open (extrapolation)

dx = W / width_cells  # isotropic cell size set by the box-width resolution
long_cells = round(Lx / dx)  # cells along the 7.5D migration axis (= 1.5*width_cells)
m = long_cells - 1  # MFC index = cells - 1
n = width_cells - 1
p = (width_cells - 1) if dim == 3 else 0

# Diffuse color interface ~2 cells wide at every resolution (sharp-interface limit as dx->0).
# smooth_coeff = dx/w with w = 2*dx  ->  smooth_coeff = 0.5 (constant), per the smoothing kernel
# eta = 0.5*(1 - tanh((smooth_coeff/min(dx..))*(|r_vec| - radius))).
cf_smooth_coeff = 0.5

# -- Equation of state (two IDENTICAL stiffened-gas fluids, gamma = 2; mu* = 1, k* = 1) --
# Background pressure kept well above the Laplace jump sigma0/r = 0.2 so the uniform-pressure IC is
# only mildly out of capillary equilibrium (low-Mach, stable). The flow is deeply incompressible
# (Mach ~ v_YGB/c0 ~ 1e-3), so the sound speed -- set by (p_inf, cv, p0) -- affects only the acoustic
# CFL, NOT v_YGB or the migration ratio. 3D is acoustically stiff and expensive, so it uses a SOFTER
# stiffened gas (lower c0 ~ 5 -> ~3x larger dt -> ~3x fewer steps). rho_coeff (=10) and T0 (=10) are
# held fixed in both, so the analytic density IC is byte-identical (no recompile between 2D and 3D)
# and the capillary perturbation stays small (rho'/rho = (2*sigma0/r)/(rho*c0^2) ~ 1.6% in 3D).
gam = 2.0
if dim == 3:
    p_inf, p0, cv = 10.0, 2.5, 1.25  # c0 ~ 5.0  (softer; same rho_coeff=10)
else:
    p_inf, p0, cv = 100.0, 25.0, 12.5  # c0 ~ 15.8 (the proven 2D example values)
mu = 0.1  # dynamic viscosity of both phases (Samareh: 0.1); MFC takes Re = 1/mu

# -- Imposed linear temperature field T(x) = T0 + gradT*x (centered on the drop at x=0) --
# Stiffened-gas EOS at uniform pressure p0:  T = (p0 + p_inf)/((gam-1)*rho*cv)
#   => rho(x) = (p0 + p_inf)/((gam-1)*cv*T(x)) = rho_coeff / (T0 + gradT*x).
# T0 = 10 keeps T (hence rho and sigma(T)) comfortably positive across the domain.
T0 = 10.0  # temperature at the domain center (x = 0)
gradT = 2.0 / 15.0  # |dT/dx| = 0.13333, Samareh's imposed gradient
sigma0 = 0.1  # surface tension at T_ref (Samareh)
T_ref = T0 + gradT * x_drop  # sigma = sigma0 AT THE DROP (drop sits at x_drop; = T0 in open mode)
rho_coeff = (p0 + p_inf) / ((gam - 1.0) * cv)  # = 10.0 ; rho(center) = rho_coeff/T0 = 1.0

eps = 1.0e-9  # trace volume fraction of the (identical) second phase
# Precompute (1-eps)*rho_coeff as a plain decimal: embedding eps=1e-9 in the analytic string would
# render "1e-09" and MFC would expand the `e` to Euler's number, corrupting rho to a NEGATIVE value.
rho_num = (1.0 - eps) * rho_coeff  # = 9.999999990
GRAD = "x"  # gradient / migration axis
if ts_mode:
    # Decoupled mode: hold density UNIFORM (= rho_coeff/T0 = 1.0) and impose T(x) directly on the
    # independent scalar T_s. Temperature is no longer a density proxy, so conduction cannot drive
    # the density inversion that reverses the migration. rho_expr must be a numeric value (not a
    # string) so it is written to the namelist unquoted -- a quoted numeric string would be read
    # into the real alpha_rho as a datatype mismatch.
    rho_expr = (1.0 - eps) * rho_coeff / T0  # uniform background density (plain float = 0.999999999)
    T_expr = f"{T0} + {gradT:.9f}*{GRAD}"  # T_s(x) = T0 + gradT*x imposed on the advected scalar
else:
    rho_expr = f"{rho_num:.9f}/({T0} + {gradT:.9f}*{GRAD})"  # ~ 10.0/(10.0 + 0.133333333*x)

# -- Young-Goldstein-Block terminal speed (mu* = 1, k* = 1) and the diagnostic non-dim numbers --
v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu  # = 8.88e-3 at dsigma_dT = -0.1
c0 = (gam * (p0 + p_inf)) ** 0.5  # reference sound speed: ~15.8 (2D) / ~5.0 (3D, softer EOS)
tau = (rho_coeff / T0) * r**2 / mu  # viscous time rho*r^2/mu = 2.5

# -- Bulk thermal conduction (k* = 1): alpha_T from the requested thermal Marangoni number --
# U_r = |sigma_T|*gradT*r/mu is the Marangoni interfacial velocity scale; alpha_T = U_r*r/Ma_th.
# k follows from alpha_T = k/(rho*cv*gam) (stiffened-gas cp = gam*cv) at the drop-center density
# rho_ref = rho_coeff/T0 = 1. Both fluids get the same k (the YGB k* = 1 limit).
if Ma_th > 0:
    U_r = (-dsigma_dT) * gradT * r / mu
    alpha_T = U_r * r / Ma_th
    k_therm = alpha_T * (rho_coeff / T0) * cv * gam

# -- Time stepping: acoustic-CFL limited (Ma ~ U/c0 ~ 1e-3, so cost is set by step count) --
# dt = ICFL*dx/c0 with ICFL = 0.40 (RK3+WENO5 stable; matches the proven 2D example's dt ~ 0.025*dx
# at c0 ~ 15.8). The softer 3D EOS (c0 ~ 5) automatically yields ~3x larger dt here.
# With conduction on, the explicit diffusion number is additionally capped at 0.35
# (dt <= 0.35*dx^2/(2*dim*alpha_T)). At the default Ma_th = 0.3 the cap is slack on the
# coarse 2D grids but binds ~1.6x on 2D_w256 and ~2-3x on the 3D grids (diffusion dt
# scales with dx^2 while the acoustic dt scales with dx) -- the modest price of finite-Ma
# physics; emulating Ma -> 0 outright would cost another order of magnitude.
# Run 3 viscous times tau to reach the quasi-steady migration plateau.
mydt = 0.40 * dx / c0
if Ma_th > 0:
    mydt = min(mydt, 0.35 * dx**2 / (2.0 * dim * alpha_T))
t_end = n_tau * tau  # default 3*tau = 7.5; raise SAMAREH_TAU to test long-time plateau vs drift
t_step_stop = int(round(t_end / mydt))
t_step_save = max(1, t_step_stop // 50)  # ~50 snapshots

data = {
    # Logistics
    "run_time_info": "T",
    # Computational domain: long (gradient) axis x in [-Lx/2, Lx/2]; short axes in [-W/2, W/2]
    "x_domain%beg": -Lx / 2,
    "x_domain%end": Lx / 2,
    "y_domain%beg": -W / 2,
    "y_domain%end": W / 2,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation algorithm (6-equation model; same proven WENO5/HLLC settings as the 2D example)
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
    # Boundaries: open (-3, approximates the unbounded domain) or slip walls (-2, Samareh's box)
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
    "patch_icpp(1)%geometry": 9 if dim == 3 else 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": W,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    # Patch 2 -- droplet (circle in 2D / sphere in 3D): marks c=1, smeared over ~2 cells. IDENTICAL
    # density/composition/pressure to patch 1, so the capillary stress acts purely on the c interface
    # with no real fluid/density jump (the mu*=1, k*=1 / undistorted-T limit of YGB).
    "patch_icpp(2)%geometry": 8 if dim == 3 else 2,
    "patch_icpp(2)%x_centroid": x_drop,
    "patch_icpp(2)%y_centroid": 0.0,
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

if dim == 3:
    # Third short axis z in [-W/2, W/2]; extend velocity/centroid/length to 3D.
    data.update(
        {
            "z_domain%beg": -W / 2,
            "z_domain%end": W / 2,
            "bc_z%beg": bc,
            "bc_z%end": bc,
            "patch_icpp(1)%z_centroid": 0.0,
            "patch_icpp(1)%length_z": W,
            "patch_icpp(1)%vel(3)": 0.0,
            "patch_icpp(2)%z_centroid": 0.0,
            "patch_icpp(2)%vel(3)": 0.0,
        }
    )

if Ma_th > 0:
    # Bulk Fourier conduction (harmonic mixture closure; equal k -> the YGB k* = 1 limit).
    data.update(
        {
            "thermal_conduction": "T",
            "fluid_pp(1)%k_therm": k_therm,
            "fluid_pp(2)%k_therm": k_therm,
        }
    )
    # Isothermal far-field BC: pin the x-end temperatures to T(x) = T0 + gradT*x via the Dirichlet
    # reflection T_ghost = 2*Twall - T_interior. This is a WALL boundary condition and is ONLY applied
    # with the slip-wall geometry (SAMAREH_WALL=1). Measured fact: applying it at an OPEN boundary
    # (bc_x = -3, the default) catastrophically REVERSES the drop (-2.7 v_YGB) -- the fixed-value
    # reflection fights the advective throughflow and pumps the thermal field. The open box instead
    # leaves the x-ends ADIABATIC, which migrates correctly toward hot (+0.75 v_YGB); the imposed
    # gradient is sustained by the IC over the quasi-steady window (conduction relaxes it only on the
    # slow domain^2/alpha timescale). See README Sec. 5 and the SAMAREH_NOBC diagnostic.
    if wall and not no_bc:
        data.update(
            {
                "bc_x%isothermal_in": "T",
                "bc_x%isothermal_out": "T",
                "bc_x%Twall_in": T0 + gradT * (-Lx / 2.0),  # cold (-x) wall temperature
                "bc_x%Twall_out": T0 + gradT * (Lx / 2.0),  # hot (+x) wall temperature
            }
        )

if ts_mode:
    # Independent temperature scalar: T(x) imposed directly on T_s for BOTH patches (temperature is
    # continuous across the color interface; only the color function jumps). sigma(T) reads T_s, and
    # with SAMAREH_MA > 0 the thermal_conduction block above diffuses T_s at alpha = k/(rho cp).
    # T_s_wrt exposes the field for visualization (finally a real temperature output).
    data.update(
        {
            "thermal_scalar": "T",
            "T_s_wrt": "T",
            "patch_icpp(1)%T_temp_val": T_expr,
            "patch_icpp(2)%T_temp_val": T_expr,
        }
    )

print(json.dumps(data))
