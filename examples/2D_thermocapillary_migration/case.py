#!/usr/bin/env python3
# 2D thermocapillary (Marangoni) droplet migration -- validation of the temperature-dependent
# surface tension closure sigma(T) (sigma_model = 1).
#
# PHYSICS
# A circular "droplet" (a smooth blob of the diffuse-interface color function c) sits in a
# quiescent, mechanically-equilibrated medium carrying a linear temperature field T(x) = T0 + |gradT|*x.
# The surface tension decreases with temperature,
#       sigma(T) = sigma0 + (dsigma/dT) * (T - T_ref),   dsigma/dT < 0,
# so sigma is lower on the hot (+x) side of the interface. The resulting tangential (Marangoni)
# stress drives interfacial fluid from the hot side toward the cold side, and by momentum balance
# the droplet migrates toward the HOT side (+x). This force emerges automatically from the spatial
# variation of the capillary stress tensor -- no extra source term is added.
#
# REFERENCE (Young, Goldstein & Block 1959): in the zero-gravity, low-Reynolds (Stokes),
# low-Marangoni limit the terminal migration speed is
#       U_YGB = 2 / ((2 + k*)(2 + 3*mu*)) * ( -dsigma/dT * |gradT| * a / mu_c ),
# with a = droplet radius, mu_c = continuous-phase viscosity, k* = k_drop/k_cont,
# mu* = mu_drop/mu_cont. Here both phases are identical (mu* = 1) and the imposed temperature
# field is undistorted (the k* = 1 / uniform-conductivity limit), so the prefactor is 2/15:
#       U_YGB = (2/15) * (-dsigma/dT) * |gradT| * a / mu_c.
# With identical inside/outside properties the two-phase droplet is dynamically equivalent to this
# color-function blob, which keeps the initial condition clean (single uniform-composition patch).
#
# IMPORTANT CAVEAT (no bulk conduction): MFC has no Fourier heat conduction in the bulk energy
# equation, so the imposed linear T is FROZEN as an initial condition and is slowly advected/
# distorted by the developing flow. This setup therefore validates the Marangoni STRESS COUPLING
# (the new physics) in the quasi-steady window: measure the droplet's migration velocity on the
# plateau it reaches after ~a few viscous times tau = rho*a^2/mu and before T advects appreciably
# (low Peclet). It is not a full conjugate-conduction thermocapillary validation.
#
# MEASUREMENT: track the color-function (c = 0.5) centroid x-position over time (cf_wrt = T),
# finite-difference it for the migration velocity, and compare the quasi-steady plateau to U_YGB.
#
# STARTUP-STABILITY NOTE (read before running):
# The case validates and pre_process initializes it cleanly, but the explicit-compressible
# simulation currently trips the ICFL>1 (infinite sound speed) guard at ~step 2. This was traced
# to the *initial condition*, NOT the sigma(T) implementation: a control run with sigma_model = 0
# (constant sigma) -- and even with the surface-tension force set to ~0 -- blows up identically,
# and it is independent of model_eqns (2 vs 3), reconstruction (WENO5 vs MUSCL+int_comp),
# mpp_lim, and boundary type (-2/-3/-6). The common trigger is the analytic density gradient
# rho(x) = (p0+p_inf)/((gam-1)*cv*T(x)) used to impose T(x) at uniform pressure: although a density
# gradient at uniform pressure with u = 0 is a continuum equilibrium, this discrete initialization
# is not reproducing a consistent (p, rho, E) static state here. To stabilize before production:
#   - shrink the imposed density variation (raise T0 so rho varies by <~2% across the domain while
#     keeping dsigma/dT and |gradT|, hence U_YGB, fixed), and/or
#   - initialize from a relaxed/quiescent restart, or impose T via a balanced body force, and/or
#   - first run the simpler Test B (flat interface, imposed linear T, measure the tangential force)
#     to verify the Marangoni coupling before attempting full migration.
# The sigma(T) closure itself is exercised and verified by tests/DA1AF83D and tests/BD3DF323.

import json

# -- Geometry --
a = 0.5  # droplet (color-blob) radius
L = 4.0  # square domain side; droplet centered at the origin
Nx = 99
Ny = 99
w = 0.08  # diffuse color-interface half-width (~2 cells) for a well-resolved gradient

eps = 1.0e-9  # trace volume fraction of the second (identical) phase

# -- Equation of state (identical stiffened-gas fluids: gamma = 2) --
# The background pressure is kept well above the Laplace jump sigma/a so that the (uniform-pressure)
# initial condition is only mildly out of capillary equilibrium: sigma/a ~ 2 << p0 + p_inf = 125,
# giving an acoustic density perturbation rho'/rho ~ (sigma/a)/(rho c^2) < 1% (stable, low Mach).
gam = 2.0
p_inf = 100.0
cv = 12.5  # chosen so rho ~ 1 at the reference temperature (=> rho(x) = 10/(10 + x))
p0 = 25.0  # uniform background pressure (mechanical equilibrium)

# -- Imposed linear temperature field T(x) = T0 + gradT * x --
T0 = 10.0  # temperature at the droplet center (x = 0)
gradT = 1.0  # temperature gradient magnitude (in +x, toward the hot side)

# -- sigma(T) closure --
sigma0 = 1.0  # reference surface tension at T_ref
T_ref = T0  # reference temperature for the closure
dsigma_dT = -0.05  # dsigma/dT (negative: hotter -> lower sigma)

# -- Viscosity (mu* = 1) --
mu = 0.1  # dynamic viscosity of both phases; MFC takes Re = 1/mu

# Density that reproduces T(x) exactly from the stiffened-gas EOS at uniform pressure p0:
#   T = (p0 + p_inf) / ((gam - 1) * rho * cv)  =>  rho(x) = (p0 + p_inf) / ((gam - 1)*cv*T(x))
# With the chosen constants this is rho(x) = 10 / (10 + x). A single uniform-composition patch
# carries this analytic density consistently (no smoothing of partial densities), and the droplet
# is defined entirely by the smooth color function below.
rho_coeff = (p0 + p_inf) / ((gam - 1.0) * cv)  # = 10.0
rho_expr = f"(1.0 - {eps})*{rho_coeff}/({T0} + {gradT}*x)"  # ~ 10.0/(10.0 + x)
# Smooth color blob: c ~ 1 inside r < a, ~ 0 outside, with a tanh interface of half-width w.
cf_expr = f"0.5*(1.0 - tanh((sqrt(x**2 + y**2) - {a})/{w}))"

# Expected YGB terminal migration velocity (mu* = 1, k* = 1):
U_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * a / mu  # ~ 0.0333
c0 = (gam * (p0 + p_inf)) ** 0.5  # reference sound speed ~ 15.8
# Re = rho*U*a/mu ~ 0.17, Ca = mu*U/sigma0 ~ 0.0033, Ma = U/c0 ~ 0.0021 (Stokes / low-Ma / low-Ca)

# -- Time stepping (acoustic CFL ~ 0.13; run ~3 viscous times tau = rho*a^2/mu = 2.5 to reach
#    the quasi-steady migration plateau, after the brief acoustic capillary-relaxation transient) --
dx = L / (Nx + 1)
mydt = 3.0e-4

data = {
    # Logistics
    "run_time_info": "T",
    # Computational domain (centered on the droplet)
    "x_domain%beg": -L / 2,
    "x_domain%end": L / 2,
    "y_domain%beg": -L / 2,
    "y_domain%end": L / 2,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": 25000,
    "t_step_save": 500,
    # Simulation algorithm
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
    # Open (ghost-cell extrapolation) boundaries on all sides
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "num_patches": 1,
    "num_fluids": 2,
    # Physics: viscosity + temperature-dependent surface tension
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
    # Single uniform-composition patch spanning the domain: analytic linear-T density field and a
    # smooth color blob marking the droplet (surface tension acts on the color-function gradient).
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": L,
    "patch_icpp(1)%length_y": L,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": cf_expr,
}

print(json.dumps(data))
