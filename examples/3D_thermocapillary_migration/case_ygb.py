#!/usr/bin/env python3
# 3D thermocapillary droplet migration -- a VALIDATION of MFC's variable surface tension sigma(T)
# against the Young-Goldstein-Block analytic terminal velocity. Recovering u_YGB is a convergence
# claim, not a single number: v_t/u_YGB -> 1.0 as the three deficits vanish -- CONFINEMENT (finite
# box), FINITE Ma (temperature not perfectly invariant), and GRID. Re_M = rho*v_YGB*D/mu ~ 0.018 is
# already deep Stokes, so no Reynolds sweep is needed.
#
# WHY THE DECOUPLED THERMAL SCALAR (and not the density proxy of case.py)
# The sibling case.py fakes temperature through density (rho = rho_coeff/T(y)); that proxy is a
# TRANSPORTED field, so the drop's own flow advects the gradient it is meant to hold -- the local
# gradient collapses and reverses, and the rise velocity decays. This case carries temperature as an
# INDEPENDENT advected+diffused scalar T_s (thermal_scalar = T), decoupled from density. With
# thermal_scalar = T the surface-tension closure reads T_s directly (m_surface_tension.fpp), so
# sigma(T) is driven by the true temperature field, not an EOS artifact. Density is uniform, both
# fluids identical (mu* = k* = 1), so the ONLY thing driving the drop is the sigma(T) gradient -- the
# clean YGB setup.
#
# GEOMETRY MODES (env YGB_GEOM)
#   cube     (default) -- box W^3 cubic, drop CENTERED at y=0; maximal symmetric clearance so
#                         confinement is a clean one-parameter family in YGB_W -> extrapolate W->inf.
#   samareh  -- box 5D x 5D x 7.5D, drop offset 1.5D above the cold floor; reproduces Samareh Fig 6
#               (v_t/v_YGB ~ 0.95) as a confined-box anchor that de-risks the sweep vs the literature.
#
# PARAMETERS (Samareh, Mostaghimi & Moreau 2014, Sec. 4.1.1): D=1 sphere, rho_d=rho_b=0.2,
# mu_d=mu_b=0.1, sigma0=0.1, sigma_T=-0.1, |gradT|=2/15 -> v_YGB=8.889e-3. Slip walls on all six
# faces; isothermal Dirichlet gradient walls on y (cold floor / hot ceiling) hold the imposed field.
# T is shifted up by T0=10 (Samareh's T0=0) to keep it positive for the isothermal-wall validator;
# only gradT and sigma_T (which set v_YGB) are physical.
#
# GOTCHA (analytic-IC parser): MFC expands a bare `e` to Euler's number even inside `1e-9`. Keep
# `e`-notation out of every analytic patch string (T_temp_val); eps is folded into plain decimals.

import json
import os

# -- Variant selection (env vars) --
geom = os.environ.get("YGB_GEOM", "cube")  # cube (unbounded sweep) | samareh (confined anchor)
W_in_D = float(os.environ.get("YGB_W", "8"))  # cube box width in D (confinement knob; cube geom only)
Nx = int(os.environ.get("YGB_NX", "80"))  # cells per box WIDTH (grid knob)
Ma = float(os.environ.get("YGB_MA", "0.5"))  # thermal Marangoni number (conduction strength knob)
n_tr = float(os.environ.get("YGB_TR", "3.0"))  # run length in capillary-thermal times t_r
assert Ma > 0, "3D needs conduction (YGB_MA > 0): a frozen-T 3D rise has no plateau to validate."
assert geom in ("cube", "samareh"), f"YGB_GEOM must be 'cube' or 'samareh', got {geom!r}"

# -- Geometry --
D = 1.0  # droplet diameter
r = D / 2.0  # droplet radius = 0.5
if geom == "cube":
    W = W_in_D * D  # lateral extent (x, z)
    Ly = W  # cubic: rise (gradient) axis same length as lateral
    y_drop = 0.0  # drop centered -> maximal symmetric clearance from all walls
else:  # samareh
    W = 5.0 * D
    Ly = 7.5 * D
    y_drop = -Ly / 2.0 + 1.5 * D  # 1.5D above the cold floor = -2.25

dx = W / Nx  # isotropic cell size set by the box-width resolution
m = Nx - 1  # x (short)
n = round(Ly / dx) - 1  # y (rise axis); cube -> Nx-1, samareh -> round(1.5*Nx)-1
p = Nx - 1  # z (short) -- 3D

# -- Equation of state (two IDENTICAL stiffened-gas fluids, gamma=2; mu*=1, k*=1) --
gam = 2.0
p_inf, p0 = 32.0, 8.0  # with rho=0.2: c = sqrt(gam*(p0+p_inf)/rho) = 20
rho_b = 0.2  # uniform density EVERYWHERE (no density proxy); Samareh rho_d=rho_b=0.2
mu = 0.1  # dynamic viscosity of both phases; MFC takes Re = 1/mu
cv_b = 1.0  # EOS heat capacity (arbitrary; T_s is the thermal field, sigma reads T_s)

# -- Imposed linear temperature field T(y) = T0 + gradT*y, carried by the scalar T_s --
T0 = 10.0  # baseline (keeps T > 0 for the isothermal-wall validator)
gradT = 2.0 / 15.0  # |dT/dy| = 0.13333 (Samareh)
sigma0 = 0.1  # surface tension at T_ref
sigma_T = -0.1  # dsigma/dT
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop's INITIAL center
T_expr = f"{T0} + {gradT:.9f}*y"  # plain decimals only (no `e` notation)

eps = 1.0e-9  # trace volume fraction of the (identical) second phase

# -- Young-Goldstein-Block terminal speed and diagnostic time scales --
v_YGB = (2.0 / 15.0) * (-sigma_T) * gradT * r / mu  # = 8.889e-3
U_r = (-sigma_T) * gradT * r / mu  # Marangoni interfacial velocity scale = 0.06667
t_r = mu / abs(sigma_T * gradT)  # capillary-thermal time = 7.5

# -- Bulk thermal conduction (k*=1): alpha_T from the requested thermal Marangoni number Ma --
alpha_T = U_r * r / Ma  # thermal diffusivity = 0.03333/Ma
cp_b = gam * cv_b  # specific heat at constant pressure = 2.0
k_therm = alpha_T * rho_b * cp_b  # bulk conductivity; both fluids equal

# -- Time stepping: min(acoustic CFL, 3D explicit-diffusion limit d=3) --
c_max = (gam * (p0 + p_inf) / rho_b) ** 0.5  # = 20 (uniform density -> single sound speed)
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (6.0 * alpha_T))  # dt <= 0.35*dx^2/(2*d*alpha), d=3 (3D)
t_step_stop = int(round(n_tr * t_r / mydt))
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
    # Boundaries: slip walls (-2) on all six faces; isothermal gradient walls on y set below
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + sigma(T); T carried by an independent scalar T_s
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    "thermal_scalar": "T",
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "T_s_wrt": "T",
    "parallel_io": "T",
    # Continuous phase (fluid 1)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(1)%k_therm": k_therm,
    # Second phase (fluid 2) -- identical properties (mu* = 1, k* = 1)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_b,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%k_therm": k_therm,
    # Patch 1 -- background medium (3D cuboid spanning the domain): uniform density, color c=0,
    # linear T_s(y). NO density proxy -- alpha_rho is a plain uniform value.
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
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_b,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_b,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%T_temp_val": T_expr,
    # Patch 2 -- droplet (3D sphere): color c=1, same linear T_s(y) (T continuous across interface),
    # identical fluid to patch 1. Pressure carries the Laplace jump p0 + sigma/r (no t=0 transient).
    "patch_icpp(2)%geometry": 8,  # 3D sphere
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
    "patch_icpp(2)%pres": p0 + sigma0 / r,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_b,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_expr,
    # Isothermal Dirichlet gradient walls pin T_s to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T0 + gradT * (-Ly / 2.0),  # cold floor (y%beg)
    "bc_y%Twall_out": T0 + gradT * (Ly / 2.0),  # hot ceiling (y%end)
}

print(json.dumps(data))
