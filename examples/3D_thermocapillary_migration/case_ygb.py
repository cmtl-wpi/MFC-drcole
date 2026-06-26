#!/usr/bin/env python3
# 3D thermocapillary droplet migration -- a VALIDATION of MFC's variable surface tension sigma(T)
# against the Young-Goldstein-Block analytic terminal velocity. Recovering u_YGB is a convergence
# claim, not a single number: v_t/u_YGB -> 1.0 as the three deficits vanish -- CONFINEMENT (finite
# box), FINITE Ma (temperature not perfectly invariant), and GRID. Re_M = rho*v_YGB*D/mu ~ 0.018 is
# already deep Stokes, so no Reynolds sweep is needed.
#
# TEMPERATURE = DENSITY PROXY + CONDUCTION
# MFC is compressible: at the low migration Mach here, pressure is ~uniform and the EOS locks
# T = (p+p_inf)/((gam-1)*rho*cv), so the imposed gradient T(y) is encoded as rho(y) = rho_coeff/T(y).
# A bare density proxy is a TRANSPORTED field -- the drop's own flow advects the gradient it should
# hold -- so we ALSO run bulk conduction + isothermal gradient walls, which actively restore the
# field (the finite-Ma realization; YGB_MA sets the conduction strength). Both fluids are identical,
# so one analytic patch (shared smooth sphere eta(x,y,z) driving color, the Laplace pressure jump,
# and density together) makes the EOS-recovered T exactly linear everywhere -- drop included. Honest
# caveat: density stratifies as ~1/T (a compressibility artifact absent in the incompressible
# reference; magnitude ~ dT/T, shrunk by the T0 offset).
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
# T is shifted up by T0=10 (Samareh's T0=0) to keep rho positive; only gradT and sigma_T (which set
# v_YGB) are physical.
#
# GOTCHA (analytic-IC parser): MFC expands a bare `e` to Euler's number even inside `1e-9`, and a
# bare `r`/`eps` to the patch radius/epsilon. Keep them out of every analytic patch string (use
# plain decimals); tanh/sqrt are safe.

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
mu = 0.1  # dynamic viscosity of both phases; MFC takes Re = 1/mu

# -- Imposed linear temperature field T(y) = T0 + gradT*y, encoded as rho(y) = rho_coeff/T(y) --
T0 = 10.0  # baseline (keeps rho > 0 and the isothermal-wall validator happy)
gradT = 2.0 / 15.0  # |dT/dy| = 0.13333 (Samareh)
sigma0 = 0.1  # surface tension at T_ref
sigma_T = -0.1  # dsigma/dT
rho_b = 0.2  # density at the reference temperature T0 (Samareh rho_d=rho_b=0.2)
rho_coeff = rho_b * T0  # = 2.0
cv_b = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # closes the EOS so rho(T0) = rho_b; = 20
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop's INITIAL center

eps = 1.0e-9  # trace volume fraction of the (identical) second phase

# -- Young-Goldstein-Block terminal speed and diagnostic time scales --
v_YGB = (2.0 / 15.0) * (-sigma_T) * gradT * r / mu  # = 8.889e-3
U_r = (-sigma_T) * gradT * r / mu  # Marangoni interfacial velocity scale = 0.06667
t_r = mu / abs(sigma_T * gradT)  # capillary-thermal time = 7.5

# -- Bulk thermal conduction (k*=1): alpha_T from the requested thermal Marangoni number Ma --
alpha_T = U_r * r / Ma  # thermal diffusivity = 0.03333/Ma
cp_b = gam * cv_b  # specific heat at constant pressure
k_therm = alpha_T * rho_b * cp_b  # bulk conductivity at the reference state; both fluids equal

# -- One analytic patch: smooth sphere eta(x,y,z) ~ 1 in the drop / 0 outside drives color, the
# Laplace pressure jump, and density together, so the EOS-recovered T is the imposed linear field
# everywhere (rho tracks p, the (p+p_inf) factors cancel). Hardcode center/radius as decimals. --
xc_d, yc_d, zc_d, r_d = 0.0, y_drop, 0.0, r
w_if = 0.75 * dx  # interface half-width (~3-cell transition)
laplace = sigma0 / r  # Laplace pressure jump sigma/r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2 + (z - ({zc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
cf_expr = f"({eta})"
pres_expr = f"{p0:.9f} + {laplace:.9f}*({eta})"
rho_num = (1.0 - eps) * rho_coeff
rho_expr = f"{rho_num:.9f}*({p0 + p_inf:.9f} + {laplace:.9f}*({eta}))/({p0 + p_inf:.9f}*({T0} + {gradT:.9f}*y))"

# -- Time stepping: min(acoustic CFL, 3D explicit-diffusion limit d=3) --
rho_min = rho_coeff / (T0 + gradT * Ly / 2.0)  # hot wall: lowest density, max sound speed
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
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
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(1)%k_therm": k_therm,
    # Second phase (fluid 2) -- identical properties (mu* = 1, k* = 1)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_b,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%k_therm": k_therm,
    # Single analytic patch (3D cuboid spanning the domain). The drop lives entirely in eta(x,y,z):
    # cf, pres, and alpha_rho all share it, so rho tracks the (jump-carrying) pressure and the
    # EOS-recovered T is exactly the imposed linear field everywhere -- drop included.
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
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": cf_expr,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T0 + gradT * (-Ly / 2.0),  # cold floor (y%beg)
    "bc_y%Twall_out": T0 + gradT * (Ly / 2.0),  # hot ceiling (y%end)
}

print(json.dumps(data))
