#!/usr/bin/env python3
# TC3 -- Samareh, Mostaghimi & Moreau (2014), Sec. 4.2: thermocapillary migration of a droplet at
# LARGE Marangoni number, matched to the Life & Microgravity Science (LMS) Space Shuttle experiment
# (Hadland et al.). A Fluorinert FC-75 drop (D = 10.7 mm) rises through Dow-Corning DC-200 silicon oil
# in a 60 x 45 x 45 mm cell with a 1000 K/m vertical temperature gradient (cold floor 283 K, hot ceiling
# 343 K). Re = 17.79, Ma = 1723 at the reference T0 = 313 K. Samareh's Figs 8/13 compare the rise
# velocity vs distance from the cold wall.
#
# THIS IS THE CASE THAT NEEDED mu(T). TC3 is the one Samareh scenario whose blocker was the physics, not
# the dimensionality: the silicon oil's viscosity varies substantially across the 60 K cell (the paper:
# density and cp variations are negligible, only mu matters), and that temperature-dependent drag is what
# produces the experiment's non-monotonic rise-velocity "loop". It needs THREE features, now all present:
# bulk conduction (thermal_conduction), the sigma(T) closure (sigma_model=1), and temperature-dependent
# viscosity mu(T) = exp(C + D/T) (visc_model=1, the new feature).
#
# TEMPERATURE HANDLING. The drop and bulk have DIFFERENT, ~constant densities (rho_d=1727.7, rho_b=918.3),
# so the density-proxy T encoding used by the zero-Ma cases (rho = rho_coeff/T) does not work. Instead T
# is carried by the INDEPENDENT advected+diffused scalar T_s (thermal_scalar = T), decoupled from density.
# Both sigma(T) and mu(T) read T_s in this mode (T in KELVIN, so the paper's Arrhenius C/D apply directly).
# The stiffened-gas EOS then only sets the (low-Mach, insensitive) acoustics, so pi_inf/p0 are softened
# for a tractable dt -- the migration is Mach ~ 1e-4 and independent of the sound speed.
#
# PROPERTIES (Samareh Sec. 4.2, Eq. 30). mu = exp(C + D/T):
#   silicon oil (bulk):  C = -10.17, D = 1643;  rho = 918.3,  k = 0.13389, cp = 1778.2
#   Fluorinert  (drop):  C = -11.76, D = 1540;  rho = 1727.7, k = 0.063,   cp = 1047.0
#   sigma0 = 0.007 N/m, sigma_T = -3.6e-5 N/m/K, T0 = 313 K.   (all SI: m, kg/m^3, Pa.s, W/mK, J/kgK)
#
# SCOPE. A converged Fig-8/13 comparison is a heavy 3D run (Samareh used up to 240 x 640 x 240). This
# case validates and smoke-runs immediately at a coarse grid; the headline comparison is a production run.

import json
import math
import os

# -- Variant selection --
nx_cross = int(os.environ.get("TC3_NX", "30"))  # cells across the 45 mm cross-section (Samareh used much finer)
uniform_ic = os.environ.get("TC3_UNIFORM", "0") == "1"  # sub-case A (0): drop equilibrated/linear; sub-case B (1): drop uniform 298 K, bulk stays linear
n_tr = float(os.environ.get("TC3_TR", "1"))  # run length in capillary-thermal times t_r
two_d = os.environ.get("TC3_2D", "0") == "1"  # 1 = 2D cylinder analogue (tractable Fig-8 physics, like TC1 Fig 5)
const_visc = os.environ.get("TC3_VISC", "1") == "0"  # 1 = mu(T) Arrhenius (default); 0 = constant-mu control run

# -- Geometry (LMS test cell; gradient along +y) --
D = 10.7e-3  # droplet diameter (m)
r = D / 2.0  # droplet radius = 5.35 mm
Wx = 45.0e-3  # cross-section (x, z)
Ly = 60.0e-3  # gradient / rise axis (y), cold floor -> hot ceiling
y_drop = -Ly / 2 + 15.0e-3  # released ~15 mm above the cold wall (experiment's start, Fig 8 x-axis)

dx = Wx / nx_cross  # cross-section (x, z) cell size
# Rise-axis (y) resolution is independent of the cross-section. Samareh's meshes are ~2x finer along the
# gradient (180x480x180, 240x640x240 -> y-count = 2.67x the cross-section, not the isotropic 1.33x). TC3_NY
# defaults to the isotropic count (dy = dx); set TC3_NX=180 TC3_NY=480 (or 240 / 640) to match the paper.
ny_rise = int(os.environ.get("TC3_NY", str(round(Ly / dx))))  # cells along the 60 mm rise axis
dy = Ly / ny_rise
m = nx_cross - 1  # x
p = 0 if two_d else nx_cross - 1  # z (0 in the 2D cylinder analogue)
n = ny_rise - 1  # y (gradient axis)

cf_smooth_coeff = 0.5  # ~2-cell diffuse color interface

# -- Imposed temperature field (Kelvin) --
T_c, T_h = 283.0, 343.0  # cold floor / hot ceiling wall temperatures
gradT = (T_h - T_c) / Ly  # = 1000 K/m
T0 = 313.0  # reference temperature (cell center), sigma_T_ref
sigma0 = 0.007  # surface tension at T0 (N/m)
sigma_T = -3.6e-5  # dsigma/dT (N/m/K)

# -- Per-fluid material properties (fluid 1 = silicon oil bulk, fluid 2 = Fluorinert drop) --
rho_b, rho_d = 918.3, 1727.7  # densities (kg/m^3)
k_b, k_d = 0.13389, 0.063  # thermal conductivities (W/mK)
cp_b, cp_d = 1778.2, 1047.0  # specific heats (J/kgK)
C_b, Db = -10.17, 1643.0  # silicon oil Arrhenius mu(T) = exp(C + D/T)
C_d, Dd = -11.76, 1540.0  # Fluorinert  Arrhenius mu(T)

# -- Stiffened-gas EOS: softened for a stable low-Mach state (acoustics decoupled from T via thermal_scalar) --
gam = 2.0  # ratio of specific heats (EOS knob; cp = gam*cv -> cv = cp/gam)
cv_b, cv_d = cp_b / gam, cp_d / gam  # specific heat at constant volume per fluid
c_snd = 30.0  # softened sound speed (m/s); migration Mach ~ U_r/c ~ 1e-3, so c only sets the acoustic CFL
p0 = 1.0e5  # background pressure (~1 atm); Laplace jump sigma/r ~ 1.3 Pa is negligible vs p0
p_inf_b = rho_b * c_snd**2 / gam - p0  # stiffening pressure so c = sqrt(gam*(p0+p_inf)/rho) = c_snd
p_inf_d = rho_d * c_snd**2 / gam - p0

# -- Migration scales (Samareh): U_r = |sigma_T*gradT|*r/mu_b(T0), t_r = mu_b(T0)/|sigma_T*gradT| --
mu_b0 = math.exp(C_b + Db / T0)  # silicon oil viscosity at T0 = 313 K (~7.3e-3 Pa.s)
G = abs(sigma_T * gradT)  # Marangoni stress scale
U_r = G * r / mu_b0  # ~ 26.5 mm/s
t_r = mu_b0 / G  # ~ 0.20 s

# -- Viscosity model: mu(T) Arrhenius (default) or a constant-mu CONTROL frozen at the IC reference T --
# The control freezes mu at the temperature the DROP sits in at t=0 (so the two runs are identical at t=0
# and diverge ONLY because mu(T) lets the local viscosity change as the drop migrates into warmer oil).
visc_model = 0 if const_visc else 1  # 0 = constant mu (control), 1 = mu(T) = exp(C + D/T)
T_visc_ref = 298.0 if uniform_ic else (T0 + gradT * y_drop)  # T at the drop center at t=0
mu_b_ref = math.exp(C_b + Db / T_visc_ref)  # bulk constant-mu (= mu at the drop's start)
mu_d_ref = math.exp(C_d + Dd / T_visc_ref)  # drop constant-mu

# -- Bulk conduction diffusivity (slowest fluid sets the explicit-diffusion dt cap) --
alpha_b = k_b / (rho_b * cp_b)  # ~ 8.2e-8 m^2/s (large Ma -> tiny diffusivity, dt cap not binding)
alpha_d = k_d / (rho_d * cp_d)

# -- Time stepping: acoustic CFL + 3D explicit-diffusion cap, evaluated on the SMALLEST cell so an
#    anisotropic (dy < dx) mesh stays stable (d = 3) --
dmin = min(dx, dy)
mydt = 0.35 * dmin / c_snd
mydt = min(mydt, 0.35 * dmin**2 / (6.0 * max(alpha_b, alpha_d)))
t_end = n_tr * t_r
t_step_stop = int(round(t_end / mydt))
n_save = int(os.environ.get("TC3_NSAVE", "60"))  # # of snapshots (denser sampling averages over acoustic ringing)
t_step_save = max(1, t_step_stop // n_save)

# -- Optional chunked restart: run [TC3_TSTART, TC3_TSTOP] instead of [0, t_step_stop], so a long run can be
#    completed in short pieces (each piece restarts from the previous chunk's checkpoint at TC3_TSTART, which
#    must be a saved step). Defaults span the whole run, so the normal single-shot run is unchanged.
t_step_start = int(os.environ.get("TC3_TSTART", "0"))
t_step_stop = int(os.environ.get("TC3_TSTOP", str(t_step_stop)))

# -- Temperature scalar IC. The bulk liquid ALWAYS starts on the imposed linear profile; the sub-cases
#    differ only in the DROP's initial internal temperature (the experiment's drop temperature was unknown):
#      Sub-case A (TC3_UNIFORM=0, Sec. 4.2.1): drop equilibrated to the local linear field (matches bulk).
#      Sub-case B (TC3_UNIFORM=1, Sec. 4.2.2): drop injected at a uniform 298 K; the bulk stays linear.
GRAD = "y"
T_linear = f"{T0} + {gradT:.6f}*{GRAD}"  # imposed linear field (T0 at y=0)
T_bulk = T_linear  # the bulk is always the linear gradient
T_drop = "298.0" if uniform_ic else T_linear  # uniform 298 K only in sub-case B

eps = 1.0e-9  # trace volume fraction of the other phase

data = {
    "run_time_info": "T",
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
    "z_domain%beg": -Wx / 2,
    "z_domain%end": Wx / 2,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": t_step_start,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
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
    # Slip (-2) cross-section walls; the gradient (y) walls are isothermal (set below)
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity (mu(T) per fluid) + sigma(T) + bulk conduction + independent temperature scalar
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T0,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    "thermal_scalar": "T",
    "T_s_wrt": "T",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 -- silicon oil (bulk): Arrhenius mu(T), real k/cp
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf_b / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b_ref,
    "fluid_pp(1)%k_therm": k_b,
    "fluid_pp(1)%visc_model": visc_model,
    "fluid_pp(1)%visc_c": C_b,
    "fluid_pp(1)%visc_d": Db,
    # Fluid 2 -- Fluorinert (drop): Arrhenius mu(T), real k/cp
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d_ref,
    "fluid_pp(2)%k_therm": k_d,
    "fluid_pp(2)%visc_model": visc_model,
    "fluid_pp(2)%visc_c": C_d,
    "fluid_pp(2)%visc_d": Dd,
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
    "patch_icpp(1)%T_temp_val": T_bulk,
    # Patch 2 -- Fluorinert drop (3D sphere), distinct density/properties from the bulk
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": cf_smooth_coeff,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_d,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_drop,
    # Isothermal Dirichlet walls on the gradient axis pin the cold floor / hot ceiling
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_c,
    "bc_y%Twall_out": T_h,
}

if two_d:
    # 2D cylinder analogue (Samareh's Fig-5-style plane of the 3D drop): a circle in a rectangular cell,
    # no z-direction. Same mu(T)/sigma(T)/conduction physics; tractable enough to run to migration.
    data["patch_icpp(1)%geometry"] = 3  # rectangle (was 3D cuboid, geometry 9)
    data["patch_icpp(2)%geometry"] = 2  # circle    (was 3D sphere, geometry 8)
    for key in (
        "z_domain%beg",
        "z_domain%end",
        "bc_z%beg",
        "bc_z%end",
        "patch_icpp(1)%z_centroid",
        "patch_icpp(1)%length_z",
        "patch_icpp(1)%vel(3)",
        "patch_icpp(2)%z_centroid",
        "patch_icpp(2)%vel(3)",
    ):
        data.pop(key, None)

print(json.dumps(data))
