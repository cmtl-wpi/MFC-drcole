#!/usr/bin/env python3
# Thermocapillary droplet migration, 2D finite-Marangoni -- Samareh, Mostaghimi & Moreau,
# Int. J. Heat Mass Transfer 73 (2014) 616-626 (Sec. 4.1.2 / Fig 7; orig. Nas & Tryggvason 2003).
#
# Unlike the frozen-T Fig 5 case (case.py), here Ma is finite, so the energy equation is coupled:
# the drop's motion distorts T, a thermal boundary layer forms at the interface, and sigma responds
# to the evolving local T. This is the case that exercises bulk Fourier conduction. Drop and bulk
# are genuinely different fluids (all material properties at ratio 0.5), so density cannot also
# encode T -- T is carried by the independent advected+diffused scalar T_s (thermal_scalar), and
# sigma(T) reads T_s. Targets: Re=5, Ma=20, Ca=0.01666; U* = U/U_r peaks ~0.13 near t* ~ 5.
#
# Non-dimensional inversion (length r0, velocity U_r = G*r0/mu_b, time t_r = mu_b/G, G = |sigma_T*gradT|):
#   G = Re*mu_b^2/(rho_b*r0^2), sigma0 = G*r0/Ca, alpha_b = G*r0^2/(mu_b*Ma); gradT fixes T over 0..1.
# rho_b, mu_b are arbitrary (only Re/Ma/Ca are physical), chosen for a deeply incompressible state.
#
# Env vars:
#   FIG7_NX         cells across the 2D box width (default 64; Samareh used 64, 128)
#   FIG7_TR         run length in capillary-thermal times t_r (default 15)
#   FIG7_RATIO      droplet/bulk material-property ratio (default 0.5)
#   FIG7_COND       1 = bulk conduction on (default); 0 = advect T_s only (diagnostic)
#   FIG7_ADIABATIC  1 = adiabatic walls; default = isothermal Dirichlet (faithful Samareh setup)
#   FIG7_UNBALANCED 1 = bare uniform-pressure IC (rings); default balances the t=0 Laplace jump
#
# Note: keep `e`-notation out of analytic patch strings -- MFC's IC parser reads `e` as Euler's number.

import json
import os

# Variant selection (defaults = Samareh Fig 7 at the coarse grid)
width_cells = int(os.environ.get("FIG7_NX", "64"))
n_tr = float(os.environ.get("FIG7_TR", "15"))
prop_ratio = float(os.environ.get("FIG7_RATIO", "0.5"))  # droplet/bulk material-property ratio
conduction = os.environ.get("FIG7_COND", "1") == "1"  # off -> advection of T_s only
adiabatic = os.environ.get("FIG7_ADIABATIC", "0") == "1"  # adiabatic vs isothermal Dirichlet walls
unbalanced_ic = os.environ.get("FIG7_UNBALANCED", "0") == "1"

# Samareh's Fig 7 non-dimensional targets
Re = 5.0
Ma = 20.0
Ca = 0.01666

# Geometry: 2D drop in a 2D-wide x 4D-tall box, rise/gradient axis = y
D = 1.0
r0 = D / 2.0  # length scale
Wx = 2.0 * D
Hy = 4.0 * D
y_bottom = -Hy / 2.0  # cold floor
y_drop = y_bottom + 1.0 * D  # drop center 1D above the cold floor (Samareh)

dx = Wx / width_cells
height_cells = round(Hy / dx)
m = width_cells - 1
n = height_cells - 1
p = 0

cf_smooth_coeff = 0.5  # ~2-cell diffuse color interface

# Bulk-phase reference scales (arbitrary; only Re, Ma, Ca are physical)
rho_b = 1.0
mu_b = 0.02  # small -> deeply incompressible Marangoni flow
gam = 2.0  # stiffened-gas exponent (same for both fluids)
cv_b = 1.0  # bulk cv (enters conduction via cp = gam*cv)

# Invert the non-dimensional numbers for the physical (sigma, gradT, alpha, k)
G = Re * mu_b**2 / (rho_b * r0**2)  # |sigma_T*gradT| = 0.008
gradT = 1.0 / Hy  # T runs 0 (cold) .. 1 (hot) over the 4D box = 0.25
sigma_T = -G / gradT  # dsigma/dT < 0 = -0.032
sigma0 = G * r0 / Ca  # surface tension at T_ref ~ 0.240
alpha_b = G * r0**2 / (mu_b * Ma)  # bulk thermal diffusivity = 0.005
cp_b = gam * cv_b  # = 2.0
k_b = alpha_b * rho_b * cp_b  # bulk conductivity = 0.01

# Droplet properties = prop_ratio * bulk
rho_d = prop_ratio * rho_b
mu_d = prop_ratio * mu_b
cv_d = prop_ratio * cv_b
k_d = prop_ratio * k_b

# Diagnostic scales (Samareh's normalization)
U_r = G * r0 / mu_b  # Marangoni velocity scale = 0.2
t_r = mu_b / G  # capillary-thermal time = 2.5

# Equation of state: stiffened gas; background pressure well above the Laplace jump sigma0/r0
p0, p_inf = 5.0, 20.0
rho_min = rho_d  # lowest density (droplet) sets the max sound speed / acoustic CFL
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5  # ~ 10

# Imposed field T(y) = T_base + gradT*(y - y_bottom); T_base > 0 so the isothermal-BC validator
# passes (inert for the decoupled scalar T_s, which is shift-invariant). T runs 1 (cold) .. 2 (hot).
T_base = 1.0


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the drop's initial position
T_expr = f"{T_base} + {gradT}*(y - ({y_bottom}))"  # plain decimals only (no `e` notation)
eps = 1.0e-9  # trace volume fraction of the "other" fluid in each patch

# Time stepping: min(acoustic CFL, explicit-conduction limit)
mydt = 0.35 * dx / c_max
alpha_max = max(alpha_b, k_d / (rho_d * gam * cv_d))  # fastest-diffusing phase sets the conduction dt
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_max))  # explicit diffusion number, d = 2 (2D)
t_end = n_tr * t_r
t_step_stop = int(round(t_end / mydt))
t_step_save = max(1, t_step_stop // 100)  # ~100 snapshots

data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Hy / 2,
    "y_domain%end": Hy / 2,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation Algorithm Parameters
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
    # Boundaries: closed slip-wall box (Samareh / Nas & Tryggvason)
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + sigma(T); T carried by an independent scalar T_s
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_scalar": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "T_s_wrt": "T",
    "parallel_io": "T",
    # Continuous phase (fluid 1 = bulk)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    # Dispersed phase (fluid 2 = droplet), all material properties at prop_ratio * bulk
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    # Patch 1: bulk medium, density rho_b, color c=0, linear T_s(y)
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Hy,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_b,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_d,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%T_temp_val": T_expr,
    # Patch 2: droplet (circle), fluid 2 at rho_d, color c=1, same T_s(y) as the bulk (T is
    # continuous across the interface). Pressure is the Laplace overpressure p0 + sigma/r unless
    # unbalanced, removing the t=0 acoustic ring.
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%radius": r0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": cf_smooth_coeff,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 if unbalanced_ic else p0 + sigma0 / r0,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_d,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_expr,
}

if conduction:
    # Finite Ma: solve energy transport for T_s at alpha = k/(rho cp)
    data.update(
        {
            "thermal_conduction": "T",
            "fluid_pp(1)%k_therm": k_b,
            "fluid_pp(2)%k_therm": k_d,
        }
    )
    if not adiabatic:
        # Faithful setup: isothermal Dirichlet on top/bottom walls pins T to the imposed gradient,
        # so the walls sink the drop's thermal wake and migration reaches a clean plateau
        data.update(
            {
                "bc_y%isothermal_in": "T",
                "bc_y%isothermal_out": "T",
                "bc_y%Twall_in": T_of_y(-Hy / 2),
                "bc_y%Twall_out": T_of_y(Hy / 2),
            }
        )

print(json.dumps(data))
