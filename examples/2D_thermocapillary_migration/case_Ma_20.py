#!/usr/bin/env python3
# Thermocapillary migration of a 2D drop at low Marangoni number -- Samareh, Mostaghimi & Moreau,
# Int. J. Heat Mass Transfer 73 (2014) 616-626, Sec. 4.1.2 / Fig. 7 (test case of Nas & Tryggvason).
# A drop of diameter D = 1 starts 1D above the cold floor of a 2D x 4D slip-wall box. With Ma finite
# the energy equation is coupled: the drop's motion distorts the temperature field, a thermal
# boundary layer forms at the interface, and sigma responds to the evolving local T. Drop and bulk
# are genuinely different fluids (all material properties at ratio 0.5), so density cannot also
# encode T -- temperature is carried by an independent advected/diffused scalar T_s, and sigma(T)
# reads T_s. Non-dimensional targets: Re = 5, Ma = 20, Ca = 0.01666; the migration velocity
# U* = U/U_r peaks near 0.13 at t* = t/t_r ~ 5.
#
# rho_b and mu_b are arbitrary (only Re, Ma, Ca are physical); they are chosen for a deeply
# incompressible state. Keep `e`-notation out of analytic strings (MFC reads `e` as Euler's number).

import json

# Non-dimensional targets (Samareh Fig. 7)
Re = 5.0
Ma = 20.0
Ca = 0.01666
prop_ratio = 0.5  # droplet / bulk material-property ratio

# Geometry: D = 1 drop in a 2D wide x 4D tall box; drop 1D above the cold floor, gradient axis = y
D = 1.0
r = D / 2.0
Wx = 2.0 * D
Hy = 4.0 * D
y_bottom = -Hy / 2.0
y_drop = y_bottom + 1.0 * D

Nx = 64  # cells across the box width (Samareh used 64, 128)
dx = Wx / Nx
Ny = round(Hy / dx)

# Bulk reference scales (arbitrary; only Re, Ma, Ca are physical)
rho_b = 1.0
mu_b = 0.02
gam = 2.0
cv_b = 1.0

# Invert the non-dimensional numbers for the physical surface-tension and conduction properties
G = Re * mu_b**2 / (rho_b * r**2)  # |sigma_T*gradT| = 0.008
gradT = 1.0 / Hy  # T runs 0 (cold) .. 1 (hot) across the box height = 0.25
sigma_T = -G / gradT  # dsigma/dT < 0 = -0.032
sigma0 = G * r / Ca  # surface tension at T_ref ~ 0.240
alpha_b = G * r**2 / (mu_b * Ma)  # bulk thermal diffusivity = 0.005
cp_b = gam * cv_b
k_b = alpha_b * rho_b * cp_b  # bulk conductivity = 0.01

# Droplet properties = prop_ratio * bulk
rho_d = prop_ratio * rho_b
mu_d = prop_ratio * mu_b
cv_d = prop_ratio * cv_b
k_d = prop_ratio * k_b

U_r = G * r / mu_b  # Marangoni velocity scale = 0.2
t_r = mu_b / G  # capillary-thermal time = 2.5

# Stiffened-gas EOS; background pressure well above the Laplace jump sigma0/r
p0, p_inf = 5.0, 20.0
c_max = (gam * (p0 + p_inf) / rho_d) ** 0.5  # droplet (lowest density) sets the max sound speed ~ 10

# Imposed field T(y) = T_base + gradT*(y - y_bottom); T_base > 0 so the isothermal-wall validator
# passes (inert for the shift-invariant scalar T_s). T runs 1 (cold) .. 2 (hot).
T_base = 1.0


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the drop's initial position
T_expr = f"{T_base} + {gradT}*(y - ({y_bottom}))"  # plain decimals only (no `e` notation)

eps = 1.0e-9  # trace volume fraction of the "other" fluid in each patch

# Time stepping: min(acoustic CFL, explicit-conduction limit)
mydt = 0.35 * dx / c_max
alpha_max = max(alpha_b, k_d / (rho_d * gam * cv_d))  # fastest-diffusing phase sets the diffusion dt
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_max))  # 2D explicit diffusion number
t_step_stop = round(15.0 * t_r / mydt)  # 15 capillary-thermal times
t_step_save = max(1, t_step_stop // 100)

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
    # Closed slip-wall box (Samareh / Nas & Tryggvason); isothermal gradient walls set below
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
    # Fluid 1 (bulk)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (droplet), all material properties at prop_ratio * bulk
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    "fluid_pp(2)%k_therm": k_d,
    # Patch 1: bulk medium, density rho_b, color c = 0, linear T_s(y)
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
    # Patch 2: droplet (circle), fluid 2 at rho_d, color c = 1, same T_s(y) as the bulk (T is
    # continuous across the interface). Pressure carries the Laplace jump p0 + sigma/r.
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 + sigma0 / r,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_d,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_expr,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling)
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Hy / 2),
    "bc_y%Twall_out": T_of_y(Hy / 2),
}

print(json.dumps(data))
