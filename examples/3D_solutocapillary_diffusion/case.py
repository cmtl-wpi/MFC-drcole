#!/usr/bin/env python3
# 3D VALIDATION: surfactant surface diffusion on a SPHERE vs the exact Laplace-Beltrami rate.
#
# A passive insoluble surfactant (sigma_model = 0 -> constant sigma, no Marangoni; the drop stays
# static) is seeded on a sphere of radius R with the l=1 spherical-harmonic mode Y_1 ~ z/R:
#   Gamma(z) = Gamma0 * (1 + eps * z/R).
# Under pure surface diffusion an l-mode decays exactly as exp(-l(l+1) D_s/R^2 t), so for l=1 the
# amplitude M1(t) = sum(Gamma_tilde * z) decays at rate 2 D_s/R^2. This is the canonical 3D
# surface-diffusion benchmark (the genuinely 2D surface Laplacian, curved). The measured rate
# approaches 2 D_s/R^2 as the interface is resolved -- run at several NX (see run_convergence.sh).
#
# Resolution is set by the MFC_NX environment variable so one build serves the whole sweep (the
# analytic IC z/R does not depend on NX). Keep `e`-notation out of analytic strings.
import json
import math
import os

R = 0.5
W = 3.0  # box = 6 R; sphere well separated from the slip walls
Nx = int(os.environ.get("MFC_NX", "64"))
dx = W / Nx
Ny = Nz = Nx

gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1
rho0 = 0.2
sigma0 = 0.1
laplace = 2.0 * sigma0 / R  # sphere Laplace jump 2 sigma/R (suppresses the startup ring)
D_s = 0.2  # surface diffusivity; l=1 decay time tau = R^2/(2 D_s) = 0.625
eps = 0.3
surf_expr = f"1.0 + {eps / R:.6f}*z"  # Gamma0=1, l=1 mode ~ 1 + eps*z/R

c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
mydt = min(0.3 * dx / c_max, 0.3 * dx**2 / (6.0 * D_s))  # min(acoustic, 3D surface-diffusion) CFL
t_step_stop = round(0.9 / mydt)  # ~1.5 tau
t_step_save = max(1, t_step_stop // 60)

data = {
    "run_time_info": "T",
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -W / 2,
    "y_domain%end": W / 2,
    "z_domain%beg": -W / 2,
    "z_domain%end": W / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": Nz - 1,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
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
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "surfactant": "T",
    "surf_diff": D_s,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": 1.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": 1.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    # Patch 1: bath (cuboid spanning the domain), no surfactant, c=0
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": W,
    "patch_icpp(1)%length_z": W,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha_rho(2)": 1.0e-8,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%surf_val": 0.0,
    # Patch 2: sphere, c=1, l=1 surfactant mode ~ z/R
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p0 + laplace,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8,
    "patch_icpp(2)%alpha_rho(2)": rho0,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%surf_val": surf_expr,
}
print(json.dumps(data))
