#!/usr/bin/env python3
# 2D VALIDATION: surfactant surface diffusion on a CIRCLE vs the exact Laplace-Beltrami rate.
#
# The curved-but-cheap analog of the 3D sphere study: a passive insoluble surfactant (sigma_model = 0
# -> constant sigma, no Marangoni; the drop stays static) is seeded on a circle of radius R with the
# m=1 azimuthal mode Gamma = Gamma0*(1 + eps*x/R) (x/R = cos(theta)). On a circle the m-th mode decays
# exactly as exp(-m^2 D_s/R^2 t), so for m=1 the amplitude M1(t) = sum(Gamma_tilde*x) decays at rate
# D_s/R^2 (note: the CIRCLE value m^2 D_s/R^2, not the sphere's l(l+1) D_s/R^2). Because it is 2D, the
# interface can be resolved far more finely than the sphere. Running a resolution sweep
# (run_convergence.sh) exposes that the whole-field-moment rate estimator does NOT converge cleanly on a
# curved interface -- it brackets D_s/R^2 rather than pinning it (see README); a proper interfacial
# measurement is future work.
#
# Resolution is set by MFC_NX so one build serves the whole sweep. Keep `e`-notation out of analytic strings.
import json
import os

R = 0.5
W = 3.0  # box = 6 R
Nx = int(os.environ.get("MFC_NX", "96"))
dx = W / Nx
Ny = Nx

gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1
rho0 = 0.2
sigma0 = 0.1
laplace = sigma0 / R  # 2D circle Laplace jump sigma/R
D_s = 0.2  # surface diffusivity; m=1 decay time tau = R^2/D_s = 1.25
eps = 0.3
surf_expr = f"1.0 + {eps / R:.6f}*x"  # Gamma0=1, m=1 mode ~ 1 + eps*x/R

c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
mydt = min(0.3 * dx / c_max, 0.3 * dx**2 / (4.0 * D_s))  # min(acoustic, 2D surface-diffusion) CFL
t_step_stop = round(1.5 / mydt)  # ~1.2 tau
t_step_save = max(1, t_step_stop // 60)

data = {
    "run_time_info": "T",
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -W / 2,
    "y_domain%end": W / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
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
    # Patch 1: bath (full box), no surfactant, c=0
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": W,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha_rho(2)": 1.0e-8,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%surf_val": 0.0,
    # Patch 2: circle, c=1, m=1 surfactant mode ~ x/R
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 + laplace,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8,
    "patch_icpp(2)%alpha_rho(2)": rho0,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%surf_val": surf_expr,
}
print(json.dumps(data))
