#!/usr/bin/env python3
# DISPERSION variant of case.py: seed a SUPERPOSITION of surface modes n=1,2,3 and let them decay
# under pure surface diffusion. Because the operator is linear, each Fourier component decays
# independently at its own eigenvalue rate D_s*k_n^2 (k_n = n*2*pi/Lx), so a single run yields the
# whole rate(k)=D_s k^2 dispersion -- the canonical Laplace-Beltrami mode-decay benchmark. Measure with
#   ./mfc.sh run sweep.py -n 1 -t pre_process simulation
#   python3 measure.py sweep      # -> figures/dispersion.png
import json
import math

Lx, Ly = 2.0, 2.0
Nx = 48
dx = Lx / Nx
Ny = Nx
k1 = 2.0 * math.pi / Lx  # fundamental

gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1
rho0 = 0.2
sigma0 = 0.1
D_s = 0.2
eps = 0.2
# 1 + eps*(cos k1 x + cos k2 x + cos k3 x); each mode decays at D_s*k_n^2
surf_expr = f"1.0 + {eps:.6f}*(cos({k1:.9f}*x) + cos({2*k1:.9f}*x) + cos({3*k1:.9f}*x))"

c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
mydt = min(0.3 * dx / c_max, 0.3 * dx**2 / (4.0 * D_s))
t_step_stop = round(0.9 / mydt)  # ~1.8 tau of the slowest (n=1) mode
t_step_save = max(1, t_step_stop // 80)  # dense enough to resolve the fast (n=3) decay too

data = {
    "run_time_info": "T",
    "x_domain%beg": -Lx / 2,
    "x_domain%end": Lx / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
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
    "bc_x%beg": -1,
    "bc_x%end": -1,
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
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": Ly / 4,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly / 2,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha_rho(2)": 1.0e-8,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%surf_val": 0.0,
    "patch_icpp(2)%geometry": 3,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": -Ly / 4,
    "patch_icpp(2)%length_x": Lx,
    "patch_icpp(2)%length_y": Ly / 2,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8,
    "patch_icpp(2)%alpha_rho(2)": rho0,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%surf_val": surf_expr,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
}
print(json.dumps(data))
