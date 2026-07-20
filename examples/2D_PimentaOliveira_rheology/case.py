#!/usr/bin/env python3
# Pimenta & Oliveira (2021) emulsion-rheology benchmark (guide Table 3.3): a single surfactant drop in
# simple shear, from which the interfacial contribution to the bulk stress is extracted (see measure_m3.py).
# Same maintained-shear setup as the Xu M1/M2 cases (opposite moving no-slip walls, Langmuir EOS) at a
# GENTLE Ca=0.1 (small-deformation rheology regime, stable), coverage X swept by MFC_SURF (0, 0.1, 0.3).
# Non-dim: Ca=mu*gdot*R/sigma0=0.1, Pe=gdot*R^2/D_s=10, Re=rho*gdot*R^2/mu=1 (finite-Re; Pimenta use
# Re=1e-2, unreachable by an explicit compressible solver -- see README). Box [-5,5]x[-2,2].
import json
import math
import os

W, Hbox = 10.0, 4.0
R = 1.0
Nx = int(os.environ.get("MFC_NX", "128"))
dx = W / Nx
Ny = int(Hbox / dx)
gdot = 1.0
Hwall = Hbox / 2

X = float(os.environ.get("MFC_SURF", "0.1"))
E_el = 0.2
surf_max = 1.0
Ca = 0.1
Pe = 10.0
rho0 = 0.2  # Re = rho*gdot*R^2/mu = 1
mu = 0.2
sigma0 = mu * gdot * R / Ca
D_s = gdot * R**2 / Pe
gam, p_inf, p0 = 2.0, 32.0, 8.0

sig_init = sigma0 * (1.0 + E_el * math.log(max(1.0 - X / surf_max, 1e-6)))
laplace = sig_init / R
c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
nu = mu / rho0
mydt = min(0.3 * dx / c_max, 0.2 * dx**2 / nu)
t_step_stop = round(5.0 / mydt)
t_step_save = max(1, t_step_stop // 60)
vx = f"{gdot:.4f}*y"

data = {
    "run_time_info": "T",
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -Hwall,
    "y_domain%end": Hwall,
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
    "weno_avg": "F",
    "weno_Re_flux": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -16,
    "bc_y%end": -16,
    "bc_y%vb1": -gdot * Hwall,
    "bc_y%ve1": gdot * Hwall,
    "num_patches": 2,
    "num_fluids": 2,
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "surfactant": "T",
    "surf_diff": D_s,
    "sigma_model": 3,
    "sigma_El": E_el,
    "surf_max": surf_max,
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
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Hbox,
    "patch_icpp(1)%vel(1)": vx,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha_rho(2)": 1.0e-8,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%surf_val": 0.0,
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": vx,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p0 + laplace,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8,
    "patch_icpp(2)%alpha_rho(2)": rho0,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%surf_val": X,
}
print(json.dumps(data))
