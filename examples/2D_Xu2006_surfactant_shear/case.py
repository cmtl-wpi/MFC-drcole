#!/usr/bin/env python3
# Xu, Li, Lowengrub & Zhao (2006) single insoluble-surfactant drop in simple shear (guide Table 3.1).
# Domain [-4,4]x[-2,2], drop R=1 centered; simple shear u=(gammadot*y,0) via opposite moving no-slip
# walls; nonlinear Langmuir sigma(Gamma)=sigma0*(1+E*ln(1-Gamma/Gamma_inf)) with E=0.2. Coverage
# X=Gamma_initial/Gamma_inf swept via MFC_SURF (constant surf_val=X, Gamma_inf=surf_max=1). Non-dim
# groups: Ca=mu*gammadot*R/sigma0=0.7, Pe=gammadot*R^2/D_s=10, Re=rho*gammadot*R^2/mu kept LOW (~0.14,
# near the Stokes limit Xu uses) via low density; low-Mach via a stiff EOS. Gate: deformation D
# increases with coverage X; surfactant sweeps to the drop tips; sigma is minimal there.
import json
import os

W, Hbox = 8.0, 4.0  # [-4,4] x [-2,2]
R = 1.0
Nx = int(os.environ.get("MFC_NX", "128"))
dx = W / Nx
Ny = int(Hbox / dx)
gdot = 1.0
Hwall = Hbox / 2

X = float(os.environ.get("MFC_SURF", "0.1"))  # coverage Gamma_initial/Gamma_inf (0, 0.1, 0.3)
E_el = 0.2  # surfactant elasticity (Xu 2006)
surf_max = 1.0  # Gamma_inf (max packing)
# MFC is explicit-compressible: the Stokes limit (Re->0) needs a huge kinematic viscosity nu=mu/rho,
# whose viscous timestep ~ dx^2/nu vanishes, and very low rho blows up the acoustic CFL (c ~ 1/sqrt(rho)).
# So run the ACHIEVABLE regime: rho=0.2 (acoustically safe), Re=1 (finite-Re approx to Xu's Stokes),
# and Ca=0.4 (a 2D diffuse-interface drop over-deforms/breaks toward Ca=0.7 at finite Re). Faithful:
# geometry [-4,4]x[-2,2], R=1, Langmuir E=0.2, Pe=10. See README for these documented deviations.
rho0 = 0.2
mu = 0.2  # Re = rho*gdot*R^2/mu = 1
Ca = 0.3  # moderate deformation: sharp tips at higher Ca under-resolve and, with the
sigma0 = mu * gdot * R / Ca  # Langmuir saturation, destabilize the surfactant. See README.
D_s = 0.1  # Pe = gdot*R^2/D_s = 10
gam, p_inf, p0 = 2.0, 32.0, 8.0
import math

sig_init = sigma0 * (1.0 + E_el * math.log(max(1.0 - X / surf_max, 1e-6)))  # initial uniform sigma
laplace = sig_init / R
c_max = (gam * (p0 + p_inf) / rho0) ** 0.5  # ~20 -> Mach ~ gdot*R/c ~ 0.05
nu = mu / rho0
mydt = min(0.3 * dx / c_max, 0.2 * dx**2 / nu)  # min(acoustic, viscous) CFL
t_step_stop = round(5.0 / mydt)  # ~5 shear times (quasi-steady before tip instability)
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
