#!/usr/bin/env python3
# M0 (guide sec 3): interface-confined surfactant transport under a KNOWN velocity, decoupled from the
# flow solver. Uniform translation u=(U,0) is exactly maintained in a compressible solver (Galilean),
# so it is the cleanest prescribed-velocity test in MFC. A passively-coated drop (sigma_model=0, no
# Marangoni) with a cos(theta) surfactant pattern translates once around the periodic domain; the
# surfactant should ride rigidly with the interface and stay confined (no bulk smearing). surf_diff=0
# => pure advection, isolating advection-driven leakage (the question that decides whether the Jain
# sharpening flux must be extended to advection). Drop starts at x=-1 and returns there at t=W/U, so
# t=0 and the final frame are directly comparable.
import json
import os

R = 0.5
W = 4.0
Nx = int(os.environ.get("MFC_NX", "128"))
dx = W / Nx
xc0 = -1.0                     # initial drop-center x
U = 1.0                        # translation speed
gam, p_inf, p0, mu, rho0, sigma0 = 2.0, 32.0, 8.0, 0.1, 0.2, 0.1
laplace = sigma0 / R
# cos(theta) pattern about the drop center: Gamma = 1 + amp*(x-xc0)/R (=1+amp*cos(theta) on the band)
amp = 0.5
surf_expr = f"1.0 + {amp / R:.6f}*(x - ({xc0}))"

c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
mydt = 0.3 * dx / c_max        # acoustic CFL (no diffusion, no viscous strain in uniform flow)
t_cross = W / U                # one full periodic lap -> drop back at xc0
t_step_stop = round(t_cross / mydt)
t_step_save = max(1, t_step_stop // 60)

data = {
    "run_time_info": "T", "x_domain%beg": -W/2, "x_domain%end": W/2, "y_domain%beg": -W/2, "y_domain%end": W/2,
    "m": Nx-1, "n": Nx-1, "p": 0, "cyl_coord": "F", "dt": mydt, "t_step_start": 0,
    "t_step_stop": t_step_stop, "t_step_save": t_step_save, "model_eqns": 3, "alt_soundspeed": "F",
    "mixture_err": "T", "mpp_lim": "F", "time_stepper": 3, "weno_order": 5, "weno_eps": 1e-16,
    "mapped_weno": "T", "null_weights": "F", "mp_weno": "T", "weno_avg": "T", "weno_Re_flux": "T",
    "riemann_solver": 2, "wave_speeds": 1, "avg_state": 2,
    "bc_x%beg": -1, "bc_x%end": -1, "bc_y%beg": -1, "bc_y%end": -1,   # periodic
    "num_patches": 2, "num_fluids": 2, "viscous": "T", "surface_tension": "T", "sigma": sigma0,
    "surfactant": "T", "surf_diff": 0.0,
    "format": 1, "precision": 2, "prim_vars_wrt": "T", "cons_vars_wrt": "T", "cf_wrt": "T", "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0/(gam-1.0), "fluid_pp(1)%pi_inf": gam*p_inf/(gam-1.0), "fluid_pp(1)%cv": 1.0, "fluid_pp(1)%Re(1)": 1.0/mu,
    "fluid_pp(2)%gamma": 1.0/(gam-1.0), "fluid_pp(2)%pi_inf": gam*p_inf/(gam-1.0), "fluid_pp(2)%cv": 1.0, "fluid_pp(2)%Re(1)": 1.0/mu,
    # Patch 1: bath (full box), moving at U, no surfactant
    "patch_icpp(1)%geometry": 3, "patch_icpp(1)%x_centroid": 0.0, "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W, "patch_icpp(1)%length_y": W, "patch_icpp(1)%vel(1)": U, "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0, "patch_icpp(1)%alpha_rho(1)": rho0, "patch_icpp(1)%alpha_rho(2)": 1.0e-8,
    "patch_icpp(1)%alpha(1)": 1.0, "patch_icpp(1)%alpha(2)": 0.0, "patch_icpp(1)%cf_val": 0.0, "patch_icpp(1)%surf_val": 0.0,
    # Patch 2: drop at (xc0,0), moving at U, cos(theta) surfactant
    "patch_icpp(2)%geometry": 2, "patch_icpp(2)%x_centroid": xc0, "patch_icpp(2)%y_centroid": 0.0, "patch_icpp(2)%radius": R,
    "patch_icpp(2)%smoothen": "T", "patch_icpp(2)%smooth_patch_id": 1, "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": U, "patch_icpp(2)%vel(2)": 0.0, "patch_icpp(2)%pres": p0+laplace,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8, "patch_icpp(2)%alpha_rho(2)": rho0, "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0, "patch_icpp(2)%cf_val": 1.0, "patch_icpp(2)%surf_val": surf_expr,
}
print(json.dumps(data))
