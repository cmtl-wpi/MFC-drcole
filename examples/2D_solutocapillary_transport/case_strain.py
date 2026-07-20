#!/usr/bin/env python3
# M0 STRAIN test (guide sec 3): the stretching term. A drop with UNIFORM surfactant sits in an
# extensional flow u=(eps*x, -eps*y) (divergence-free, so low-Mach friendly). The flow elongates the
# drop along x; the interfacial surfactant transport (surface convection - Gamma*div_s(u_s)) should
# sweep/redistribute surfactant, concentrating it toward the x-tips (elongation axis) and depleting the
# y-poles -- the classic Stone & Leal trend -- while TOTAL interfacial surfactant is conserved. This
# exercises the stretching term, which the conservative Gamma_tilde=Gamma|grad c| transport captures.
# The strain flow is imposed only as an IC (MFC cannot freeze it); run a moderate advective time (t~eps^-1).
import json
import os

R = 0.5
W = 4.0
Nx = int(os.environ.get("MFC_NX", "128"))
dx = W / Nx
epsr = 1.0                     # strain rate; Ca = mu*eps*R/sigma = 0.5 (moderate deformation)
gam, p_inf, p0, mu, rho0, sigma0 = 2.0, 32.0, 8.0, 0.1, 0.2, 0.1
laplace = sigma0 / R
vx = f"{epsr:.4f}*x"
vy = f"{-epsr:.4f}*y"

c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
mydt = 0.3 * dx / c_max
t_stop_phys = 0.6              # ~0.6/eps advective times
t_step_stop = round(t_stop_phys / mydt)
t_step_save = max(1, t_step_stop // 60)

data = {
    "run_time_info": "T", "x_domain%beg": -W/2, "x_domain%end": W/2, "y_domain%beg": -W/2, "y_domain%end": W/2,
    "m": Nx-1, "n": Nx-1, "p": 0, "cyl_coord": "F", "dt": mydt, "t_step_start": 0,
    "t_step_stop": t_step_stop, "t_step_save": t_step_save, "model_eqns": 3, "alt_soundspeed": "F",
    "mixture_err": "T", "mpp_lim": "F", "time_stepper": 3, "weno_order": 5, "weno_eps": 1e-16,
    "mapped_weno": "T", "null_weights": "F", "mp_weno": "T", "weno_avg": "T", "weno_Re_flux": "T",
    "riemann_solver": 2, "wave_speeds": 1, "avg_state": 2,
    "bc_x%beg": -3, "bc_x%end": -3, "bc_y%beg": -3, "bc_y%end": -3,   # non-reflecting (extensional in/outflow)
    "num_patches": 2, "num_fluids": 2, "viscous": "T", "surface_tension": "T", "sigma": sigma0,
    "surfactant": "T", "surf_diff": 0.0,
    "format": 1, "precision": 2, "prim_vars_wrt": "T", "cons_vars_wrt": "T", "cf_wrt": "T", "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0/(gam-1.0), "fluid_pp(1)%pi_inf": gam*p_inf/(gam-1.0), "fluid_pp(1)%cv": 1.0, "fluid_pp(1)%Re(1)": 1.0/mu,
    "fluid_pp(2)%gamma": 1.0/(gam-1.0), "fluid_pp(2)%pi_inf": gam*p_inf/(gam-1.0), "fluid_pp(2)%cv": 1.0, "fluid_pp(2)%Re(1)": 1.0/mu,
    "patch_icpp(1)%geometry": 3, "patch_icpp(1)%x_centroid": 0.0, "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W, "patch_icpp(1)%length_y": W, "patch_icpp(1)%vel(1)": vx, "patch_icpp(1)%vel(2)": vy,
    "patch_icpp(1)%pres": p0, "patch_icpp(1)%alpha_rho(1)": rho0, "patch_icpp(1)%alpha_rho(2)": 1.0e-8,
    "patch_icpp(1)%alpha(1)": 1.0, "patch_icpp(1)%alpha(2)": 0.0, "patch_icpp(1)%cf_val": 0.0, "patch_icpp(1)%surf_val": 0.0,
    "patch_icpp(2)%geometry": 2, "patch_icpp(2)%x_centroid": 0.0, "patch_icpp(2)%y_centroid": 0.0, "patch_icpp(2)%radius": R,
    "patch_icpp(2)%smoothen": "T", "patch_icpp(2)%smooth_patch_id": 1, "patch_icpp(2)%smooth_coeff": 0.5,
    "patch_icpp(2)%vel(1)": vx, "patch_icpp(2)%vel(2)": vy, "patch_icpp(2)%pres": p0+laplace,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8, "patch_icpp(2)%alpha_rho(2)": rho0, "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0, "patch_icpp(2)%cf_val": 1.0, "patch_icpp(2)%surf_val": 1.0,
}
print(json.dumps(data))
