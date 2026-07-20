#!/usr/bin/env python3
# FEASIBILITY TEST: can MFC hold solid-body rotation well enough for an exact-reference M0 transport check?
# Solid-body rotation u = (-w*y, w*x) is a steady Euler solution IF the pressure balances the centrifugal
# force: p = p0 + 0.5*rho*w^2*r^2. Matched density in both fluids -> no differential centrifugal force ->
# drop should stay circular. Inviscid: rigid rotation has ZERO strain rate, so viscosity is irrelevant.
# Surfactant is PASSIVE (sigma_model=0, surf_diff=0) -> pure advection. Exact solution on the interface:
#   Gamma(theta,t) = 1 + 0.5*cos(theta - w*t)   -> returns to the IC after one period T = 2*pi/w.
# The square box has NO boundary condition that exactly supports rotation (fluid crosses it), so we use
# ghost-cell extrapolation (-3, least constraining) and MEASURE the drift. That drift is the whole question.
import json
import os

W = 3.0  # [-1.5,1.5]^2
R = 0.6  # drop radius (R/dx ~ 26 at Nx=128)
Nx = int(os.environ.get("MFC_NX", "128"))
dx = W / Nx
Ny = Nx
w = 1.0  # angular velocity -> period T = 2*pi
rho0 = 1.0  # MATCHED density (both fluids) -> drop stays circular under rigid rotation
sigma0 = 0.5
gam, p_inf, p0 = 2.0, 100.0, 8.0  # stiff EOS -> c ~ 14.7 -> Mach ~ 0.14 at the corners

laplace = sigma0 / R
c_max = (gam * (p0 + 0.5 * rho0 * w**2 * (W / 2) ** 2 * 2 + laplace + p_inf) / rho0) ** 0.5
u_max = w * (W / 2) * 2**0.5  # corner speed
mydt = 0.3 * dx / (c_max + u_max)
T_rot = 2.0 * 3.141592653589793 / w
frac = float(os.environ.get("MFC_TFRAC", "1.0"))  # fraction of one rotation period to run
t_step_stop = round(frac * T_rot / mydt)
t_step_save = max(1, t_step_stop // 60)

# analytic ICs: rigid rotation + centrifugal-balanced pressure; cos(theta) surfactant pattern (regularized
# at the origin, which is interior to the drop and never sampled on the interface band)
vx = f"-{w:.4f}*y"
vy = f"{w:.4f}*x"
# NOTE: analytic exprs are parsed as PYTHON ast -> power is '**', not '^' ('^' is bitwise-xor and rejected)
pres_amb = f"{p0:.4f} + {0.5 * rho0 * w**2:.4f}*(x**2 + y**2)"
pres_drop = f"{p0 + laplace:.4f} + {0.5 * rho0 * w**2:.4f}*(x**2 + y**2)"
surf_pat = "1.0 + 0.5*x/sqrt(x**2 + y**2 + 0.0001)"  # eps only regularizes r=0 (interior); ~1e-4 err at r=R

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
    "weno_avg": "F",
    "weno_Re_flux": "F",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -3,
    "bc_y%end": -3,
    "num_patches": 2,
    "num_fluids": 2,
    "viscous": "F",
    "surface_tension": "T",
    "sigma": sigma0,
    "surfactant": "T",
    "surf_diff": 0.0,
    "sigma_model": 0,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": 1.0,
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": 1.0,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": W,
    "patch_icpp(1)%vel(1)": vx,
    "patch_icpp(1)%vel(2)": vy,
    "patch_icpp(1)%pres": pres_amb,
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
    "patch_icpp(2)%vel(2)": vy,
    "patch_icpp(2)%pres": pres_drop,
    "patch_icpp(2)%alpha_rho(1)": 1.0e-8,
    "patch_icpp(2)%alpha_rho(2)": rho0,
    "patch_icpp(2)%alpha(1)": 0.0,
    "patch_icpp(2)%alpha(2)": 1.0,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%surf_val": surf_pat,
}
print(json.dumps(data))
