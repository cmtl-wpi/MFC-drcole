#!/usr/bin/env python3
# 2D solutocapillary Marangoni flow driven by an insoluble interfacial surfactant.
#
# A single drop (fluid 2) sits at rest in an otherwise quiescent bath (fluid 1); the two
# phases are identical (same EOS, density, viscosity), so there is no buoyancy, no density
# jump, and no external forcing -- the ONLY thing that can move the fluid is the surfactant.
# The drop is seeded with a NON-UNIFORM interfacial surfactant coverage (linear in x, higher
# on the +x side). Surface tension falls where surfactant is high (sigma_model = 2:
# sigma(Gamma) = sigma + sigma_dGamma*Gamma, sigma_dGamma < 0), so a tangential gradient of
# sigma along the interface appears and drives a Marangoni flow from the low-sigma (surfactant-
# rich) pole toward the high-sigma (clean) pole. The surfactant is insoluble: it is transported
# along the interface and its total is conserved (slip walls close the box), redistributing as
# the flow develops -- a transient Marangoni convection that any nonzero max|u| makes visible,
# since the drop starts from rest.
#
# The interfacial surfactant is carried as the conserved smeared area-density Gamma*|grad c|;
# the surface concentration recovered for the closure is Gamma = (Gamma*|grad c|)/|grad c| on
# the interface band. This recovered Gamma is resolution-dependent (it scales with the diffuse
# interface width), so sigma_dGamma here is calibrated to give a clear, stable effect on this
# grid rather than to match a specific physical coefficient. Surface (tangential) diffusion of
# the surfactant is not modeled (infinite surface Peclet limit).
#
# Keep `e`-notation out of analytic strings: MFC's IC parser reads a bare `e` as Euler's number.

import json

# Geometry: D = 1 drop centered in a 5D x 5D box
D = 1.0
r = D / 2.0
W = 5.0 * D

Nx = 128  # cells across the box (weno5 needs >= ~25 cells per rank-block, so <= 4 ranks per axis)
dx = W / Nx
Ny = Nx

# Stiffened-gas EOS: two identical fluids (gamma = 2); low-Mach, c ~ 20 at rho = 0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1  # dynamic viscosity of both phases (MFC takes Re = 1/mu)
rho0 = 0.2

# Surface tension: constant base sigma0 plus a linear solutocapillary correction sigma(Gamma).
sigma0 = 0.1
sigma_dGamma = -0.05  # dsigma/dGamma < 0: surfactant lowers surface tension
Gamma0 = 1.0  # reference interfacial surfactant coverage
laplace = sigma0 / r  # Laplace pressure jump, initialized to suppress the startup acoustic ring

# Non-uniform interfacial surfactant: more surfactant on the +x side of the drop.
surf_expr = f"{Gamma0:.6f}*(1.0 + 0.5*x)"

# Time stepping: acoustic-CFL limited; run long enough to develop a clear Marangoni flow.
c_max = (gam * (p0 + p_inf) / rho0) ** 0.5
mydt = 0.3 * dx / c_max
t_step_stop = 2000
t_step_save = t_step_stop // 50

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
    # Simulation algorithm
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
    # Slip walls on all four faces: closed box, total surfactant conserved
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + solutocapillary surface tension sigma(Gamma)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "surfactant": "T",
    "sigma_model": 2,
    "sigma_dGamma": sigma_dGamma,
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluids (identical)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": 1.0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": 1.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    # Patch 1: bath (full box), no surfactant, color c = 0
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
    # Patch 2: drop (circle), color c = 1, non-uniform interfacial surfactant
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": r,
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
