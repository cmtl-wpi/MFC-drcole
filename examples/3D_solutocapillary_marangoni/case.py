#!/usr/bin/env python3
# 3D solutocapillary Marangoni flow driven by an insoluble interfacial surfactant.
#
# The 3D sibling of ../2D_solutocapillary_marangoni: a single spherical drop (fluid 2) at rest in
# an identical quiescent bath (fluid 1). No buoyancy, no density jump, no external forcing -- the
# drop's interface is seeded with a NON-UNIFORM insoluble surfactant coverage (linear in x, higher
# on the +x side). Surface tension falls where surfactant is high (sigma_model = 2:
# sigma(Gamma) = sigma + sigma_dGamma*Gamma, sigma_dGamma < 0), so a tangential gradient of sigma
# over the sphere drives a Marangoni flow from the surfactant-rich pole toward the clean pole.
# The surfactant is insoluble (transported along the interface, total conserved by the slip-wall
# box) and redistributes as the flow develops -- a transient 3D Marangoni convection visible as any
# nonzero max|u|, since the drop starts from rest.
#
# The interfacial surfactant is carried as the conserved smeared area-density Gamma*|grad c|; the
# surface concentration recovered for the closure is Gamma = (Gamma*|grad c|)/|grad c| on the
# interface band. That recovered Gamma is resolution-dependent (scales with the diffuse interface
# width), so sigma_dGamma is calibrated for a clear, stable effect on this grid rather than a
# specific physical coefficient. Surface (tangential) diffusion is not modeled (infinite surface
# Peclet). Keep `e`-notation out of analytic strings: MFC reads a bare `e` as Euler's number.

import json

# Geometry: D = 1 sphere centered in a 5D cube
D = 1.0
r = D / 2.0
W = 5.0 * D

Nx = 64  # cells per axis (64^3 ~ 2.6e5 cells); weno5 needs >= ~25 cells per rank-block
dx = W / Nx
Ny = Nx
Nz = Nx

# Stiffened-gas EOS: two identical fluids (gamma = 2); low-Mach, c ~ 20 at rho = 0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1
rho0 = 0.2

sigma0 = 0.1
sigma_dGamma = -0.05  # dsigma/dGamma < 0: surfactant lowers surface tension
Gamma0 = 1.0
laplace = 2.0 * sigma0 / r  # 3D Laplace jump (2*sigma/r for a sphere) suppresses the startup ring

# Non-uniform interfacial surfactant: more surfactant on the +x side of the drop.
surf_expr = f"{Gamma0:.6f}*(1.0 + 0.5*x)"

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
    # Slip walls on all six faces: closed box, total surfactant conserved
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
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
    # Patch 1: bath (3D cuboid spanning the domain), no surfactant, color c = 0
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
    # Patch 2: drop (sphere), color c = 1, non-uniform interfacial surfactant
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": r,
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
