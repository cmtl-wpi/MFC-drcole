#!/usr/bin/env python3
# 3D thermocapillary droplet migration -- Samareh, Mostaghimi & Moreau, Int. J. Heat Mass Transfer
# 73 (2014) 616-626, Sec. 4.1.1 / Fig. 6: a fully 3D SPHERE in an imposed linear temperature
# gradient, target v_t/v_YGB ~ 0.95. This is the 3D sibling of ../2D_thermocapillary_migration,
# which reproduces Samareh's Fig 5 (the planar 2D cylinder -> ~0.80). Same physics and EOS
# realization; the only changes are the third dimension (a sphere patch, p > 0) and -- crucially --
# BULK CONDUCTION IS MANDATORY here.
#
# Why conduction is required in 3D (and optional in 2D): in the no-conduction (frozen-T) limit the
# 2D rise reaches a quasi-steady plateau, but the 3D rise does NOT -- the toroidal internal
# circulation continuously steepens the frozen interfacial gradient, so the velocity drifts past
# v_YGB without saturating (finer grid -> faster), leaving no validatable 3D number. Bulk Fourier
# conduction (the thermal_conduction feature, set by Ma > 0) diffuses the temperature, tames that
# runaway, and restores a steady plateau comparable to Samareh's 0.95. So this case runs with
# thermal_conduction = T and isothermal gradient y-walls by construction.
#
# Like the 2D case, T is imposed through density: at uniform pressure the EOS gives T =
# (p+p_inf)/((gam-1)*rho*cv), so rho(y) = rho_coeff/T(y) encodes the field; a finite conductivity
# (small Ma) holds T near the imposed gradient against the drop's advection. T is shifted up by
# T0 = 10 (Samareh's T0 = 0) so rho stays positive; only gradT and sigma_T (which set v_YGB) are
# physical. Keep `e`-notation out of analytic strings: MFC's IC parser reads `e` as Euler's number.

import json

# Marangoni number: conduction is REQUIRED in 3D (the frozen-T 3D rise has no plateau), so Ma > 0
# here by construction (smaller -> closer to Samareh's invariant-T limit, but dt ~ Ma so longer).
Ma = 1.0

# Geometry: D = 1 sphere, 5D x 5D wide x 7.5D tall box; drop 1.5D above the cold floor, rise axis = y
D = 1.0
r = D / 2.0
W = 5.0 * D
Ly = 7.5 * D
y_drop = -Ly / 2.0 + 1.5 * D

Nx = 64  # cells across the box width (Samareh 3D used 64, 128)
dx = W / Nx
Ny = round(Ly / dx)  # cells along the 7.5D rise axis (= 1.5*Nx)
Nz = Nx  # cube cross-section (W x W)

# Stiffened-gas EOS: two identical fluids (gamma = 2); low-Mach, c ~ 20 at rho = 0.2
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1  # dynamic viscosity of both phases (MFC takes Re = 1/mu)

# Imposed linear field T(y) = T0 + gradT*y, encoded as rho(y) = rho_coeff/T(y) at uniform p0
T0 = 10.0
gradT = 2.0 / 15.0  # |dT/dy| = 0.1333 (Samareh)
sigma0 = 0.1
sigma_T = -0.1  # dsigma/dT
G = abs(sigma_T * gradT)  # Marangoni stress scale = 0.013333
rho_drop = 0.2  # Samareh: rho_d = rho_b = 0.2 (density at the drop's reference temperature T0)
rho_coeff = rho_drop * T0
cv = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # closes the EOS so rho(T0) = rho_drop
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop's initial position
v_YGB = (2.0 / 15.0) * (-sigma_T) * gradT * r / mu  # YGB terminal speed = 8.889e-3 (mu* = k* = 1)

# Conduction properties from Ma: smaller Ma = larger diffusivity = T held closer to the gradient.
# k = alpha*rho*cp at the reference state; the local diffusivity k/(rho*cp) varies mildly with rho.
alpha_b = G * r**2 / (mu * Ma)  # bulk thermal diffusivity
cp = gam * cv
k_b = alpha_b * rho_drop * cp  # bulk conductivity

eps = 1.0e-9  # trace volume fraction of the (identical) second phase

# Two patches share the analytic linear-T density rho(y) = rho_coeff/T(y); only the color function
# differs (background c = 0, sphere c = 1), so the capillary stress acts purely on the color
# interface (the mu* = k* = 1 YGB limit, no real density jump). Fold eps into a plain decimal --
# never embed a bare `e` (1e-9) token in an analytic string.
rho_num = (1.0 - eps) * rho_coeff
rho_expr = f"{rho_num:.9f}/({T0} + {gradT:.9f}*y)"  # ~ 2.0/(10.0 + 0.133333333*y)

# Time stepping: min(acoustic CFL, explicit-conduction limit). 3D diffusion number uses d = 3.
rho_min = rho_coeff / (T0 + gradT * Ly / 2.0)  # hot wall: lowest density, max sound speed
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
t_r = mu / G  # capillary-thermal time = 7.5
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (6.0 * alpha_b))  # 3D explicit diffusion number (dt <= 0.35*dx^2/(2*d*alpha), d = 3)
t_step_stop = round(2.0 * t_r / mydt)  # 2 capillary-thermal times (enough to reach the plateau)
t_step_save = max(1, t_step_stop // 60)

# Configuration case dictionary
data = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain: rise (gradient) axis y in [-Ly/2, Ly/2]; short axes x,z in [-W/2, W/2]
    "x_domain%beg": -W / 2,
    "x_domain%end": W / 2,
    "y_domain%beg": -Ly / 2,
    "y_domain%end": Ly / 2,
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
    # Simulation Algorithm
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
    # Slip walls on all six faces (Samareh's box); isothermal gradient walls on y set below
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + temperature-dependent surface tension sigma(T)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    # Database Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Fluid 1 (continuous phase)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(1)%k_therm": k_b,
    # Fluid 2 (identical properties)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%k_therm": k_b,
    # Patch 1: background medium (3D cuboid spanning the domain). Analytic linear-T density, color c = 0.
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": W,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": W,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    # Patch 2: droplet (3D sphere), color c = 1 smeared over ~2 cells. Identical density/composition/
    # pressure to patch 1, so the capillary stress acts purely on the color interface.
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": r,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.5,  # diffuse color interface ~2 cells wide (half-width w = dx/coeff)
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p0,
    "patch_icpp(2)%alpha_rho(1)": rho_expr,
    "patch_icpp(2)%alpha_rho(2)": eps,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1.0,
    # Isothermal Dirichlet gradient walls pin T to the imposed gradient (cold floor / hot ceiling).
    # This is what tames the frozen-T 3D runaway into a validatable plateau.
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T0 + gradT * (-Ly / 2.0),
    "bc_y%Twall_out": T0 + gradT * (Ly / 2.0),
}

print(json.dumps(data))
