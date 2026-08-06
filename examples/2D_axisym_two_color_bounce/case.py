#!/usr/bin/env python3
# Axisymmetric head-on droplet collision in the bouncing regime (Qian & Law
# case b scale: We = 0.5, R = 131 um), with each droplet carrying its own
# color-function marker (num_colors = 2) so the interfaces keep full surface
# tension in the gas film and cannot numerically coalesce. Geometry follows
# the 2D_axisym_droplet_coalescence example: x is axial, y is radial,
# bc_y%beg = -2 is the axis, and a circle patch on the axis is a sphere.
# The liquid sound speed is artificial (100 m/s, gas matched) to keep the
# acoustic time step affordable -- this is a low-Mach marker demo, not a
# compressible-film study.

import json
import math

We = 0.5
Re = 23.6

rho_l = 763.0
rho_g = rho_l / 666
sigma = 0.0266
gamma_l = 3.7
gamma_g = 1.4

R = 131e-6
D = 2.0 * R

c_l = 100.0
p_gas = c_l**2 * rho_g / gamma_g  # gas sound speed matched to the liquid one
pi_inf_l = rho_l * c_l**2 / gamma_l - p_gas
p_liq = p_gas + 2 * sigma / R

Ur = math.sqrt(We * sigma / (rho_l * D))
Ud = Ur / 2.0

mu_l = rho_l * Ur * D / Re
mu_g = mu_l / 119

eps = 1e-9
sep = 0.68 * D

t_c = D / Ur
t_end = 2.5 * t_c

data = {
    "run_time_info": "T",
    "x_domain%beg": -1.5 * D,
    "x_domain%end": 1.5 * D,
    "y_domain%beg": 0.0,
    "y_domain%end": 1.0 * D,
    "m": 239,
    "n": 79,
    "p": 0,
    "cyl_coord": "T",
    "n_start": 0,
    "cfl_adap_dt": "T",
    "cfl_target": 0.2,
    "t_stop": t_end,
    "t_save": t_end / 100,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "low_Mach": 0,
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": 3,
    "recon_type": 2,
    "muscl_order": 2,
    "muscl_lim": 4,
    "int_comp": 1,
    "ic_beta": 1.6,
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "bc_x%beg": -8,
    "bc_x%end": -8,
    "bc_y%beg": -2,
    "bc_y%end": -8,
    "num_patches": 3,
    "num_fluids": 2,
    "surface_tension": "T",
    "num_colors": 2,
    "sigma": sigma,
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    # background gas: both colors zero
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.5 * D,
    "patch_icpp(1)%length_x": 4.0 * D,
    "patch_icpp(1)%length_y": 3.0 * D,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val(1)": 0,
    "patch_icpp(1)%cf_val(2)": 0,
    # left droplet (sphere on the axis): marker 1
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": -sep,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.99,
    "patch_icpp(2)%vel(1)": Ud,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val(1)": 1,
    "patch_icpp(2)%cf_val(2)": 0,
    # right droplet: marker 2
    "patch_icpp(3)%geometry": 2,
    "patch_icpp(3)%x_centroid": sep,
    "patch_icpp(3)%y_centroid": 0.0,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%smoothen": "T",
    "patch_icpp(3)%smooth_patch_id": 1,
    "patch_icpp(3)%smooth_coeff": 0.99,
    "patch_icpp(3)%vel(1)": -Ud,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%pres": p_liq,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(3)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%cf_val(1)": 0,
    "patch_icpp(3)%cf_val(2)": 1,
}

print(json.dumps(data))
