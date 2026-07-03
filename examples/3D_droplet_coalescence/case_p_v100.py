#!/usr/bin/env python3
#
# V100-16GB max-grid variant of examples/3D_droplet_coalescence/case.py, case p.
#
# Identical physics and stretch geometry to the production case.py; ONLY the cell
# counts are reduced so the run fits on a single Tesla V100 (16 GB). Profiled on
# nighthawk 2026-07-02 (--gpu acc, 1 rank, case p config):
#   244x160x160 = 6.35M cells -> 15.75 GB used, 0.63 GB free (RC=0)   <- this grid
#   245x161x161 = 6.46M cells -> 16.13 GB used, 0.25 GB free (hard ceiling)
#   248x162x162               -> CUDA OOM
# 244x160x160 is the largest grid that keeps real headroom for a multi-day run.
#
# Fine-zone base resolution here ~D/77 (x), ~D/81 (yz) -- about half the linear
# resolution of the 45.2M-cell production grid (469x309x309, ~D/147-155); a single
# 16 GB V100 holds ~1/7 the production cell count. Expect the interface to be
# under-resolved relative to production -- this is a V100-scale run, not a
# converged one. See memory coalescence-v100-max-grid.

import argparse
import json
import math

CASES = {
    "a": (0.2, 14.8, 0.20, 120),
    "b": (0.5, 23.6, 0.10, 131),
    "c": (8.6, 105.9, 0.08, 153),
    "d": (15.2, 139.8, 0.08, 151),
    "e": (19.4, 158.0, 0.05, 151),
    "f": (32.8, 210.8, 0.08, 159),
    "g": (37.2, 228.0, 0.01, 164),
    "h": (61.4, 296.5, 0.06, 168),
    "i": (61.3, 295.3, 0.11, 167),
    "j": (56.3, 288.9, 0.13, 174),
    "k": (70.8, 327.7, 0.25, 178),
    "l": (48.1, 270.1, 0.39, 178),
    "m": (60.1, 302.8, 0.55, 179),
    "n": (65.1, 320.3, 0.49, 185),
    "o": (60.8, 313.7, 0.68, 190),
    "p": (64.9, 312.8, 0.71, 177),
    "q": (48.8, 260.3, 0.72, 163),
    "r": (14.5, 149.1, 0.34, 180),
}

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default="{}")
parser.add_argument("--case", type=str, default="p", choices=CASES.keys())
parser.add_argument("--n_start", type=int, default=0)
args = parser.parse_args()

We, Re, B, R_um = CASES[args.case]

rho_l = 763.0
rho_g = rho_l / 666
sigma = 0.0266
gamma_l = 3.7
gamma_g = 1.4

R = R_um * 1e-6
D = 2.0 * R

c_l_phys = 100

cg_over_cl = 1.0
p_gas = (cg_over_cl * c_l_phys) ** 2 * rho_g / gamma_g
p_ref = sigma / (0.72 * D)  # original reference pressure, retained for context
pi_inf_l = rho_l * c_l_phys**2 / gamma_l - p_gas
pi_inf_g = 0.0
p_liq = p_gas + 2 * sigma / R

Ur = math.sqrt(We * sigma / (rho_l * D))
Ud = Ur / 2.0

mu_l = rho_l * Ur * D / Re
mu_g = mu_l / 119

# Production stretch geometry (see case.py for the derivation). fine_* is the
# half-width of the uniform fine zone; the tanh sponge (a_x=8, a_yz=7) pushes the
# boundary out with few cells beyond it.
fine_x = 1.45 * D
fine_yz = 0.90 * D

x0, x1 = -1.60 * D, 1.60 * D
y0, y1 = -1.0 * D, 1.0 * D
z0, z1 = -1.0 * D, 1.0 * D

Nx = 244  # 245 cells over pre-stretch x = 3.2D  -> base ~D/77   (V100-16GB max)
Ny = 160  # 161 cells over pre-stretch y = 2.0D  -> base ~D/81
Nz = 160  # 161 cells over pre-stretch z = 2.0D  -> base ~D/81

eps = 1e-9

sep = 0.68 * D
b = B * D

t_end = 3e-3
t_save = t_end / 300  # 10 us cadence, 300 frames

data = {
    "run_time_info": "T",
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "stretch_x": "T",
    "a_x": 8.0,
    "x_a": -fine_x,
    "x_b": fine_x,
    "loops_x": 2,
    "stretch_y": "T",
    "a_y": 7.0,
    "y_a": -fine_yz,
    "y_b": fine_yz,
    "loops_y": 2,
    "stretch_z": "T",
    "a_z": 7.0,
    "z_a": -fine_yz,
    "z_b": fine_yz,
    "loops_z": 2,
    "cyl_coord": "F",
    "n_start": args.n_start,
    "cfl_adap_dt": "T",
    "cfl_target": 0.1,
    "t_stop": t_end,
    "t_save": t_save,
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
    "bc_y%beg": -8,
    "bc_y%end": -8,
    "bc_z%beg": -8,
    "bc_z%end": -8,
    "num_patches": 3,
    "num_fluids": 2,
    "surface_tension": "T",
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "vel_wrt(3)": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    "sigma": sigma,
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": 16.0 * D,
    "patch_icpp(1)%length_y": 12.0 * D,
    "patch_icpp(1)%length_z": 12.0 * D,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0,
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%x_centroid": -sep,
    "patch_icpp(2)%y_centroid": -b / 2.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.99,
    "patch_icpp(2)%vel(1)": Ud,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1,
    "patch_icpp(3)%geometry": 8,
    "patch_icpp(3)%x_centroid": sep,
    "patch_icpp(3)%y_centroid": b / 2.0,
    "patch_icpp(3)%z_centroid": 0.0,
    "patch_icpp(3)%radius": R,
    "patch_icpp(3)%alter_patch(1)": "T",
    "patch_icpp(3)%smoothen": "T",
    "patch_icpp(3)%smooth_patch_id": 1,
    "patch_icpp(3)%smooth_coeff": 0.99,
    "patch_icpp(3)%vel(1)": -Ud,
    "patch_icpp(3)%vel(2)": 0.0,
    "patch_icpp(3)%vel(3)": 0.0,
    "patch_icpp(3)%pres": p_liq,
    "patch_icpp(3)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(3)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(3)%alpha(1)": 1.0 - eps,
    "patch_icpp(3)%alpha(2)": eps,
    "patch_icpp(3)%cf_val": 1,
}

print(json.dumps(data))
