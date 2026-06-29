#!/usr/bin/env python3

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
parser.add_argument("--case", type=str, default="b", choices=CASES.keys())
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

# Ambient gas pressure.
# The previous setup used p_gas = sigma/(0.72*D) ~ 116 Pa for case f -- a
# near-vacuum that sits FAR below the collision dynamic pressure
# rho_l*Ur^2 ~ 2.7 kPa. At impact the tensile side of the collision drove the
# liquid to negative absolute pressure and cavitated (rho_l -> 0). The stiffened
# liquid EOS has c_l = sqrt(gamma_l*(p + pi_inf_l)/rho_l), which is UNBOUNDED as
# rho_l -> 0, so the local acoustic speed spiked to ~1.6e4 m/s (~166x c_l) and
# the adaptive-CFL controller drove dt -> 0 (run stalled near t ~ 1.5 D/Ur).
#
# Fix: raise the ambient so the gas acoustic speed matches the (artificial)
# liquid one, c_g = sqrt(gamma_g*p_gas/rho_g) = c_l. This is the HIGHEST ambient
# that does not tighten the acoustic dt -- the gas never limits below the liquid,
# so dt is unchanged from the healthy 2.12e-9 s -- while giving ~3x headroom over
# the dynamic pressure. If cavitation still appears, raise cg_over_cl above 1.0
# (cost: dt shrinks by ~1/cg_over_cl, since the gas then sets the CFL).
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

# Grid stretching (hyperbolic tangent, see src/pre_process/m_grid.f90).
# Strategy: define a SMALL fine domain at the target resolution, then let MFC
# expand the cells outside the fine zone outward (a sponge/buffer layer). The
# x_domain below is the PRE-stretch box; the actual simulated domain after
# stretching is larger. In the uniform fine zone, dx ~ (pre-stretch length)/N.
#
# Grid is 599 x 399 = 600 x 400 cells = 240k. The pre-stretch boxes are
# x = 3.0D over 600 cells and y = 2.0D over 400 cells, so the uniform spacing is
# D/200 in both directions; stretching makes the fine zone slightly finer and
# the outer buffer coarser. Empirically the fine-zone spacing is dx ~ dy ~ D/200
# (the initial acoustic dt = 0.1*dx/c_l scales with dx, so expect it ~1/2.2 of
# the D/90 value). NOTE: a diffuse interface is ~3-4 cells thick, so the
# inter-droplet gas film is now distinct down to ~6-8 cells ~ D/30, resolving
# the drainage/rupture phase substantially better than the prior D/90.
#
# fine_*  : half-width of the uniform fine zone the droplets live in (x_a/x_b).
# fine_x = 1.25D still encloses the initial droplet outer edge (centroid 0.68D +
# R 0.5D = 1.18D) with margin; the droplets then advect inward. The pre-stretch
# domain is set slightly wider than the fine zone to leave room for the stretch
# ramp before the outer (coarsening) buffer takes over.
fine_x = 1.25 * D
fine_y = 0.75 * D

x0, x1 = -1.5 * D, 1.5 * D
y0, y1 = -1.0 * D, 1.0 * D

Nx = 899  # 900 cells over x = 3.0D  -> D/300
Ny = 599  # 600 cells over y = 2.0D  -> D/300

eps = 1e-9

sep = 0.68 * D
b = B * D

t_end = 1e-3
t_save = t_end / 200

data = {
    "run_time_info": "T",
    "x_domain%beg": x0,
    "x_domain%end": x1,
    "y_domain%beg": y0,
    "y_domain%end": y1,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "stretch_x": "T",
    # Stretch factor maximized (was 2.0): a_x=2 left only ~10% of the declared
    # +/-fine_x zone truly uniform -- coarsening began almost immediately. a_x=40
    # makes ~85% of the fine zone genuinely D/200 while holding the cell-to-cell
    # growth ratio at the ~1.05 smoothness limit. Fine-zone dx (D/200) and extent
    # (x_a/x_b) are unchanged; only the buffer transition sharpens.
    "a_x": 40.0,
    "x_a": -fine_x,
    "x_b": fine_x,
    "loops_x": 2,
    "stretch_y": "T",
    "a_y": 40.0,
    "y_a": -fine_y,
    "y_b": fine_y,
    "loops_y": 2,
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
    "num_patches": 3,
    "num_fluids": 2,
    "surface_tension": "T",
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    "sigma": sigma,
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    # Background must cover the full STRETCHED grid, not just the pre-stretch
    # box; oversize generously (rectangle is clipped to the grid, droplets overwrite).
    "patch_icpp(1)%length_x": 16.0 * D,
    "patch_icpp(1)%length_y": 12.0 * D,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0,
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": -sep,
    "patch_icpp(2)%y_centroid": -b / 2.0,
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
    "patch_icpp(2)%cf_val": 1,
    "patch_icpp(3)%geometry": 2,
    "patch_icpp(3)%x_centroid": sep,
    "patch_icpp(3)%y_centroid": b / 2.0,
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
    "patch_icpp(3)%cf_val": 1,
}

print(json.dumps(data))
