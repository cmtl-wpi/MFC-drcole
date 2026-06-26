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
parser.add_argument("--case", type=str, default="f", choices=CASES.keys())
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
# Strategy: a uniform-ish D/150 fine zone sized to the region the interface
# actually sweeps, then an aggressively-stretched sponge that pushes the boundary
# far out with few cells. The x_domain below is the PRE-stretch box; the actual
# simulated domain after stretching is larger. In the fine zone, dx ~ length/N.
#
# This replaces the prior 96M-cell D/200 grid (599x399x399). That grid was
# inefficient: in 3D only ~23% of its cells sat on the droplet (0.65*0.60*0.60
# per-axis core fractions), with ~77% spent on far-field buffer, and its gentle
# a=2 stretch let resolution sag to ~D/108 right where the interface peaks
# (x=+/-1.17D, y/z=+/-0.69D, measured from run 2026-06-06_162725). Here we firm
# the core with a higher stretch factor (a_x=8, a_yz=7) and size the fine zone
# just past the interface reach, so the cells land where the physics is.
#
# Validated against the MFC tanh formula (470 x 310 x 310 = 45.2M cells, 2.1x
# fewer than 96M):
#   x : center D/147, dx@1.17D = D/112, domain +/-2.28D, peak growth 1.9%/cell
#   yz: center D/155, dx@0.69D = D/121, domain +/-1.47D, peak growth 2.7%/cell
# Interface resolution (D/112-121) is thus slightly FINER than the old D/108 at
# ~1/2.9 the compute cost (fewer cells x larger dt, dx_min ~D/150 vs D/200).
#
# fine_*  : half-width of the uniform fine zone the droplets live in (x_a/x_b),
# set ~0.25D beyond the interface reach so the stretch ramp stays outside the
# droplet. Keep cell-to-cell growth <~5%; SMOKE-TEST the sponge (pre_process +
# a few hundred steps, check symmetry_violation / pressure field) before a full
# run, since the buffer geometry differs from the proven D/200 grid.
fine_x = 1.45 * D
fine_yz = 0.90 * D

x0, x1 = -1.60 * D, 1.60 * D
y0, y1 = -1.0 * D, 1.0 * D
z0, z1 = -1.0 * D, 1.0 * D

Nx = 469  # 470 cells over pre-stretch x = 3.2D  -> base ~D/147
Ny = 309  # 310 cells over pre-stretch y = 2.0D  -> base ~D/155
Nz = 309  # 310 cells over pre-stretch z = 2.0D  -> base ~D/155

eps = 1e-9

sep = 0.68 * D
b = B * D

# 3.0 ms: case f needs the full Qian & Law timeline (their last frame is 2.88 ms);
# +0.12 ms past it covers the sim/experiment contact offset (~0.035 ms) with margin.
# Grid adequacy verified against the D/90 3 ms run (2026-06-02_140913_case-f):
# over 0-3 ms the interface stays within x=+/-1.18D and y=+/-0.90D (the drop
# coalesces, n_fragments=1, and OSCILLATES rather than stretching/separating), so
# the fine zone above (fine_x=1.45D, fine_yz=0.90D) already contains it -- x with
# 0.27D margin, y exactly at the fine-zone edge (still resolved, zero spare margin).
t_end = 3e-3
t_save = t_end / 300  # 10 us cadence, 300 frames (resolves the Qian & Law frames)

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
    # Background must cover the full STRETCHED grid, not just the pre-stretch
    # box; oversize generously (cuboid is clipped to the grid, droplets overwrite).
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
