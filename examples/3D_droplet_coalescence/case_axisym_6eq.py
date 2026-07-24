#!/usr/bin/env python3
# 2D axisymmetric droplet coalescence -- 6-equation vs 5-equation film test.
#
# PURPOSE. The 3D 5-equation runs fuse on contact regardless of grid (case s, the
# deepest Regime II bounce point, still merged at D/200), which is why film AMR was
# being considered. Before spending effort on AMR plumbing, this case tests whether
# the INTERFACE MODEL is the limiter: the 6-equation model (model_eqns = 3) carries
# separate per-phase internal energies and relaxes them to a common pressure, rather
# than assuming pressure equilibrium pointwise as the 5-equation (Allaire) model does.
# If the gas film survives longer under 6-eq, the model is the lever; if it fuses at
# the same time, resolution/AMR is not the missing piece either and the answer is a
# film-drainage (disjoining-pressure) source term, which MFC does not have.
#
# RUN BOTH -- the comparison is the point, a single run answers nothing:
#   ./mfc.sh run examples/3D_droplet_coalescence/case_axisym_6eq.py -n 8 -- --model 5
#   ./mfc.sh run examples/3D_droplet_coalescence/case_axisym_6eq.py -n 8 -- --model 6
# Both use IDENTICAL numerics (see SCHEME note below), so any difference is the model.
#
# TWO IDEALIZATIONS, both deliberate -- neither reproduces case b exactly:
#   1. HEAD-ON. Axisymmetry forces impact parameter b = 0; case b has B = 0.10.
#      A 2D axisymmetric grid cannot represent an off-axis collision at all. Case b
#      is a defensible choice because its B is already small, but this run answers
#      "does 6-eq hold a film in a head-on collision", not "does it reproduce case b".
#      If you want the most axisymmetry-faithful case in the table, use --case g
#      (B = 0.01).
#   2. NO INTERFACE COMPRESSION. The 3D case runs MUSCL + THINC (int_comp = 1), but
#      int_comp is PROHIBITED under model_eqns = 3 ("THINC does not update per-fluid
#      internal energies, leaving thermodynamically inconsistent face states"). So
#      both models here run WENO5 with no compression. WENO5 is used rather than bare
#      MUSCL because losing THINC already costs interface sharpness and 5th order
#      recovers some of it. Expect a wider interface than the 3D runs either way --
#      compare 5-eq-here against 6-eq-here, NOT against the 3D results.

import argparse
import json
import math

# (We, Re, B, R_um) -- Qian & Law, same table as case.py
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
parser.add_argument("--model", type=int, default=6, choices=(5, 6), help="5-eq (model_eqns=2) or 6-eq (model_eqns=3)")
parser.add_argument("--npd", type=int, default=200, help="cells per droplet diameter in the film")
parser.add_argument("--n_start", type=int, default=0)
args = parser.parse_args()

We, Re, B, R_um = CASES[args.case]

# Fluid properties and the artificial-acoustics setup are copied verbatim from case.py
# so the 2D and 3D runs share an EOS. See case.py for why p_gas is raised to match the
# gas and liquid sound speeds (a near-vacuum ambient cavitates at impact and stalls dt).
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
pi_inf_l = rho_l * c_l_phys**2 / gamma_l - p_gas
pi_inf_g = 0.0
p_liq = p_gas + 2 * sigma / R

Ur = math.sqrt(We * sigma / (rho_l * D))
Ud = Ur / 2.0

mu_l = rho_l * Ur * D / Re
mu_g = mu_l / 119

sep = 0.68 * D
# b = B * D is DROPPED: axisymmetry requires a head-on collision (see header).

# Grid. Axial (x) is the collision axis, radial (y) starts at the axis.
#
# Axial: the droplets occupy 0.18D <= |x| <= 1.18D, so the uniform fine zone must
# cover 1.18D; the pre-stretch box is only slightly larger and a sharp ramp (a_x=16)
# keeps the stretch out of the droplet. Verified against the MFC tanh map:
#   reach +/-2.38D, film dx = D/npd, dx at the outer cap (1.18D) = D/134,
#   peak cell-to-cell growth 4.9% parked at 1.47D (outside the droplet).
#
# Radial: MFC's axisymmetric grid (grid_geometry = 2, y_domain%beg = 0) puts a
# HALF-WIDTH cell at the axis -- dy = y_end/(2n+1) for cell 0 and 2*dy for every cell
# beyond it. So the interior spacing is 2*y_end/(2n+1), which is what n is solved for
# below; the axis cell comes out at D/(2*npd), finer still. loops_y = 1 (not 3): the
# one-sided map at y_a = 0 overshoots badly with more loops (reach 5.8D, 20% growth).
#   reach 2.20D, interior dy = D/198, dy at 0.95D = D/200, growth 6.5% at 1.16D.
x_box = 1.40 * D
fine_x = 1.30 * D
y_end = 1.50 * D
fine_y = 1.15 * D

dx_target = D / args.npd
Nx = int(round(2 * x_box / dx_target))  # 560 cells at npd=200
# interior radial spacing is 2*y_end/(2n+1); solve for n
Ny = int(round(y_end / dx_target - 0.5)) + 1  # 301 cells at npd=200

eps = 1e-9

# 3.0 ms matches case.py: the default 1 ms is only ~1 tau and leaves the drop
# mid-stretch. Contact happens at t ~ (2*0.68 - 1)*D/Ur, well inside this window.
t_end = 3e-3
t_save = t_end / 300

# WENO5 for BOTH models -- int_comp is prohibited under model_eqns = 3, so the only
# way to compare the two models fairly is to run neither with compression.
scheme = {
    "recon_type": 1,
    "weno_order": 5,
    "weno_eps": 1e-6,
    "mapped_weno": "F",
    "wenoz": "T",
    "mp_weno": "T",
    "teno": "F",
}

model_eqns = 2 if args.model == 5 else 3

data = {
    "run_time_info": "T",
    "x_domain%beg": -x_box,
    "x_domain%end": x_box,
    "y_domain%beg": 0.0,
    "y_domain%end": y_end,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "T",  # 2D axisymmetric: x axial, y radial
    "stretch_x": "T",
    "a_x": 16.0,
    "x_a": -fine_x,
    "x_b": fine_x,
    "loops_x": 3,
    "stretch_y": "T",
    "a_y": 40.0,
    "y_a": 0.0,
    "y_b": fine_y,
    "loops_y": 1,
    "n_start": args.n_start,
    "cfl_adap_dt": "T",
    "cfl_target": 0.1,
    "t_stop": t_end,
    "t_save": t_save,
    "model_eqns": model_eqns,
    "alt_soundspeed": "F",  # 5-eq only; must stay F so both models share this deck
    "low_Mach": 0,
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": 3,
    **scheme,
    "avg_state": 2,  # required by model_eqns = 3
    "riemann_solver": 2,  # HLLC, required by model_eqns = 3
    "wave_speeds": 1,  # required by model_eqns = 3
    "viscous": "T",
    # Axial: nonreflecting subsonic outflow (as in the 3D case).
    "bc_x%beg": -8,
    "bc_x%end": -8,
    # Radial: reflective at the axis (y=0), nonreflecting outflow at y_max.
    "bc_y%beg": -2,
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
    # Patch 1: gas background. Oversized on purpose -- the post-stretch domain is
    # larger than the box above, and the rectangle is clipped to the grid.
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 3.0 * D,
    "patch_icpp(1)%length_x": 12.0 * D,
    "patch_icpp(1)%length_y": 12.0 * D,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0,
    # Patches 2/3: the two droplets, on the axis, approaching head-on.
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
    "patch_icpp(2)%cf_val": 1,
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
    "patch_icpp(3)%cf_val": 1,
}

# Case b at npd=200:  560 x 301 = 169k cells,  Ur = 0.258 m/s,  tau = D/Ur = 1.02 ms,
# contact at t ~ 0.37 ms,  dt ~ 1.3e-9 s (acoustic, c_l = 100 m/s) -> ~2.3M steps.
# No per-phase internal-energy IC is needed for the 6-eq model: pre_process sets the
# partial pressures to the mixture pressure (m_assign_variables.fpp).

print(json.dumps(data))
