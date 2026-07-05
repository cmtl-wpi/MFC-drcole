#!/usr/bin/env python3
# Film-zoom validation probe - axisymmetric head-on droplet vs symmetry plane.
#
# Purpose (docs/gas_film_drainage_bouncing/DESIGN.md): produce an in-house DNS
# reference for the sub-grid film model's DIMPLE / EARLY-DRAINAGE phase, by
# resolving the gas film down to h ~ 2-3 cells (~0.3-0.45 um at the pilot
# resolution) between a deforming drop and the collision symmetry plane.
# NOT a rupture-resolving run: continuum breaks below ~66 nm (Kn > 1) and the
# validation window deliberately ends when h_min ~ 2-3 dx (expect NUMERICAL
# merge with the mirror image shortly after - that event is recorded, not
# trusted). Targets: Chen & Yang (2020) Figs. 13/15 shapes, Klaseboer (2000)
# h_min(t) (prior Phase-1 goal, +/-15%), and the DESIGN.md compliance closure.
#
# Physics: Chen & Yang Fig. 6 head-on bounce point - We=8.6, Re=105.9, B=0,
# D=306 um, tetradecane in N2. EOS family identical to the 2D case_s runs
# (ARTIFICIAL liquid sound speed c_l=100 m/s for affordability; ambient raised
# so c_g ~ c_l - see examples/2D_droplet_coalescence/runs/case_s_nf3/case.py).
#
# Artificial-Mach caveat (squeeze-film radial outflow, u_r ~ a*V_n/(2h)):
#     a=20um V_n=0.50 h=0.50um -> u_r ~ 10 m/s   (early film)
#     a=40um V_n=0.30 h=0.30um -> u_r ~ 20 m/s   (mid drainage)
#     a=60um V_n=0.10 h=0.20um -> u_r ~ 15 m/s   (late, decelerated)
# vs c_artificial = 100 m/s: worst-case artificial Ma ~ 0.2. GATE-3: monitor
# max gas |u| in the film band; ABORT the validation window if it exceeds
# 30 m/s (artificial Ma > 0.3 - compressibility errors ~Ma^2 no longer small).
#
# Modes:
#   --mode probe  : moving drop (Ud = Ur/2 toward the plane), t_stop 4e-4 s
#   --mode static : GATE-1 - drop AT REST, surface 1 um from the plane, 5e-5 s.
#                   Pass criteria: (i) Laplace jump p_drop - p_gas = 2*sigma/R
#                   = 347.7 Pa +/- 10% (this MEASURES whether MFC's surface
#                   tension carries the axisymmetric geometric terms - no
#                   shipped example validates cyl_coord + surface_tension);
#                   (ii) max spurious |u| < 0.025 m/s (5% of Ud), not growing.
#   --check       : print derived quantities + internal consistency checks, exit.
#
# Cost (pilot dx=150nm, cfl 0.1, c=100): dt ~ 1.5e-10 s, 2.7e6 steps to 0.4 ms,
# ~6-11 h on 64 ranks at 278k cells. Production dx=75nm: ~5.4e6 steps, 599k
# cells. Restart/extend: --n_start <step> after raising --t-stop.

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default="{}")
parser.add_argument("--mode", type=str, default="probe", choices=["probe", "static"])
parser.add_argument("--dx-fine", type=float, default=150e-9, dest="dx_fine")
parser.add_argument("--cfl", type=float, default=0.1)
parser.add_argument("--t-stop", type=float, default=None, dest="t_stop")
parser.add_argument("--n_start", type=int, default=0)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

# Physics - Chen & Yang Fig. 6 point, case_s EOS family
We, Re = 8.6, 105.9
R = 153e-6
D = 2.0 * R

rho_l = 763.0
rho_g = rho_l / 666.0
sigma = 0.0266
gamma_l = 3.7
gamma_g = 1.4
c_l = 100.0  # artificial liquid sound speed (affordability; see header)

p_gas = c_l**2 * rho_g / gamma_g
pi_inf_l = rho_l * c_l**2 / gamma_l - p_gas
p_liq = p_gas + 2.0 * sigma / R

Ur = math.sqrt(We * sigma / (rho_l * D))  # relative (mirror-pair) velocity
Ud = Ur / 2.0  # single-drop speed toward the plane
mu_l = rho_l * Ur * D / Re
mu_g = mu_l / 119.0

# Geometry per mode
if args.mode == "probe":
    gap0 = 20e-6  # drop surface to plane: inertial approach + film formation
    vel_x = -Ud
    t_stop = args.t_stop if args.t_stop is not None else 4.0e-4
    t_save = 2.0e-6  # 200 frames
else:  # static (GATE-1)
    gap0 = 1e-6  # surface 1 um from the plane, inside the fine band
    vel_x = 0.0
    t_stop = args.t_stop if args.t_stop is not None else 5.0e-5
    t_save = 2.5e-6  # 20 frames
x0 = R + gap0  # drop centroid (axial); drop on the axis (y = 0)

# Grid - parameters verified by grid_check.py (MFC tanh stretch replicated).
# Axial x: half-domain [0, Lx_pre], symmetry plane at x=0. Slope of the stretch
# map at x=0 is exactly 1, so dx(plane) = Lx_pre/(m+1) = dx_fine. Fine band
# +/-4 um; two loops expand the far field to ~1.5 mm while the drop stays
# resolved at <= 3 um. Radial y: axisym half-cell-offset grid, interior
# dy = 2*Ly_pre/(2n+1) = 0.50 um uniform through the film-disk/rim zone
# (r_h grows to ~0.2D = 61 um), one loop coarsens to the far field.
Lx_pre = 80e-6
m = int(round(Lx_pre / args.dx_fine)) - 1
Ly_pre = 260e-6
n = 520

data = {
    "run_time_info": "T",
    # grid
    "x_domain%beg": 0.0,
    "x_domain%end": Lx_pre,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly_pre,
    "m": m,
    "n": n,
    "p": 0,
    "cyl_coord": "T",  # x axial, y radial (2D axisymmetric)
    "stretch_x": "T",
    "a_x": 10.0,
    "x_a": -4e-6,
    "x_b": 4e-6,
    "loops_x": 2,
    "stretch_y": "T",
    "a_y": 8.0,
    "y_a": -90e-6,
    "y_b": 90e-6,
    "loops_y": 1,
    # boundary conditions
    "bc_x%beg": -2,  # collision symmetry plane (free-slip reflective)
    "bc_x%end": -6,  # far field (shipped-axisym-example precedent)
    "bc_y%beg": -2,  # axis (validator-required for 2D cylindrical)
    "bc_y%end": -6,
    # time
    "n_start": args.n_start,
    "cfl_adap_dt": "T",
    "cfl_target": args.cfl,
    "t_stop": t_stop,
    "t_save": t_save,
    # model / numerics (case_s family)
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
    "viscous": "T",  # MANDATORY - film drainage IS gas viscosity
    "num_patches": 2,
    "num_fluids": 2,
    "surface_tension": "T",
    "sigma": sigma,
    # output
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    # probes: pressure time series at the plane center + three rim radii
    "fd_order": 2,  # required by probe_wrt
    "probe_wrt": "T",
    "num_probes": 4,
    "probe(1)%x": 0.5e-6,
    "probe(1)%y": 0.5e-6,
    "probe(2)%x": 0.5e-6,
    "probe(2)%y": 20e-6,
    "probe(3)%x": 0.5e-6,
    "probe(3)%y": 40e-6,
    "probe(4)%x": 0.5e-6,
    "probe(4)%y": 60e-6,
    # fluids
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    # patch 1: nitrogen background (oversized; clipped to the stretched grid)
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 8e-3,
    "patch_icpp(1)%length_y": 8e-3,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": 1e-9 * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - 1e-9) * rho_g,
    "patch_icpp(1)%alpha(1)": 1e-9,
    "patch_icpp(1)%alpha(2)": 1.0 - 1e-9,
    "patch_icpp(1)%cf_val": 0,
    # patch 2: tetradecane drop on the axis
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": x0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.99,
    "patch_icpp(2)%vel(1)": vel_x,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - 1e-9) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": 1e-9 * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - 1e-9,
    "patch_icpp(2)%alpha(2)": 1e-9,
    "patch_icpp(2)%cf_val": 1,
}

if args.check:
    dx0 = Lx_pre / (m + 1)
    dy_int = 2.0 * Ly_pre / (2 * n + 1)
    dt = args.cfl * min(dx0, dy_int) / (c_l + 1.0)
    checks = {
        "We recovered": abs(rho_l * Ur**2 * D / sigma - We) < 1e-9,
        "Re recovered": abs(rho_l * Ur * D / mu_l - Re) < 1e-9,
        "alphas sum patch1": abs(data["patch_icpp(1)%alpha(1)"] + data["patch_icpp(1)%alpha(2)"] - 1.0) < 1e-15,
        "alphas sum patch2": abs(data["patch_icpp(2)%alpha(1)"] + data["patch_icpp(2)%alpha(2)"] - 1.0) < 1e-15,
        "dx0 == dx_fine": abs(dx0 - args.dx_fine) / args.dx_fine < 0.01,
        "static drop in fine band": (args.mode != "static") or (x0 - R <= 4e-6),
        "drop moves toward plane": (args.mode != "probe") or (vel_x < 0),
        "viscous on": data["viscous"] == "T",
        "symmetry plane": data["bc_x%beg"] == -2,
        "axis bc": data["bc_y%beg"] == -2,
    }
    print(f"mode={args.mode}  We={We} Re={Re}  Ur={Ur:.5f} Ud={Ud:.5f} m/s")
    print(f"mu_l={mu_l:.4e} mu_g={mu_g:.4e}  sigma={sigma}  2sigma/R={2*sigma/R:.1f} Pa")
    print(f"p_gas={p_gas:.1f} p_liq={p_liq:.1f} Pa  pi_inf_l={pi_inf_l:.4e}")
    print(f"grid m={m} n={n} cells={(m+1)*(n+1)}  dx0={dx0*1e9:.1f}nm dy_int={dy_int*1e6:.2f}um")
    print(f"dt~{dt:.3e}s  steps to t_stop~{t_stop/dt:.3e}  x0={x0*1e6:.1f}um gap0={gap0*1e6:.1f}um")
    bad = [k for k, v in checks.items() if not v]
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    sys.exit(1 if bad else 0)

print(json.dumps(data))
