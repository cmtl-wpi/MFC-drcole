#!/usr/bin/env python3
# Recovering-sphere validation for axisymmetric surface tension.
#
# Purpose: validate the cyl_coord + surface_tension capillary geometric source
# (the azimuthal/hoop capillary stress added in m_surface_tension.fpp). Two
# modes, each with a sharp analytic target:
#
#   --mode oscillate (default): a volume-preserving prolate spheroid (axial
#     semi-axis 1.08*R_eff) relaxes toward a sphere, oscillating in the l=2
#     Rayleigh-Lamb mode. Frequency target:
#         omega^2 = 24*sigma / ((3*rho_l + 2*rho_g) * R^3)   (sphere, l=2)
#     The DIAGNOSTIC power: a solver that is secretly planar (missing the
#     azimuthal curvature) oscillates instead at the 2D-drop frequency
#         omega^2 ~ 6*sigma / ((rho_l + rho_g) * R^3)
#     i.e. a ~15% longer period - unambiguously measurable.
#   --mode static: a resting sphere must hold the Laplace jump 2*sigma/R
#     (time-averaged). The unfixed code holds ~sigma/R and/or erupts (see the
#     film-zoom GATE-1 record: instant 8.8 m/s spurious flow, dt -> 0 abort).
#
# EOS family: artificial liquid sound speed c_l = 100 m/s (as in the
# 2D_droplet_coalescence case-s runs) so the runs take minutes, not days.
# R = 50 um: T_sphere = 2*pi/omega ~ 133 us; two periods ~ 266 us.

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--mfc", type=json.loads, default="{}")
parser.add_argument("--mode", type=str, default="oscillate", choices=["oscillate", "static"])
parser.add_argument("--cfl", type=float, default=0.1)
parser.add_argument("--t-stop", type=float, default=None, dest="t_stop")
parser.add_argument("--n_start", type=int, default=0)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

# Fluids (tetradecane / N2, artificial-c_l family)
rho_l = 763.0
rho_g = rho_l / 666.0
sigma = 0.0266
gamma_l = 3.7
gamma_g = 1.4
c_l = 100.0

p_gas = c_l**2 * rho_g / gamma_g
pi_inf_l = rho_l * c_l**2 / gamma_l - p_gas

# Drop
R = 50e-6
p_liq = p_gas + 2.0 * sigma / R  # spherical Laplace jump

# l=2 Rayleigh-Lamb target (sphere with outer fluid)
omega2 = 24.0 * sigma / ((3.0 * rho_l + 2.0 * rho_g) * R**3)
T_sphere = 2.0 * math.pi / math.sqrt(omega2)
# planar-drop (cylinder) frequency, for the discrimination bound
omega2_cyl = 6.0 * sigma / ((rho_l + rho_g) * R**3)
T_cyl = 2.0 * math.pi / math.sqrt(omega2_cyl)

if args.mode == "oscillate":
    stretch_ax = 1.08  # axial semi-axis a0 = 1.08 R; radial b0 = R/sqrt(1.08)
    t_stop = args.t_stop if args.t_stop is not None else 3.0e-4  # ~2.25 periods
    t_save = 2.0e-6
else:
    stretch_ax = 1.0
    t_stop = args.t_stop if args.t_stop is not None else 1.0e-4
    t_save = 5.0e-6
a0 = stretch_ax * R
b0 = R / math.sqrt(stretch_ax)

# Grid: UNIFORM axisymmetric (no stretching - keep the validation clean).
# x in [-4R, 4R] axial, y in [0, 4R] radial; dx = dy_interior ~ R/40 = 1.25 um.
Lx = 4.0 * R
m = 319   # dx = 8R/320 = 1.25 um
n = 160   # axisym offset grid: dy_interior = 2*(4R)/(2n+1) ~ 1.246 um

mu_l = 2.1823e-3  # same liquid viscosity as the case-s family (Oh ~ 0.05)
mu_g = mu_l / 119.0

data = {
    "run_time_info": "T",
    "x_domain%beg": -Lx,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Lx,
    "m": m,
    "n": n,
    "p": 0,
    "cyl_coord": "T",
    "bc_x%beg": -6,
    "bc_x%end": -6,
    "bc_y%beg": -2,  # axis
    "bc_y%end": -6,
    "n_start": args.n_start,
    "cfl_adap_dt": "T",
    "cfl_target": args.cfl,
    "t_stop": t_stop,
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
    "num_patches": 2,
    "num_fluids": 2,
    "surface_tension": "T",
    "sigma": sigma,
    "format": 1,
    "precision": 2,
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "pres_wrt": "T",
    "parallel_io": "T",
    "fd_order": 2,
    "probe_wrt": "T",
    "num_probes": 2,
    "probe(1)%x": 0.0,      # drop centre: internal pressure (Laplace check)
    "probe(1)%y": 1.0e-6,
    "probe(2)%x": 3.0 * R,  # far gas reference
    "probe(2)%y": 3.0 * R,
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    # patch 1: N2 background
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 4.0 * Lx,
    "patch_icpp(1)%length_y": 4.0 * Lx,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": 1e-9 * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - 1e-9) * rho_g,
    "patch_icpp(1)%alpha(1)": 1e-9,
    "patch_icpp(1)%alpha(2)": 1.0 - 1e-9,
    "patch_icpp(1)%cf_val": 0,
    # patch 2: the drop - ellipse (geometry 5), axisym ellipse = spheroid;
    # degenerate to a circle in static mode via equal radii
    "patch_icpp(2)%geometry": 5,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radii(1)": a0,
    "patch_icpp(2)%radii(2)": b0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.99,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - 1e-9) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": 1e-9 * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - 1e-9,
    "patch_icpp(2)%alpha(2)": 1e-9,
    "patch_icpp(2)%cf_val": 1,
}

if args.check:
    dx = 2.0 * Lx / (m + 1)
    dy = 2.0 * Lx / (2 * n + 1)
    dt = args.cfl * min(dx, dy) / (c_l + 1.0)
    print(f"mode={args.mode}  R={R*1e6:.0f}um  a0={a0*1e6:.2f} b0={b0*1e6:.2f}um")
    print(f"T_sphere(l=2)={T_sphere*1e6:.2f}us  T_cylinder={T_cyl*1e6:.2f}us  ratio={T_cyl/T_sphere:.4f}")
    print(f"2sigma/R={2*sigma/R:.1f} Pa  p_gas={p_gas:.1f}  p_liq={p_liq:.1f}")
    print(f"grid {m+1}x{n+1}={(m+1)*(n+1)} cells  dx={dx*1e6:.3f} dy_int={dy*2*1e6:.3f}um")
    print(f"dt~{dt:.3e}  steps~{t_stop/dt:.2e}  t_stop={t_stop*1e6:.0f}us  frames={int(t_stop/t_save)}")
    checks = {
        "volume preserved": abs(a0 * b0**2 - R**3) < 1e-25,
        "alphas sum p1": abs(data["patch_icpp(1)%alpha(1)"] + data["patch_icpp(1)%alpha(2)"] - 1) < 1e-15,
        "alphas sum p2": abs(data["patch_icpp(2)%alpha(1)"] + data["patch_icpp(2)%alpha(2)"] - 1) < 1e-15,
        "spherical Laplace IC": abs((p_liq - p_gas) - 2 * sigma / R) < 1e-12,
        "axis bc": data["bc_y%beg"] == -2,
        "axisym": data["cyl_coord"] == "T",
    }
    bad = [k for k, v in checks.items() if not v]
    for k, v in checks.items():
        print(f"  [{'PASS' if v else 'FAIL'}] {k}")
    sys.exit(1 if bad else 0)

print(json.dumps(data))
