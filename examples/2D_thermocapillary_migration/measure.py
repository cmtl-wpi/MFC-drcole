#!/usr/bin/env python3
"""Measure a thermocapillary drop's migration velocity and compare it to the Young-Goldstein-Block
(YGB) speed, for both geometry modes of case.py.

The drop is tracked by its color function c (1 inside, 0 outside). The reported velocity depends
on the boundary conditions, read from simulation.inp:
  * OPEN box (bc -3, drop centered): the YGB reference is an infinite quiescent domain, but the
    finite open box develops a small background return flow. We therefore report the
    DRIFT-CORRECTED speed U = u_drop - u_far, where u_drop is the color-weighted mean x-velocity
    and u_far the mean x-velocity far from the drop (|x| > 0.75*Lx/2). The correction is ~1-2% of
    v_YGB; it is not perfectly uniform (the drop's own flow field reaches the band), so treat the
    last percent of the ratio as uncertain.
  * SLIP-WALL box (bc -2, Samareh's geometry): the closed box defines the rest frame, so the
    lab-frame u_drop is the right measure (Samareh's own convention). No drift correction.

The migration / temperature-gradient axis is x. Each conserved-variable field is reshaped to
(-1, nx): x is MFC's fastest-varying index, so collapsing y (and z, in 3D) into rows makes this
script dimension-agnostic -- 2D and 3D are handled identically.

IMPORTANT: the reported "U_window" is the MEAN OVER THE FINAL VISCOUS TIME tau of the run -- a
windowed average, NOT necessarily a terminal velocity. The JSON also reports the final
instantaneous ratio and the slope of the ratio over the window (per tau): a run is only plateaued
if that slope is small. Default 3*tau runs are still rising at ~ +0.05-0.11 v_YGB/tau; see the
6*tau runs (SAMAREH_TAU=6) and README "Long-time behaviour" for true plateaus.

Run-dependent constants come from simulation.inp / pre_process.inp so this can't silently disagree
with the data. Only analysis choices (drop radius, imposed gradient, far-field cutoff) are set by
hand and MUST match case.py.

Usage:  python3 measure.py [case_dir]
Writes: <case_dir>/viz/migration_velocity.png  and prints a JSON summary line (tag: RESULT_JSON).
"""

import glob
import json
import os
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

case_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))


def read_namelist(path):
    """Parse a Fortran namelist file's plain "name = value" lines into a dict (lowercase keys)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                name, value = line.split("=", 1)
                out[name.strip().lower()] = value.strip().rstrip(",")
    return out


params = read_namelist(os.path.join(case_dir, "simulation.inp"))
patches = read_namelist(os.path.join(case_dir, "pre_process.inp"))


def param(name, src=params):
    """Look up a parameter by name and return it as a float."""
    return float(src[name.lower()])


nx = int(param("m")) + 1  # interior cells along the gradient axis x
dt = param("dt")
Lx = param("x_domain%end") - param("x_domain%beg")  # gradient-axis extent (7.5*D)
mu = 1.0 / param("fluid_pp(1)%re(1)")  # MFC takes Re = 1/mu
dsigma_dT = param("sigma_dtdt")  # sigma(T) slope
dim = 3 if int(param("p")) > 0 else 2
wall = int(param("bc_x%beg")) == -2  # slip-wall (Samareh geometry) vs open box
x_drop0 = param("patch_icpp(2)%x_centroid", patches)  # initial drop position

# Not in the namelists -- analysis choices. MUST match case.py.
r = 0.5  # droplet radius (D = 1)
gradT = 2.0 / 15.0  # imposed |dT/dx| (paper text rounds to 0.13; 1/(7.5D) = 2/15)
far_x = 0.75 * (Lx / 2.0)  # a cell is "far field" once |x| exceeds this (open mode)

tau = r**2 / mu  # viscous time rho*r^2/mu (rho(drop) = 1 by construction in case.py)
v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu  # YGB terminal speed (mu* = k* = 1)

# Conserved-variable layout (model_eqns=3, num_fluids=2, surface_tension):
#   indices 0,1 = partial densities, index 2 = x-momentum, last index = color c.
restart_dir = os.path.join(case_dir, "restart_data")

# Cell-center x positions from the boundary file (last nx+1 boundaries are the interior).
xb = np.fromfile(os.path.join(restart_dir, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
x = 0.5 * (xb[:-1] + xb[1:])
is_far = np.abs(x) > far_x

# Available snapshots: lustre_<step>.dat (skip the lustre_x_cb / y_cb / z_cb coordinate files).
steps = []
for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")):
    match = re.search(r"lustre_(\d+)\.dat$", f)
    if match:
        steps.append(int(match.group(1)))
steps.sort()
if not steps:
    sys.exit(f"no lustre_<step>.dat snapshots in {restart_dir!r} -- run the case first")

# Each snapshot is nvars conserved fields of `cells` values; infer nvars from the file size.
cells = (int(param("m")) + 1) * (int(param("n")) + 1) * (int(param("p")) + 1)
total = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size
nvars = total // cells
assert nvars * cells == total, f"file size {total} not divisible by cell count {cells}"


def field(snapshot, i):
    """Variable i of a snapshot, reshaped to (rows, nx) with x as the fast (column) axis.

    Collapsing y (and z) into rows makes the x-centroid / x-velocity / far-field math identical in
    2D and 3D, since x is MFC's innermost (contiguous) index.
    """
    return snapshot[i * cells : (i + 1) * cells].reshape(-1, nx)


# Walk the snapshots, building time series of the centroid position and the velocities.
times, x_centroid, u_lab, u_far = [], [], [], []
for s in steps:
    snap = np.fromfile(os.path.join(restart_dir, f"lustre_{s}.dat"), np.float64)
    vx = field(snap, 2) / (field(snap, 0) + field(snap, 1))  # x-velocity (rho = sum of partials)
    c = np.clip(field(snap, nvars - 1), 0.0, None)  # color function (last variable)
    times.append(s * dt)
    x_centroid.append((c * x[None, :]).sum() / c.sum())  # color-weighted x-centroid
    u_lab.append((c * vx).sum() / c.sum())  # color-weighted drop velocity
    u_far.append(vx[:, is_far].mean())  # far-field x-velocity (open-box drift)
times, x_centroid, u_lab, u_far = map(np.array, (times, x_centroid, u_lab, u_far))
U = u_lab if wall else u_lab - u_far  # lab frame in the closed box, drift-corrected in the open one

# Sanity: the first snapshot's centroid must sit at the patch position (catches a wrong reshape).
assert abs(x_centroid[0] - x_drop0) < 0.1 * r, f"drop not at patch centroid at t=0 ({x_centroid[0]:.3g} vs {x_drop0:.3g}) -- check layout"

# Mean and trend of U/v_YGB over the FINAL viscous time of the run (a window, not a terminal value).
in_window = times >= times[-1] - tau
U_window = U[in_window].mean()
ratio_window = U_window / v_YGB
ratio_final = U[-1] / v_YGB
slope_per_tau = np.polyfit(times[in_window] / tau, U[in_window] / v_YGB, 1)[0]

# Print the time series and the headline numbers.
mode = "slip-wall box (lab-frame U)" if wall else "open box (drift-corrected U)"
print(f"dim={dim}D  nx={nx}  cells={cells}  nvars={nvars}  mode: {mode}")
print(f"tau={tau:.3f}  v_YGB={v_YGB:.6f}  run length={times[-1] / tau:.1f} tau\n")
print(f"{'step':>7} {'t':>8} {'t/tau':>6} {'x_drop':>9} {'u_lab':>9} {'u_far':>9} {'U':>9} {'U/vYGB':>8}")
for i, t in enumerate(times):
    print(f"{steps[i]:>7} {t:>8.4f} {t / tau:>6.2f} {x_centroid[i]:>9.5f} {u_lab[i]:>9.5f} {u_far[i]:>9.5f} {U[i]:>9.5f} {U[i] / v_YGB:>+8.2f}")
print(f"\nmigration toward +x (hot side): {x_centroid[-1] > x_centroid[0]}")
print(f"U/v_YGB over the last tau: mean = {ratio_window:+.3f}, final = {ratio_final:+.3f}, trend = {slope_per_tau:+.3f}/tau")
if abs(slope_per_tau) > 0.02:
    print(f"  NOTE: still rising at {slope_per_tau:+.3f} v_YGB/tau -- this window mean is NOT a terminal velocity")

# Figure: U(t) against the v_YGB ceiling, with the averaging window shaded.
fig, ax = plt.subplots(figsize=(7.0, 4.5))
ax.plot(times / tau, U / v_YGB, "o-", color="C0", ms=3.5, label="MFC " + ("(lab frame)" if wall else "(drift-corrected)"))
ax.axhline(1.0, ls=":", color="0.45", lw=1.3)
ax.text(0.05, 1.0, r" $v_{\mathrm{YGB}}$ (zero-Ma Stokes, sphere)", va="bottom", ha="left", color="0.4", fontsize=9)
ax.axvspan(times[-1] / tau - 1.0, times[-1] / tau, color="0.92", zorder=0)
ax.annotate(
    rf"window mean $= {ratio_window:.2f}\,v_{{\mathrm{{YGB}}}}$" + (rf" (rising {slope_per_tau:+.2f}/$\tau$)" if abs(slope_per_tau) > 0.02 else " (plateau)"),
    xy=(times[-1] / tau, U_window / v_YGB),
    xytext=(0.35 * times[-1] / tau, 0.35),
    color="C0",
    fontsize=10,
    arrowprops=dict(arrowstyle="->", color="C0", lw=1.1),
)
ax.set_xlabel(r"$t/\tau$  ($\tau = \rho r^2/\mu$, viscous time)")
ax.set_ylabel(r"$U / v_{\mathrm{YGB}}$")
ax.set_xlim(left=0.0)
ax.set_ylim(0.0, max(1.1, 1.05 * (U / v_YGB).max()))
ax.set_title(f"{dim}D thermocapillary migration ({'Samareh slip-wall box' if wall else 'open box'})")
ax.grid(alpha=0.3)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
viz_dir = os.path.join(case_dir, "viz")
os.makedirs(viz_dir, exist_ok=True)
out_png = os.path.join(viz_dir, "migration_velocity.png")
fig.savefig(out_png, dpi=150)
print(f"saved figure -> {out_png}")

# Machine-readable summary for the validation driver.
summary = {
    "dim": dim,
    "nx_width": int(param("n")) + 1,
    "cells_per_D": (int(param("n")) + 1) / 5.0,  # box width = 5D
    "cells": cells,
    "mu": mu,
    "dsigma_dT": dsigma_dT,
    "gradT": gradT,
    "r": r,
    "v_YGB": v_YGB,
    "wall": wall,
    "drift_corrected": not wall,
    "U_window": float(U_window),
    "ratio_window": float(ratio_window),
    "ratio_final": float(ratio_final),
    "slope_per_tau": float(slope_per_tau),
    "u_lab_final": float(u_lab[-1]),
    "u_far_final": float(u_far[-1]),
    "migrates_hot": bool(x_centroid[-1] > x_centroid[0]),
    "t_end": float(times[-1]),
    "t_end_tau": float(times[-1] / tau),
}
print("RESULT_JSON " + json.dumps(summary))
