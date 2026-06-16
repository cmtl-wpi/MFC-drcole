#!/usr/bin/env python3
"""Measure a 3D thermocapillary sphere's RISE velocity and compare it to the Young-Goldstein-Block
(YGB) speed (Samareh 2014, Sec. 4.1.1 / Fig 6, the fully-3D sphere -> v_t/v_YGB ~ 0.95).

Dimension-agnostic: each field is reshaped to (nz, ny, nx), so the color-weighted centroid / rise
velocity / far-field math is identical to the 2D example; only the reference ratio differs (0.95 vs 0.80).

The drop rises along +y (the 7.5D gradient axis) and is tracked by its color function c
(1 inside, 0 outside). The reported velocity depends on the boundary conditions, read from
simulation.inp:
  * SLIP-WALL box (bc -2, Samareh's Fig 5/6 geometry, the default): the closed box defines the
    rest frame, so the lab-frame rise velocity u_drop is the right measure (Samareh's own
    convention). No drift correction. Compare against Samareh ~0.80 (2D) / ~0.95 (3D).
  * OPEN box (bc -3, drop centered): the YGB reference is an infinite quiescent domain, but the
    finite open box develops a small background return flow. We therefore report the
    DRIFT-CORRECTED speed U = u_drop - u_far, where u_far is the mean y-velocity far from the drop
    (|y| > 0.75*Ly/2). The 2D anchor is then the unbounded-cylinder analytic 15/16*v_YGB.

Each conserved-variable field is reshaped to (nz, ny, nx) with x MFC's fastest-varying index, so
the y-resolved centroid / rise velocity / far-field math is identical in 2D (nz=1) and 3D.

IMPORTANT: the reported "U_window" is the MEAN OVER THE FINAL VISCOUS TIME tau of the run -- a
windowed average, NOT necessarily a terminal velocity. The JSON also reports the final
instantaneous ratio and the slope of the ratio over the window (per t_r): a run is only plateaued
if that slope is small.

Run-dependent constants come from simulation.inp / pre_process.inp so this can't silently disagree
with the data. Only analysis choices (drop radius, imposed gradient, far-field cutoff) are set by
hand and MUST match case.py.

Usage:  python3 measure.py [case_dir]
Writes: <case_dir>/viz/rise_velocity.png  and prints a JSON summary line (tag: RESULT_JSON).
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
# thermal_scalar appends an independent temperature scalar (eqn_idx%T_s) AFTER the color function,
# so in that mode the color function is the second-to-last conserved variable, not the last.
ts_mode = str(params.get("thermal_scalar", "F")).upper().strip(". ").startswith("T")
# Bulk Fourier conduction tames the frozen-T runaway; when it's on (this 3D case by construction)
# the "frozen-T drift" caveat -- a 2D no-conduction artifact -- does not apply.
cond_mode = str(params.get("thermal_conduction", "F")).upper().strip(". ").startswith("T")


def param(name, src=params):
    """Look up a parameter by name and return it as a float."""
    return float(src[name.lower()])


nx = int(param("m")) + 1  # short-axis cells (x)
ny = int(param("n")) + 1  # rise-axis cells (y)
nz = int(param("p")) + 1  # short-axis cells (z); 1 in 2D
dt = param("dt")
Ly = param("y_domain%end") - param("y_domain%beg")  # rise-axis extent (7.5*D)
mu = 1.0 / param("fluid_pp(1)%re(1)")  # MFC takes Re = 1/mu
dsigma_dT = param("sigma_dtdt")  # sigma(T) slope
dim = 3 if int(param("p")) > 0 else 2
wall = int(param("bc_y%beg")) == -2  # slip-wall (Samareh geometry) vs open box
y_drop0 = param("patch_icpp(2)%y_centroid", patches)  # initial drop position

# Not in the namelists -- analysis choices. MUST match case.py.
r = 0.5  # droplet radius (D = 1)
gradT = 2.0 / 15.0  # imposed |dT/dy| (paper text rounds to 0.13; 1/(7.5D) = 2/15)
far_y = 0.75 * (Ly / 2.0)  # a cell is "far field" once |y| exceeds this (open mode)

t_r = mu / abs(dsigma_dT * gradT)  # capillary-thermal time mu/|sigma_T*gradT| (Samareh's time scale)
v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu  # YGB terminal speed (mu* = k* = 1)

# Samareh's converged ratio: Fig 6 (3D sphere) ~ 0.95; Fig 5 (2D cylinder) ~ 0.80. This example is 3D.
samareh_ratio = 0.95 if dim == 3 else 0.80

# Conserved-variable layout (model_eqns=3, num_fluids=2, surface_tension):
#   indices 0,1 = partial densities, 2 = x-momentum, 3 = y-momentum, color c is last --
#   or second-to-last when thermal_scalar appends T_s after it.
restart_dir = os.path.join(case_dir, "restart_data")

# Cell-center y positions from the boundary file (last ny+1 boundaries are the interior).
yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
y = 0.5 * (yb[:-1] + yb[1:])  # shape (ny,)
is_far = np.abs(y) > far_y

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
cells = nx * ny * nz
total = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size
nvars = total // cells
assert nvars * cells == total, f"file size {total} not divisible by cell count {cells}"


def field(snapshot, i):
    """Variable i of a snapshot, reshaped to (nz, ny, nx) -- x is MFC's fast (innermost) index."""
    return snapshot[i * cells : (i + 1) * cells].reshape(nz, ny, nx)


# Walk the snapshots, building time series of the y-centroid and the rise velocities.
times, y_centroid, u_lab, u_far, rho_drop_t = [], [], [], [], []
yb3 = y[None, :, None]  # broadcast y over (nz, ny, nx)
for s in steps:
    snap = np.fromfile(os.path.join(restart_dir, f"lustre_{s}.dat"), np.float64)
    rho = field(snap, 0) + field(snap, 1)  # total density (sum of partial densities)
    vy = field(snap, 3) / rho  # rise (y) velocity
    c_idx = nvars - 2 if ts_mode else nvars - 1  # T_s is appended after the color function in ts_mode
    c = np.clip(field(snap, c_idx), 0.0, None)  # color function
    csum = c.sum()
    times.append(s * dt)
    y_centroid.append((c * yb3).sum() / csum)  # color-weighted y-centroid
    u_lab.append((c * vy).sum() / csum)  # color-weighted drop rise velocity
    u_far.append(vy[:, is_far, :].mean())  # far-field y-velocity (open-box drift)
    rho_drop_t.append((c * rho).sum() / csum)  # color-weighted drop density
times, y_centroid, u_lab, u_far = map(np.array, (times, y_centroid, u_lab, u_far))
U = u_lab if wall else u_lab - u_far  # lab frame in the closed box, drift-corrected in the open one

# Drop density and viscous time tau come from the data (t=0 drop), not an assumption.
rho_drop = float(rho_drop_t[0])
tau = rho_drop * r**2 / mu  # viscous time rho*r^2/mu (~0.5 for Samareh's rho=0.2)

# Sanity: the first snapshot's centroid must sit at the patch position (catches a wrong reshape).
assert abs(y_centroid[0] - y_drop0) < 0.1 * r, f"drop not at patch centroid at t=0 ({y_centroid[0]:.3g} vs {y_drop0:.3g}) -- check layout"

# QUASI-STEADY PLATEAU (the value to compare with Samareh). With no bulk conduction the frozen
# linear-T field is slowly advected, so the rise velocity does NOT stay flat like Samareh's
# invariant-T case: it ramps up, OVERSHOOTS, settles to a quasi-steady plateau, then slowly DRIFTS
# (upward) as T distorts. The faithful comparison to Samareh's plateau is that settled value -- the
# post-overshoot minimum of the smoothed curve, where the internal circulation is fully developed
# but T is not yet appreciably distorted -- NOT the drifted endpoint. We also report the endpoint
# and the late-time drift slope as diagnostics of how far the frozen-IC approximation has decayed.
ratio_t = U / v_YGB
dt_snap = times[1] - times[0] if len(times) > 1 else tau
smooth_w = max(1, int(round(tau / dt_snap)))  # rolling mean over ~ one viscous time
if smooth_w > 1:  # edge-corrected boxcar: divide by the window's actual coverage so the ends aren't
    k = np.ones(smooth_w)  # depressed by zero-padding (which would put a spurious minimum at t_end)
    sm = np.convolve(ratio_t, k, mode="same") / np.convolve(np.ones_like(ratio_t), k, mode="same")
else:
    sm = ratio_t
peak_i = int(np.argmax(np.where(times < 0.6 * t_r, sm, -np.inf)))  # overshoot peak (first 0.6 t_r)
plateau_i = peak_i + int(np.argmin(sm[peak_i:]))  # quasi-steady plateau = post-overshoot minimum
ratio_plateau = float(sm[plateau_i])
t_plateau_tr = float(times[plateau_i] / t_r)
overshoot = float(sm[peak_i])
ratio_final = float(ratio_t[-1])
# Late-time drift: slope of the ratio over the final viscous time (small => frozen-IC still holds).
in_tail = times >= times[-1] - tau
slope_per_tr = float(np.polyfit(times[in_tail] / t_r, ratio_t[in_tail], 1)[0]) if in_tail.sum() > 1 else 0.0

# The "frozen-T drift" caveat applies ONLY without bulk conduction. With conduction on, a positive
# late-time slope at t < ~t_r is just the unfinished ramp, not drift. A run shorter than ~1 t_r has
# not had time to overshoot and settle, so its "plateau" is really the current (still-ramping) value.
drift_is_frozenT = (not cond_mode) and abs(slope_per_tr) > 0.05
short_run = times[-1] < 1.0 * t_r

# Print the time series and the headline numbers.
mode = "slip-wall box (lab-frame U)" if wall else "open box (drift-corrected U)"
print(f"dim={dim}D  nx={nx} ny={ny} nz={nz}  cells={cells}  nvars={nvars}  mode: {mode}")
print(f"rho_drop={rho_drop:.4f}  tau={tau:.3f}  t_r={t_r:.3f}  v_YGB={v_YGB:.6f}  run length={times[-1] / t_r:.2f} t_r\n")
print(f"{'step':>7} {'t':>8} {'t/t_r':>6} {'y_drop':>9} {'u_lab':>9} {'u_far':>9} {'U':>9} {'U/vYGB':>8}")
for i, t in enumerate(times):
    print(f"{steps[i]:>7} {t:>8.4f} {t / t_r:>6.2f} {y_centroid[i]:>9.5f} {u_lab[i]:>9.5f} {u_far[i]:>9.5f} {U[i]:>9.5f} {U[i] / v_YGB:>+8.2f}")
print(f"\nrising (+y, toward hot top): {y_centroid[-1] > y_centroid[0]}")
print(f"quasi-steady plateau  v_t/v_YGB = {ratio_plateau:+.3f}  (at t/t_r = {t_plateau_tr:.2f})   [Samareh {dim}D ~ {samareh_ratio:.2f}]")
print(f"overshoot peak = {overshoot:+.3f}   endpoint = {ratio_final:+.3f}   late-time slope = {slope_per_tr:+.3f}/t_r")
if drift_is_frozenT:
    print(f"  NOTE: endpoint is frozen-T-drift contaminated ({slope_per_tr:+.3f}/t_r) -- the plateau, not the endpoint, is the Samareh comparison")
elif short_run:
    print(f"  NOTE: run is only {times[-1] / t_r:.2f} t_r -- the rise is still ramping, not yet a settled plateau (need >~1 t_r for the Samareh comparison)")

# Figure: U(t)/v_YGB vs t/t_r against the v_YGB ceiling and Samareh's converged ratio.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(times / t_r, ratio_t, "o-", color="C0", ms=3.0, lw=1.0, label="MFC " + ("(lab frame)" if wall else "(drift-corrected)"))
ax.axhline(1.0, ls=":", color="0.45", lw=1.3)
ax.text(0.02, 1.0, r" $v_{\mathrm{YGB}}$ (zero-Ma Stokes, sphere)", va="bottom", ha="left", color="0.4", fontsize=9)
ax.axhline(samareh_ratio, ls="--", color="C3", lw=1.3, label=rf"Samareh {dim}D $\approx$ {samareh_ratio:.2f}")
if short_run:
    ax.plot(times[-1] / t_r, ratio_final, "*", color="k", ms=14, zorder=5, label=rf"current value = {ratio_final:.2f} (still ramping, {times[-1] / t_r:.2f} $t_r$)")
else:
    ax.plot(t_plateau_tr, ratio_plateau, "*", color="k", ms=14, zorder=5, label=rf"quasi-steady plateau = {ratio_plateau:.2f}")
ax.annotate(
    "frozen-$T$ drift\n(no bulk conduction)" if drift_is_frozenT else "",
    xy=(times[-1] / t_r, ratio_final),
    xytext=(0.55 * times[-1] / t_r, min(1.08, ratio_final + 0.12)),
    color="0.4",
    fontsize=8,
    arrowprops=dict(arrowstyle="->", color="0.5", lw=1.0) if drift_is_frozenT else None,
)
ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
ax.set_ylabel(r"normalized rise velocity  $v / v_{\mathrm{YGB}}$")
ax.set_xlim(left=0.0)
ax.set_ylim(0.0, max(1.1, 1.05 * ratio_t.max()))
ax.set_title(f"{dim}D thermocapillary rise ({'Samareh slip-wall box' if wall else 'open box'})")
ax.grid(alpha=0.3)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
viz_dir = os.path.join(case_dir, "viz")
os.makedirs(viz_dir, exist_ok=True)
out_png = os.path.join(viz_dir, "rise_velocity.png")
fig.savefig(out_png, dpi=150)
print(f"saved figure -> {out_png}")

# Machine-readable summary for the validation driver.
summary = {
    "dim": dim,
    "nx_width": nx,
    "cells_per_D": nx / 5.0,  # box width = 5D
    "cells": cells,
    "mu": mu,
    "rho_drop": rho_drop,
    "dsigma_dT": dsigma_dT,
    "gradT": gradT,
    "r": r,
    "v_YGB": v_YGB,
    "tau": tau,
    "t_r": t_r,
    "wall": wall,
    "drift_corrected": not wall,
    "samareh_ratio": samareh_ratio,
    "ratio_plateau": ratio_plateau,
    "t_plateau_tr": t_plateau_tr,
    "overshoot": overshoot,
    "ratio_final": ratio_final,
    "slope_per_tr": slope_per_tr,
    "u_lab_final": float(u_lab[-1]),
    "u_far_final": float(u_far[-1]),
    "rises": bool(y_centroid[-1] > y_centroid[0]),
    "t_end": float(times[-1]),
    "t_end_tr": float(times[-1] / t_r),
}
print("RESULT_JSON " + json.dumps(summary))
