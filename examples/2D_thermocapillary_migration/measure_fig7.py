#!/usr/bin/env python3
"""Measure the finite-Marangoni thermocapillary migration of case_fig7.py and compare it to
Samareh-2014 Fig 7 (their Sec. 4.1.2, the Nas & Tryggvason Re=5/Ma=20/Ca=0.01666 test).

Samareh report the migration as U* = U/U_r versus t* = t/t_r, where
    U_r = |sigma_T gradT| r_0 / mu_b   (Marangoni velocity scale)
    t_r = mu_b / |sigma_T gradT|       (capillary-thermal time scale)
The published curve ramps from rest, OVERSHOOTS to U* ~ 0.13 near t* ~ 5, then relaxes toward a
terminal value (~0.10 by t* ~ 20); Samareh's fine grid matches Nas & Tryggvason to within 1.7%.

The drop is a REAL second fluid (all material properties 0.5x the bulk) tracked by its color
function c (1 inside, 0 outside). The box is a CLOSED slip-wall box, so the box frame is the rest
frame and the lab-frame color-weighted rise velocity U is the right measure (no drift correction).

Run-dependent constants (grid, dt, mu_b, sigma_T, the imposed gradient gradT) are read from
simulation.inp so this can't silently disagree with the case. gradT is recovered from the imposed
wall temperatures, (Twall_out - Twall_in)/Ly, NOT hard-coded.

Usage:  python3 measure_fig7.py [case_dir]
Writes: <case_dir>/viz/fig7_migration.png  and prints a JSON summary line (tag: RESULT_JSON).
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

# Samareh Fig 7 reference: overshoot peak of U* (Nas & Tryggvason test). Terminal ~0.10 is qualitative
# (we do not have the digitized N&T points), so the headline comparison is the peak and the 64-vs-128
# grid agreement (Samareh: 1.2% between grids, 1.7% vs Nas & Tryggvason on the fine grid).
NAS_TRYGGVASON_PEAK = 0.13


def read_namelist(path):
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
    return float(src[name.lower()])


nx = int(param("m")) + 1  # box-width cells (x)
ny = int(param("n")) + 1  # rise-axis cells (y)
nz = int(param("p")) + 1  # 1 in 2D
dt = param("dt")
Wx = param("x_domain%end") - param("x_domain%beg")  # box width (= 2D)
Ly = param("y_domain%end") - param("y_domain%beg")  # rise-axis extent (= 4D)
mu_b = 1.0 / param("fluid_pp(1)%re(1)")  # bulk viscosity (MFC takes Re = 1/mu)
sigma_T = param("sigma_dTdT")  # dsigma/dT
y_drop0 = param("patch_icpp(2)%y_centroid", patches)  # initial drop position

# Analysis choices (MUST match case_fig7.py).
D = 1.0  # droplet diameter sets the length scale
r0 = D / 2.0  # droplet radius = 0.5
# Imposed gradient: prefer the isothermal wall temperatures (the conduction case); fall back to the
# IC construction T = T_base + (1/Ly)*(y - y_bottom), i.e. gradT = 1/Ly, when there is no isothermal
# BC (the advection-only FIG7_COND=0 diagnostic). Both give the same value by construction.
if "bc_y%twall_out" in params and "bc_y%twall_in" in params:
    gradT = abs(param("bc_y%Twall_out") - param("bc_y%Twall_in")) / Ly
else:
    gradT = 1.0 / Ly

G = abs(sigma_T * gradT)  # |sigma_T gradT| (Marangoni stress magnitude)
U_r = G * r0 / mu_b  # Samareh's velocity scale
t_r = mu_b / G  # Samareh's time scale

restart_dir = os.path.join(case_dir, "restart_data")
yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
y = 0.5 * (yb[:-1] + yb[1:])

steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
if not steps:
    sys.exit(f"no lustre_<step>.dat snapshots in {restart_dir!r} -- run the case first")

cells = nx * ny * nz
nvars = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size // cells
# thermal_scalar appends the independent temperature scalar T_s after the color function, so the
# color function is the second-to-last conserved variable here.
c_idx = nvars - 2


def field(snapshot, i):
    return snapshot[i * cells : (i + 1) * cells].reshape(nz, ny, nx)


# Walk the snapshots: color-weighted drop centroid and lab-frame rise velocity.
yb3 = y[None, :, None]
times, y_centroid, u_lab = [], [], []
for s in steps:
    snap = np.fromfile(os.path.join(restart_dir, f"lustre_{s}.dat"), np.float64)
    rho = field(snap, 0) + field(snap, 1)  # total density (sum of partial densities)
    vy = field(snap, 3) / rho  # rise (y) velocity
    c = np.clip(field(snap, c_idx), 0.0, 1.0)  # color function (drop = 1)
    csum = c.sum()
    times.append(s * dt)
    y_centroid.append((c * yb3).sum() / csum)
    u_lab.append((c * vy).sum() / csum)
times = np.array(times)
y_centroid = np.array(y_centroid)
U = np.array(u_lab)
t_star = times / t_r
U_star = U / U_r

# Sanity: the first snapshot's centroid must sit at the patch position.
assert abs(y_centroid[0] - y_drop0) < 0.1 * r0, f"drop not at patch centroid at t=0 ({y_centroid[0]:.3g} vs {y_drop0:.3g})"

# Overshoot peak and terminal value (mean over the final t_r window).
peak_i = int(np.argmax(U_star))
peak = float(U_star[peak_i])
t_peak = float(t_star[peak_i])
in_tail = times >= times[-1] - t_r
terminal = float(U_star[in_tail].mean())

print(f"nx={nx} ny={ny} nz={nz}  cells={cells}  cells/D={nx / (Wx / D):.1f}  nvars={nvars}  color idx={c_idx}")
print(f"mu_b={mu_b:.4f}  sigma_T={sigma_T:.5f}  gradT={gradT:.5f}  U_r={U_r:.5f}  t_r={t_r:.4f}  run={t_star[-1]:.2f} t_r\n")
print(f"{'step':>8} {'t*':>7} {'y_drop':>9} {'U':>10} {'U*':>9}")
for i in range(len(times)):
    print(f"{steps[i]:>8} {t_star[i]:>7.2f} {y_centroid[i]:>9.5f} {U[i]:>10.6f} {U_star[i]:>9.4f}")
print(f"\nrising (+y, toward hot top): {y_centroid[-1] > y_centroid[0]}")
print(f"overshoot peak  U* = {peak:.4f}  at t* = {t_peak:.2f}   [Nas & Tryggvason / Samareh Fig 7 peak ~ {NAS_TRYGGVASON_PEAK:.2f}]")
print(f"terminal (final t_r window)  U* = {terminal:.4f}")

# Figure: U*(t*) against the Nas & Tryggvason peak band.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(t_star, U_star, "o-", color="C0", ms=3.0, lw=1.2, label=f"MFC ({nx}/box width, {nx / (Wx / D):.0f}/$D$)")
ax.axhline(NAS_TRYGGVASON_PEAK, ls="--", color="C3", lw=1.3, label=rf"Nas & Tryggvason peak $\approx$ {NAS_TRYGGVASON_PEAK:.2f}")
ax.plot(t_peak, peak, "*", color="k", ms=14, zorder=5, label=rf"MFC peak = {peak:.3f}")
ax.set_xlabel(r"$t^* = t / t_r$  ($t_r = \mu_b / |\sigma_T \nabla T|$)")
ax.set_ylabel(r"$U^* = U / U_r$  ($U_r = |\sigma_T \nabla T|\, r_0 / \mu_b$)")
ax.set_title("Fig 7: finite-Ma thermocapillary migration (Re=5, Ma=20, Ca=0.01666)")
ax.set_xlim(left=0.0)
ax.set_ylim(bottom=min(0.0, 1.05 * U_star.min()))
ax.grid(alpha=0.3)
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
viz_dir = os.path.join(case_dir, "viz")
os.makedirs(viz_dir, exist_ok=True)
out_png = os.path.join(viz_dir, "fig7_migration.png")
fig.savefig(out_png, dpi=150)
print(f"saved figure -> {out_png}")

summary = {
    "nx_width": nx,
    "cells_per_D": nx / (Wx / D),
    "cells": cells,
    "mu_b": mu_b,
    "sigma_T": sigma_T,
    "gradT": gradT,
    "U_r": U_r,
    "t_r": t_r,
    "peak": peak,
    "t_peak_tr": t_peak,
    "terminal": terminal,
    "nas_tryggvason_peak": NAS_TRYGGVASON_PEAK,
    "rises": bool(y_centroid[-1] > y_centroid[0]),
    "t_end_tr": float(t_star[-1]),
}
print("RESULT_JSON " + json.dumps(summary))
