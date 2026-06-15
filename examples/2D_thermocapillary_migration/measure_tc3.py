#!/usr/bin/env python3
"""Measure the TC3 (large-Marangoni) drop migration of case_tc3.py and present it the way Samareh's
Figs 8/13 do: RISE VELOCITY (mm/s) vs DISTANCE FROM THE COLD WALL (mm).

The drop (Fluorinert, fluid 2) rises along +y through silicon oil; it is tracked by its color
function c (1 inside, 0 outside). The reported quantities are dimensional SI (the case is in SI):
  * distance from cold wall  = y_centroid - y_domain%beg            (m -> mm)
  * rise velocity            = color-weighted lab-frame v_y         (m/s -> mm/s)
The closed isothermal-wall cell is the rest frame, so the lab-frame velocity is the right measure.

The experiment's signature (Fig 8, linear initial T) is NON-MONOTONIC -- an initial overshoot, a dip,
then re-acceleration -- driven by the temperature-dependent viscosity mu(T) of the silicon oil (drag
falls as the drop rises into warmer, less-viscous oil). We do not have the digitized experimental
points, so this reports the MFC curve and the peak; the headline comparison is qualitative + a
converged production run.

Run-dependent constants come from simulation.inp so this can't silently disagree with the data.

Usage:  python3 measure_tc3.py [case_dir]
Writes: <case_dir>/viz/tc3_rise_velocity.png  and prints a JSON summary (tag: RESULT_JSON).
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
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


P = read_namelist(os.path.join(case_dir, "simulation.inp"))
f = lambda k: float(P[k.lower()])  # noqa: E731
nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
dt = f("dt")
y_cold = f("y_domain%beg")  # cold wall position (m)
cells = nx * ny * nz
restart = os.path.join(case_dir, "restart_data")

yb = np.fromfile(os.path.join(restart, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
y = 0.5 * (yb[:-1] + yb[1:])  # cell-center y (m)

steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(restart, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
if not steps:
    sys.exit(f"no snapshots in {restart!r} -- run case_tc3.py first")
nvars = np.fromfile(os.path.join(restart, f"lustre_{steps[0]}.dat"), np.float64).size // cells
c_idx = nvars - 2  # thermal_scalar appends T_s after the color function, so color is second-to-last

t_ms, dist_mm, vrise_mms = [], [], []
yb3 = y[None, :, None]
for s in steps:
    snap = np.fromfile(os.path.join(restart, f"lustre_{s}.dat"), np.float64)
    fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
    vy = fld(3) / (fld(0) + fld(1))
    c = np.clip(fld(c_idx), 0.0, 1.0)
    csum = c.sum()
    t_ms.append(s * dt)
    dist_mm.append(((c * yb3).sum() / csum - y_cold) * 1e3)  # distance from cold wall (mm)
    vrise_mms.append(((c * vy).sum() / csum) * 1e3)  # rise velocity (mm/s)
t_ms, dist_mm, vrise_mms = map(np.array, (t_ms, dist_mm, vrise_mms))

peak = float(vrise_mms.max()) if len(vrise_mms) else 0.0
print(f"nx={nx} ny={ny} nz={nz}  cells={cells}  nvars={nvars}  snapshots={len(steps)}")
print(f"run length = {t_ms[-1] * 1e3:.2f} ms   drop rose {dist_mm[0]:.1f} -> {dist_mm[-1]:.1f} mm from cold wall")
print(f"peak rise velocity = {peak:.3f} mm/s   [experiment Fig 8 peak ~ 2-3 mm/s]")

fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(dist_mm, vrise_mms, "o-", color="C0", ms=3, lw=1.0, label="MFC")
ax.set_xlabel("Distance from cold wall (mm)")
ax.set_ylabel("Rise velocity (mm/s)")
ax.set_title("TC3: large-Ma migration (Samareh Fig 8/13) -- MFC")
ax.grid(alpha=0.3)
ax.legend(loc="best", fontsize=9)
fig.tight_layout()
viz = os.path.join(case_dir, "viz")
os.makedirs(viz, exist_ok=True)
out = os.path.join(viz, "tc3_rise_velocity.png")
fig.savefig(out, dpi=150)
print(f"saved figure -> {out}")

summary = {
    "nx": nx,
    "ny": ny,
    "nz": nz,
    "cells": cells,
    "snapshots": len(steps),
    "t_end_ms": float(t_ms[-1] * 1e3),
    "dist_start_mm": float(dist_mm[0]),
    "dist_end_mm": float(dist_mm[-1]),
    "peak_rise_velocity_mms": peak,
    "experiment_peak_mms": "2-3 (Fig 8)",
}
print("RESULT_JSON " + json.dumps(summary))
