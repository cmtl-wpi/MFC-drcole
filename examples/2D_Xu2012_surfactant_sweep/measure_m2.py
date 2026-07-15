#!/usr/bin/env python3
# Single-run M2 diagnostic: quasi-steady Taylor deformation D=(L-B)/(L+B), inclination theta, surfactant
# mass conservation, and a surfactant NON-UNIFORMITY metric (P90/median of Gamma on the interface band --
# high = concentrated at tips, low = spread out). Averages the last few saved frames for the quasi-steady
# value. Prints one JSON object. Argv[1] = case dir. Gamma = surf/|grad c| mirrors the solver.
import glob
import json
import os
import re
import sys

import numpy as np

C = sys.argv[1]


def inp(n):
    m = re.search(rf"^\s*{re.escape(n)}\s*=\s*([^\s,]+)", open(os.path.join(C, "simulation.inp")).read(), re.M)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


nx = int(inp("m")) + 1
ny = int(inp("n")) + 1
dt = inp("dt")
x0, y0 = inp("x_domain%beg"), inp("y_domain%beg")
dx = (inp("x_domain%end") - x0) / nx
dy = (inp("y_domain%end") - y0) / ny
X = np.tile(x0 + (np.arange(nx) + 0.5) * dx, ny).reshape(ny, nx)
Y = np.repeat(y0 + (np.arange(ny) + 0.5) * dy, nx).reshape(ny, nx)
fs = sorted(
    [f for f in glob.glob(C + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
    key=lambda f: int(os.path.basename(f)[7:-4]),
)
ss = os.path.getsize(fs[0]) // 8 // (nx * ny)
m0 = np.fromfile(fs[0], "<f8").reshape(ss, ny, nx)[ss - 1].sum()


def frame(f):
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    c, surf = a[ss - 2], a[ss - 1]
    w = np.clip(c, 0, 1)
    wn = w / w.sum()
    mx, my = (X * wn).sum(), (Y * wn).sum()
    Ixx = ((X - mx) ** 2 * wn).sum()
    Iyy = ((Y - my) ** 2 * wn).sum()
    Ixy = ((X - mx) * (Y - my) * wn).sum()
    tr = Ixx + Iyy
    disc = np.sqrt(max(tr * tr / 4 - (Ixx * Iyy - Ixy**2), 0.0))
    a_, b_ = np.sqrt(tr / 2 + disc), np.sqrt(max(tr / 2 - disc, 0.0))
    D = (a_ - b_) / (a_ + b_) if (a_ + b_) > 0 else 0.0
    theta = 0.5 * np.degrees(np.arctan2(2 * Ixy, Ixx - Iyy))
    # surfactant non-uniformity on the interface band: Gamma = surf/|grad c|
    gx = np.gradient(c, dx, axis=1)
    gy = np.gradient(c, dy, axis=0)
    normc = np.sqrt(gx**2 + gy**2)
    core = normc > 0.1 * normc.max()  # interface core, where Gamma = surf/|grad c| is physical
    gam = surf[core] / normc[core]
    gam = gam[np.isfinite(gam) & (gam > 0)]
    nonunif = float(np.percentile(gam, 90) / np.median(gam)) if gam.size else float("nan")
    return D, theta, surf.sum() / m0, nonunif


# Select the latest quasi-steady window where surfactant mass is still conserved (<2% drift): D has
# plateaued by then and this excludes the late-time tip-instability onset (see README). Robust to the
# onset time varying across sweep points.
allf = [(f, frame(f)) for f in fs]
valid = [(f, v) for f, v in allf if v[2] < 1.02]
use = valid[-4:] if len(valid) >= 4 else (valid or allf[-4:])
vals = np.array([v for _, v in use])
D, theta, mass, nonunif = vals.mean(axis=0)
t_meas = int(os.path.basename(use[-1][0])[7:-4]) * dt
print(json.dumps({"t_meas": round(float(t_meas), 3), "D": round(float(D), 5), "theta_deg": round(float(theta), 3), "mass_ratio": round(float(mass), 6), "surf_nonunif": round(float(nonunif), 4)}))
