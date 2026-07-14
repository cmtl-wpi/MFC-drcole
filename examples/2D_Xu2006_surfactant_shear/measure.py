#!/usr/bin/env python3
# Taylor deformation D=(L-B)/(L+B) and inclination angle theta (major axis vs x, in degrees) of the drop
# vs time, from the color function's mass-weighted second-moment (inertia) tensor -- the Xu 2006 shear
# diagnostics. Also interfacial surfactant mass (conservation). Argv[1] = case dir.
import glob
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
x0 = inp("x_domain%beg")
y0 = inp("y_domain%beg")
dx = (inp("x_domain%end") - x0) / nx
dy = (inp("y_domain%end") - y0) / ny
X = np.tile(x0 + (np.arange(nx) + 0.5) * dx, ny).reshape(ny, nx)
Y = np.repeat(y0 + (np.arange(ny) + 0.5) * dy, nx).reshape(ny, nx)
fs = sorted(
    [f for f in glob.glob(C + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
    key=lambda f: int(os.path.basename(f)[7:-4]),
)
ss = os.path.getsize(fs[0]) // 8 // (nx * ny)

print(f"{'t':>7} {'D':>8} {'theta_deg':>9} {'surfmass/m0':>11}")
m0 = None
for f in fs:
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    c, surf = a[ss - 2], a[ss - 1]
    w = np.clip(c, 0, 1)
    sm = w.sum()
    w = w / sm
    mx = (X * w).sum()
    my = (Y * w).sum()
    Ixx = ((X - mx) ** 2 * w).sum()
    Iyy = ((Y - my) ** 2 * w).sum()
    Ixy = ((X - mx) * (Y - my) * w).sum()
    tr, det = Ixx + Iyy, Ixx * Iyy - Ixy**2
    l1 = tr / 2 + np.sqrt(max(tr * tr / 4 - det, 0.0))
    l2 = tr / 2 - np.sqrt(max(tr * tr / 4 - det, 0.0))
    a_, b_ = np.sqrt(l1), np.sqrt(max(l2, 0.0))
    D = (a_ - b_) / (a_ + b_) if (a_ + b_) > 0 else 0.0
    # major-axis eigenvector -> inclination from x-axis
    theta = 0.5 * np.degrees(np.arctan2(2 * Ixy, Ixx - Iyy))
    if m0 is None:
        m0 = surf.sum()
    print(f"{int(os.path.basename(f)[7:-4]) * dt:7.3f} {D:8.5f} {theta:9.3f} {surf.sum() / m0:11.6f}")
