#!/usr/bin/env python3
# Sub-grid Taylor deformation D vs time from the color function's mass-weighted second moments
# (D=(a-b)/(a+b), a>=b the principal semi-axes). Grid-cell-extent measures of D quantize too coarsely
# for these mild transient deformations; the moment estimate resolves the coverage dependence. Argv[1]
# = case dir. Compares against the interfacial-tension sigma = sigma0 + sigma_dGamma*surf_val printed
# from simulation.inp, so the coverage->deformation trend is legible.
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
dx = (inp("x_domain%end") - x0) / nx
xcc = x0 + (np.arange(nx) + 0.5) * dx
X = np.tile(xcc, ny).reshape(ny, nx)
Y = np.repeat(xcc, nx).reshape(ny, nx)
fs = sorted(
    [f for f in glob.glob(C + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
    key=lambda f: int(os.path.basename(f)[7:-4]),
)
ss = os.path.getsize(fs[0]) // 8 // (nx * ny)

print(f"{'t':>6} {'deformD':>9}")
for f in fs:
    c = np.fromfile(f, "<f8").reshape(ss, ny, nx)[ss - 2]
    w = np.clip(c, 0, 1)
    w = w / w.sum()
    mx = (X * w).sum()
    my = (Y * w).sum()
    Ixx = ((X - mx) ** 2 * w).sum()
    Iyy = ((Y - my) ** 2 * w).sum()
    Ixy = ((X - mx) * (Y - my) * w).sum()
    tr, det = Ixx + Iyy, Ixx * Iyy - Ixy**2
    l1 = tr / 2 + np.sqrt(tr * tr / 4 - det)
    l2 = tr / 2 - np.sqrt(tr * tr / 4 - det)
    a, b = np.sqrt(l1), np.sqrt(l2)
    print(f"{int(os.path.basename(f)[7:-4]) * dt:6.3f} {(a - b) / (a + b):9.5f}")
