#!/usr/bin/env python3
"""
Measure the l=1 surfactant surface-diffusion decay rate on the sphere and compare to the exact
Laplace-Beltrami rate 2 D_s/R^2. Prints "R/dx rate exact" (one line) for the convergence sweep.

    python3 measure.py            # print the rate for the run currently in this directory
Reads grid/dt/D_s from simulation.inp; R is the example constant below.
"""

import glob
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = 0.5  # sphere radius (fixed in case.py)


def inp(name, path):
    m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([^\s,]+)", open(path).read(), re.M)
    if not m:
        raise KeyError(name)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


sim = os.path.join(HERE, "simulation.inp")
nx = int(inp("m", sim)) + 1
ny = int(inp("n", sim)) + 1
nz = int(inp("p", sim)) + 1
dt = inp("dt", sim)
D_s = inp("surf_diff", sim)
z0 = inp("z_domain%beg", sim)
dz = (inp("z_domain%end", sim) - z0) / nz
dx = (inp("x_domain%end", sim) - inp("x_domain%beg", sim)) / nx

# z of every cell (Fortran block order: x fastest, then y, then z)
z_cc = z0 + (np.arange(nz) + 0.5) * dz
Z = np.repeat(z_cc, nx * ny)

files = [f for f in glob.glob(os.path.join(HERE, "restart_data", "lustre_*.dat")) if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])]
files.sort(key=lambda f: int(os.path.basename(f)[7:-4]))
ss = os.path.getsize(files[0]) // 8 // (nx * ny * nz)  # surfactant is the last conserved variable

ts, M1, M0 = [], [], []
for f in files:
    surf = np.fromfile(f, "<f8").reshape(ss, nx * ny * nz)[ss - 1]
    ts.append(int(os.path.basename(f)[7:-4]) * dt)
    M1.append((surf * Z).sum())
    M0.append(surf.sum())
ts, M1, M0 = np.array(ts), np.array(M1), np.array(M0)
r = M1 / M1[0]
w = (r > 0.15) & (r < 0.95)
rate = -np.polyfit(ts[w], np.log(r[w]), 1)[0]
exact = 2.0 * D_s / R**2  # l=1 -> l(l+1)=2

# one machine-readable line for the sweep: R/dx  measured  exact  drift%
print(f"{R/dx:.3f} {rate:.5f} {exact:.5f} {100*(M0[-1]/M0[0]-1):+.3f}")
