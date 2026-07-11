#!/usr/bin/env python3
"""
Measure the m=1 surfactant surface-diffusion decay rate on the circle and compare to the exact
Laplace-Beltrami rate D_s/R^2 -- with TWO estimators of the mode amplitude, to expose the bias of the
whole-field moment on a curved (staircased) interface:
  - full-field: M1 = sum(Gamma_tilde * x)           (biased LOW -- plateaus below exact)
  - band-only:  M1 = sum(Gamma_tilde * x, |grad c|>) (biased HIGH)
The exact value sits BETWEEN the two, so the moment brackets rather than pins the rate. Prints one line
"R/dx  full-field  band-only  exact  drift%" for the convergence sweep.
Reads grid/dt/D_s from simulation.inp; R is the example constant below.
"""

import glob
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = 0.5  # circle radius (fixed in case.py)


def inp(name, path):
    m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([^\s,]+)", open(path).read(), re.M)
    if not m:
        raise KeyError(name)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


sim = os.path.join(HERE, "simulation.inp")
nx = int(inp("m", sim)) + 1
ny = int(inp("n", sim)) + 1
dt = inp("dt", sim)
D_s = inp("surf_diff", sim)
x0 = inp("x_domain%beg", sim)
dx = (inp("x_domain%end", sim) - x0) / nx
x_cc = x0 + (np.arange(nx) + 0.5) * dx
X = np.tile(x_cc, ny).reshape(ny, nx)

files = [f for f in glob.glob(os.path.join(HERE, "restart_data", "lustre_*.dat")) if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])]
files.sort(key=lambda f: int(os.path.basename(f)[7:-4]))
ss = os.path.getsize(files[0]) // 8 // (nx * ny)  # surfactant = last, color = second-to-last

ts, Mf, Mb, M0 = [], [], [], []
for f in files:
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    surf, c = a[ss - 1], a[ss - 2]
    gx = np.zeros_like(c)
    gy = np.zeros_like(c)
    gx[:, 1:-1] = c[:, 2:] - c[:, :-2]
    gy[1:-1, :] = c[2:, :] - c[:-2, :]
    band = np.sqrt(gx**2 + gy**2) > 0.05 * np.sqrt(gx**2 + gy**2).max()
    ts.append(int(os.path.basename(f)[7:-4]) * dt)
    Mf.append((surf * X).sum())
    Mb.append((surf * X * band).sum())
    M0.append(surf.sum())
ts = np.array(ts)


def rate_of(M):
    r = np.array(M) / M[0]
    w = (r > 0.15) & (r < 0.95)
    return -np.polyfit(ts[w], np.log(r[w]), 1)[0]


exact = D_s / R**2  # circle m=1 -> m^2 = 1
print(f"{R/dx:.3f} {rate_of(Mf):.5f} {rate_of(Mb):.5f} {exact:.5f} {100*(M0[-1]/M0[0]-1):+.3f}")
