#!/usr/bin/env python3
"""
Measure the flat-interface surfactant surface-diffusion mode-decay rate and compare to the exact
value D_s * k^2, then write figures/decay.png.

Reads grid/dt/D_s from the run's simulation.inp (never hard-coded) and the conserved surfactant
density (the last conserved variable) from the parallel-IO restart files. Run after:
    ./mfc.sh run case.py -n 1 -t pre_process simulation
    python3 measure.py
"""

import glob
import math
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def inp(name, path):
    txt = open(path).read()
    m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([^\s,]+)", txt, re.M)
    if not m:
        raise KeyError(name)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


sim = os.path.join(HERE, "simulation.inp")
nx = int(inp("m", sim)) + 1
ny = int(inp("n", sim)) + 1
dt = inp("dt", sim)
D_s = inp("surf_diff", sim)
Lx = inp("x_domain%end", sim) - inp("x_domain%beg", sim)
k = 2.0 * math.pi / Lx
dx = Lx / nx
x_cc = inp("x_domain%beg", sim) + (np.arange(nx) + 0.5) * dx
COS = np.tile(np.cos(k * x_cc), ny)
SIN = np.tile(np.sin(k * x_cc), ny)

files = [f for f in glob.glob(os.path.join(HERE, "restart_data", "lustre_*.dat")) if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])]
files.sort(key=lambda f: int(os.path.basename(f)[7:-4]))
ss = os.path.getsize(files[0]) // 8 // (nx * ny)  # surfactant is the last conserved variable

ts, A, M0 = [], [], []
for f in files:
    surf = np.fromfile(f, "<f8").reshape(ss, ny * nx)[ss - 1]
    ts.append(int(os.path.basename(f)[7:-4]) * dt)
    A.append(math.hypot((surf * COS).sum(), (surf * SIN).sum()))
    M0.append(surf.sum())
ts, A, M0 = np.array(ts), np.array(A), np.array(M0)
r = A / A[0]
mask = r > 0.15
rate = -np.polyfit(ts[mask], np.log(r[mask]), 1)[0]
ana = D_s * k**2

print(f"total surfactant drift : {100*(M0[-1]/M0[0]-1):+.3f}%")
print(f"measured decay rate    : {rate:.4f}")
print(f"exact  D_s k^2         : {ana:.4f}")
print(f"relative error         : {100*(rate/ana-1):+.2f}%")

fig, ax = plt.subplots(figsize=(6.2, 4.4))
tt = np.linspace(0, ts[-1], 200)
ax.semilogy(tt, np.exp(-ana * tt), "-", color="#333", lw=2, label=rf"exact  $e^{{-D_s k^2 t}}$   ($D_s k^2={ana:.3f}$)")
ax.semilogy(ts, r, "o", ms=6, mfc="#d1495b", mec="k", mew=0.5, label=rf"MFC  (fit ${rate:.3f}$, err ${100*(rate/ana-1):+.1f}\%$)")
ax.set_xlabel("time  $t$")
ax.set_ylabel(r"mode amplitude  $A(t)/A(0)$")
ax.set_title("Insoluble-surfactant surface diffusion: flat-interface mode decay")
ax.grid(True, which="both", alpha=0.25)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "decay.png"), dpi=130)
print("wrote figures/decay.png")
