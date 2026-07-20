#!/usr/bin/env python3
# Measure the stretching-term response in the extensional-flow test. Drop stays centered at origin.
# Report vs time: total interfacial surfactant (mass; must be conserved), the m=2 mode a2 of the
# recovered concentration Gamma(theta)=Gamma_tilde/|grad c| on the band (a2>0 => surfactant concentrated
# at the x-tips theta=0,pi, i.e. the elongation axis -- the expected Stone&Leal trend), and the drop
# Taylor deformation D=(L-B)/(L+B) from the color function (drop elongating along x => D>0).
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
TH = np.arctan2(Y, X)
fs = sorted(
    [f for f in glob.glob(C + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
    key=lambda f: int(os.path.basename(f)[7:-4]),
)
ss = os.path.getsize(fs[0]) // 8 // (nx * ny)

print(f"{'t':>6} {'mass/m0':>8} {'a2(tips)':>9} {'deformD':>8}")
m0 = None
for f in fs:
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    surf, c = a[ss - 1], a[ss - 2]
    t = int(os.path.basename(f)[7:-4]) * dt
    gx = np.zeros_like(c); gy = np.zeros_like(c)
    gx[:, 1:-1] = c[:, 2:] - c[:, :-2]; gy[1:-1, :] = c[2:, :] - c[:-2, :]
    g = np.hypot(gx, gy)
    band = g > 0.05 * g.max()
    tot = surf.sum()
    if m0 is None:
        m0 = tot
    th = TH[band]
    G = surf[band] / g[band]
    B = np.column_stack([np.ones(band.sum()), np.cos(2 * th), np.sin(2 * th)])
    coef, *_ = np.linalg.lstsq(B, G, rcond=None)
    a2 = coef[1] / coef[0]                          # relative m=2 amplitude; >0 = tips-concentrated
    inside = c > 0.5
    L = (X[inside].max() - X[inside].min()) if inside.any() else 0.0   # x-extent
    Bax = (Y[inside].max() - Y[inside].min()) if inside.any() else 0.0  # y-extent
    D = (L - Bax) / (L + Bax) if (L + Bax) > 0 else 0.0
    print(f"{t:6.3f} {tot/m0:8.5f} {a2:9.4f} {D:8.4f}")
