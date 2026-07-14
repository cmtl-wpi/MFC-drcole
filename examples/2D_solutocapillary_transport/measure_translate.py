#!/usr/bin/env python3
# Measure interface confinement of the surfactant under uniform translation. The drop center is known
# kinematically (xc0 + U*t, periodic). Report vs time: total Gamma_tilde (mass), band radial width
# (std of r about the moving center, weighted by Gamma_tilde -- constant = confined, growing = smearing
# off the interface), off-band leakage fraction (Gamma_tilde where |grad c| is small), and the recovered
# cos(theta) pattern amplitude (should ride along unchanged).
import glob
import os
import re
import sys

import numpy as np

C = sys.argv[1]
R, W, U, xc0, amp = 0.5, 4.0, 1.0, -1.0, 0.5


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

print(f"{'t':>6} {'mass/m0':>8} {'bandwidth':>9} {'leak%':>7} {'cosAmp':>7}")
m0 = None
for f in fs:
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    surf, c = a[ss - 1], a[ss - 2]
    t = int(os.path.basename(f)[7:-4]) * dt
    xc = ((xc0 + U * t + W / 2) % W) - W / 2               # periodic drop center
    dxp = ((X - xc + W / 2) % W) - W / 2                    # periodic x-distance to center
    rloc = np.hypot(dxp, Y)
    gx = np.zeros_like(c); gy = np.zeros_like(c)
    gx[:, 1:-1] = c[:, 2:] - c[:, :-2]; gy[1:-1, :] = c[2:, :] - c[:-2, :]
    g = np.hypot(gx, gy)
    band = g > 0.05 * g.max()
    tot = surf.sum()
    if m0 is None:
        m0 = tot
    w = surf / tot
    mu = (rloc * w).sum()
    bw = np.sqrt((((rloc - mu) ** 2) * w).sum())
    leak = 100 * surf[~band].sum() / tot                    # Gamma_tilde sitting off the interface band
    # cos(theta) amplitude of recovered concentration on the band
    th = np.arctan2(Y[band], dxp[band])
    G = surf[band] / g[band]
    B = np.column_stack([np.ones(band.sum()), np.cos(th), np.sin(th)])
    coef, *_ = np.linalg.lstsq(B, G, rcond=None)
    camp = np.hypot(coef[1], coef[2]) / coef[0]
    print(f"{t:6.2f} {tot/m0:8.5f} {bw:9.4f} {leak:7.3f} {camp:7.4f}")
