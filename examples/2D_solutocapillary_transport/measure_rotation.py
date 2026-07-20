#!/usr/bin/env python3
# Feasibility diagnostic for the solid-body-rotation M0 test. Separates the two error sources:
#   FLOW  : does MFC hold u=(-w*y, w*x)?  -> rel RMS velocity error vs ideal (in the drop neighbourhood),
#           plus drop circularity D (rigid rotation => D~0) and centroid drift.
#   TRANSPORT: does Gamma ride the interface?  -> m=1 phase vs the exact w*t, amplitude, and L2 error of
#           Gamma(theta) against the exact 1 + 0.5*cos(theta - w*t) on the interface band.
# If the FLOW error dominates, the test is flow-limited and cannot validate transport. Argv[1]=case dir.
import glob
import os
import re
import sys

import numpy as np

C = sys.argv[1]
w = 1.0
R = 0.6


def inp(n):
    m = re.search(rf"^\s*{re.escape(n)}\s*=\s*([^\s,]+)", open(os.path.join(C, "simulation.inp")).read(), re.M)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


nx, ny = int(inp("m")) + 1, int(inp("n")) + 1
dt = inp("dt")
x0, y0 = inp("x_domain%beg"), inp("y_domain%beg")
dx = (inp("x_domain%end") - x0) / nx
dy = (inp("y_domain%end") - y0) / ny
X = np.tile(x0 + (np.arange(nx) + 0.5) * dx, ny).reshape(ny, nx)
Y = np.repeat(y0 + (np.arange(ny) + 0.5) * dy, nx).reshape(ny, nx)
Rr = np.sqrt(X**2 + Y**2)
TH = np.arctan2(Y, X)
fs = sorted(
    [f for f in glob.glob(C + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
    key=lambda f: int(os.path.basename(f)[7:-4]),
)
ss = os.path.getsize(fs[0]) // 8 // (nx * ny)
u_ideal, v_ideal = -w * Y, w * X
near = Rr < 1.0  # drop neighbourhood: judge flow maintenance where it matters, not at the corrupt boundary

print(f"# sys_size={ss}  frames={len(fs)}  T_rot={2 * np.pi / w:.3f}")
print(f"{'t':>7} {'wt/2pi':>7} {'velerr':>8} {'D':>7} {'|cen|':>7} {'m1_amp':>7} {'phase/wt':>8} {'L2_Gam':>7} {'mass':>8}")
m0 = None
for f in fs:
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    rho = a[0] + a[1]
    u, v = a[2] / rho, a[3] / rho
    c, surf = a[ss - 2], a[ss - 1]
    t = int(os.path.basename(f)[7:-4]) * dt
    # FLOW: relative RMS velocity error vs ideal rigid rotation, in the drop neighbourhood
    num = np.sqrt(np.mean((u[near] - u_ideal[near]) ** 2 + (v[near] - v_ideal[near]) ** 2))
    den = np.sqrt(np.mean(u_ideal[near] ** 2 + v_ideal[near] ** 2))
    velerr = num / den
    # drop shape/centroid from the color function
    wgt = np.clip(c, 0, 1)
    wn = wgt / wgt.sum()
    mx, my = (X * wn).sum(), (Y * wn).sum()
    Ixx = ((X - mx) ** 2 * wn).sum()
    Iyy = ((Y - my) ** 2 * wn).sum()
    Ixy = ((X - mx) * (Y - my) * wn).sum()
    tr = Ixx + Iyy
    disc = np.sqrt(max(tr * tr / 4 - (Ixx * Iyy - Ixy**2), 0.0))
    aa, bb = np.sqrt(tr / 2 + disc), np.sqrt(max(tr / 2 - disc, 0.0))
    D = (aa - bb) / (aa + bb) if (aa + bb) > 0 else 0.0
    # TRANSPORT: Gamma = surf/|grad c| on the interface core; m=1 mode weighted by |grad c|
    gx = np.gradient(c, dx, axis=1)
    gy = np.gradient(c, dy, axis=0)
    normc = np.sqrt(gx**2 + gy**2)
    core = normc > 0.1 * normc.max()
    gam = surf[core] / normc[core]
    th, wt_ = TH[core], normc[core]
    A = 2 * np.sum(gam * np.cos(th) * wt_) / np.sum(wt_)
    B = 2 * np.sum(gam * np.sin(th) * wt_) / np.sum(wt_)
    amp = np.hypot(A, B)
    phase = np.arctan2(B, A) % (2 * np.pi)
    exact = 1.0 + 0.5 * np.cos(th - w * t)
    l2 = np.sqrt(np.sum((gam - exact) ** 2 * wt_) / np.sum(wt_))
    if m0 is None:
        m0 = surf.sum()
    ph_ratio = phase / (w * t) if t > 0.3 else float("nan")
    print(f"{t:7.3f} {w * t / (2 * np.pi):7.3f} {velerr:8.4f} {D:7.4f} {np.hypot(mx, my):7.4f} " f"{amp:7.4f} {ph_ratio:8.4f} {l2:7.4f} {surf.sum() / m0:8.5f}")
