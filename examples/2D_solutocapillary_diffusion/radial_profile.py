#!/usr/bin/env python3
# Diagnose the surface-diffusion operator's normal (off-interface) leakage on the circle: plot the
# azimuthally-averaged surfactant density Gamma_tilde vs radius, for a diffusion-ON run and a
# diffusion-OFF control. A correct surface-diffusion operator keeps Gamma_tilde a thin ring at r=R
# (only evening out AROUND the circle); leakage makes the ring spread radially into the bulk.
#
# Reproduce (from repo root), then run this with the two output dirs:
#   MFC_NX=128 MFC_DS=0.2 ./mfc.sh run examples/2D_solutocapillary_diffusion/case.py --no-debug -n 4 -t pre_process simulation
#   # ...move its restart_data aside, then the control:
#   MFC_NX=128 MFC_DS=0   ./mfc.sh run examples/2D_solutocapillary_diffusion/case.py --no-debug -n 4 -t pre_process simulation
#   python3 examples/2D_solutocapillary_diffusion/radial_profile.py <ON_dir> <OFF_dir>
import glob
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = 0.5


def load(case):
    def inp(name):
        m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([^\s,]+)", open(os.path.join(case, "simulation.inp")).read(), re.M)
        return float(m.group(1).replace("d", "e").replace("D", "e"))

    nx, ny = int(inp("m")) + 1, int(inp("n")) + 1
    dt, x0 = inp("dt"), inp("x_domain%beg")
    dx = (inp("x_domain%end") - x0) / nx
    xcc = x0 + (np.arange(nx) + 0.5) * dx
    r = np.hypot(np.tile(xcc, ny).reshape(ny, nx), np.repeat(xcc, nx).reshape(ny, nx))
    fs = sorted(
        [f for f in glob.glob(case + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
        key=lambda f: int(os.path.basename(f)[7:-4]),
    )
    ss = os.path.getsize(fs[0]) // 8 // (nx * ny)
    return nx, ny, dt, r, fs, ss


def radial(case, i):
    nx, ny, dt, r, fs, ss = load(case)
    surf = np.fromfile(fs[i], "<f8").reshape(ss, ny, nx)[ss - 1]
    edges = np.linspace(0.2, 0.8, 60)
    rc = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(r.ravel(), edges) - 1
    m = (idx >= 0) & (idx < len(rc))
    prof = np.bincount(idx[m], weights=surf.ravel()[m], minlength=len(rc)) / np.maximum(np.bincount(idx[m], minlength=len(rc)), 1)
    return rc, prof, int(os.path.basename(fs[i])[7:-4]) * dt


on_dir, off_dir = sys.argv[1], sys.argv[2]
n = len(load(on_dir)[4])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
for i, ls in [(0, "-"), (n // 6, "--"), (n // 3, "-."), (n - 1, ":")]:
    rc, p, t = radial(on_dir, i)
    ax[0].plot(rc, p, ls, label=f"t={t:.2f}")
ax[0].set_title("D_s=0.2 (surface diffusion ON): band WIDENS = leakage")
for i, ls in [(0, "-"), (n - 1, ":")]:
    rc, p, t = radial(off_dir, i)
    ax[1].plot(rc, p, ls, label=f"t={t:.2f}")
ax[1].set_title("D_s=0 (control): band STATIC")
for a in ax:
    a.axvline(R, color="gray", lw=0.6)
    a.set_xlabel("radius r")
    a.legend()
ax[0].set_ylabel(r"azimuthally-averaged $\tilde\Gamma$")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "normal_leakage.png"), dpi=120)
print("wrote figures/normal_leakage.png")
