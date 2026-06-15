#!/usr/bin/env python3
"""Diagnostic: localize WHERE the isothermal-Dirichlet wall BC reverses the Fig 7 drop.

Compares two runs at the same t* -- one with conduction + isothermal walls (reverses), one with
conduction + adiabatic walls (rises) -- by overlaying the independent temperature scalar T_s, the
vertical velocity v_y, and the drop interface (color c = 0.5). If the BC distorts the surface
temperature gradient that drives the Marangoni stress, it will show up as a different T_s pattern
around the drop and a reversed v_y near/through the drop.

Usage: python3 diag_isobc.py <iso_rundir> <adiab_rundir> [t_star]
"""

import glob
import os
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

iso_dir, adiab_dir = sys.argv[1], sys.argv[2]
t_star_target = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5


def load(case_dir, t_star):
    p = {}
    with open(os.path.join(case_dir, "simulation.inp")) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                p[k.strip().lower()] = v.strip().rstrip(",")
    f = lambda k: float(p[k.lower()])  # noqa: E731
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    dt = f("dt")
    mu_b, sigma_T = 1.0 / f("fluid_pp(1)%re(1)"), f("sigma_dtdt")
    Ly = f("y_domain%end") - f("y_domain%beg")
    gradT = 1.0 / Ly
    t_r = mu_b / abs(sigma_T * gradT)
    rd = os.path.join(case_dir, "restart_data")
    cells = nx * ny * nz
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    # snapshot nearest the requested t*
    s = min(steps, key=lambda s: abs(s * dt / t_r - t_star))
    snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
    nvars = snap.size // cells
    fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)[0]  # 2D slice  # noqa: E731
    rho = fld(0) + fld(1)
    vy = fld(3) / rho
    c = np.clip(fld(nvars - 2), 0.0, 1.0)
    Ts = fld(nvars - 1)
    xb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
    yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
    x = 0.5 * (xb[:-1] + xb[1:])
    y = 0.5 * (yb[:-1] + yb[1:])
    return dict(x=x, y=y, vy=vy, c=c, Ts=Ts, t=s * dt / t_r, ycen=(c * y[None, :].T).sum() / c.sum())


iso = load(iso_dir, t_star_target)
adi = load(adiab_dir, t_star_target)
print(f"isothermal  t*={iso['t']:.2f}  drop y_centroid={iso['ycen']:.4f}  max|vy|={np.abs(iso['vy']).max():.4e}")
print(f"adiabatic   t*={adi['t']:.2f}  drop y_centroid={adi['ycen']:.4f}  max|vy|={np.abs(adi['vy']).max():.4e}")

fig, axes = plt.subplots(2, 3, figsize=(13, 9))
for row, (tag, d) in enumerate([("isothermal (reverses)", iso), ("adiabatic (rises)", adi)]):
    X, Y = np.meshgrid(d["x"], d["y"])
    # T_s field + drop contour
    ax = axes[row, 0]
    pc = ax.pcolormesh(X, Y, d["Ts"], cmap="inferno", shading="auto")
    ax.contour(X, Y, d["c"], levels=[0.5], colors="cyan", linewidths=1.5)
    fig.colorbar(pc, ax=ax, fraction=0.046)
    ax.set_title(f"{tag}\nT_s  (t*={d['t']:.2f})")
    ax.set_aspect("equal")
    # v_y field + drop contour
    ax = axes[row, 1]
    vmax = max(np.abs(iso["vy"]).max(), np.abs(adi["vy"]).max())
    pc = ax.pcolormesh(X, Y, d["vy"], cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.contour(X, Y, d["c"], levels=[0.5], colors="k", linewidths=1.0)
    fig.colorbar(pc, ax=ax, fraction=0.046)
    ax.set_title("v_y  (red=up, blue=down)")
    ax.set_aspect("equal")
    # T_s vertical profile through the drop center column
    ax = axes[row, 2]
    jc = np.argmin(np.abs(d["x"]))  # center column
    ax.plot(d["Ts"][:, jc], d["y"], "r-", label="T_s (center column)")
    ax.axhspan(d["ycen"] - 0.5, d["ycen"] + 0.5, color="cyan", alpha=0.15, label="drop band")
    ax.set_xlabel("T_s")
    ax.set_ylabel("y")
    ax.set_title("T_s through drop center")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "diag_isobc_compare.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=130)
print(f"wrote {out}")
