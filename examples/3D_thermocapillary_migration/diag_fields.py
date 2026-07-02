#!/usr/bin/env python3
"""Cell-resolved midplane (z = Nz/2) slices of a 3D thermocapillary run, to characterize the
Ma-instability's spatial structure: a box-filling acoustic standing wave vs an interface-localized
blow-up. Lead with the raw field; extract scalars after.

Plots vy (rise velocity), |v|, and the dilatation-free pressure proxy (total energy density, an
acoustic-wave tracer) on the x-y midplane, with the color interface c=0.5 contour overlaid.

Conserved layout (model_eqns=3, num_fluids=2, surface_tension; nvars=11): 0,1 = partial densities,
2,3,4 = x/y/z-momentum, 5 = total energy, 6,7 = volume fractions, 8,9 = phasic internal energies,
10 = color. All read from each run's simulation.inp so it can't silently mis-slice.

Usage:  python3 diag_fields.py <run_dir> [step]   (step defaults to last snapshot)
"""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_inp(path):
    out = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().lower()] = v.strip().rstrip(",")
    return out


def main():
    run = sys.argv[1]
    P = read_inp(os.path.join(run, "simulation.inp"))
    f = lambda k: float(P[k.lower()])
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    cells = nx * ny * nz
    rd = os.path.join(run, "restart_data")
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat"))
                   if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    step = int(sys.argv[2]) if len(sys.argv) > 2 else steps[-1]
    snap = np.fromfile(os.path.join(rd, f"lustre_{step}.dat"), np.float64)
    nvars = snap.size // cells

    def fld(i):
        return snap[i * cells:(i + 1) * cells].reshape(nz, ny, nx)

    xb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), np.float64)[-(nx + 1):]
    yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1):]
    x = 0.5 * (xb[:-1] + xb[1:]); y = 0.5 * (yb[:-1] + yb[1:])

    rho = fld(0) + fld(1)
    kz = nz // 2
    vy = (fld(3) / rho)[kz]
    vx = (fld(2) / rho)[kz]
    vmag = np.sqrt(vx**2 + vy**2)
    E = fld(5)[kz]            # total energy density -- acoustic/compression tracer
    c = np.clip(fld(nvars - 1)[kz], 0, 1)

    dt = f("dt"); t_r = 7.5
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), constrained_layout=True)
    for ax, data, title, cmap in [
        (axes[0], vy, r"$v_y$ (rise velocity)", "RdBu_r"),
        (axes[1], vmag, r"$|v|$", "viridis"),
        (axes[2], E, "total energy density (acoustic tracer)", "magma"),
    ]:
        vmaxa = np.percentile(np.abs(data), 99.5)
        if cmap == "RdBu_r":
            im = ax.pcolormesh(x, y, data, cmap=cmap, vmin=-vmaxa, vmax=vmaxa, shading="auto")
        else:
            im = ax.pcolormesh(x, y, data, cmap=cmap, shading="auto")
        ax.contour(x, y, c, levels=[0.5], colors="lime", linewidths=1.3)
        ax.set_aspect("equal"); ax.set_title(title, fontsize=10)
        ax.set_xlabel("x"); fig.colorbar(im, ax=ax, shrink=0.8)
    axes[0].set_ylabel("y (rise axis)")
    fig.suptitle(f"{os.path.basename(run)}  step={step}  t/t_r={step*dt/t_r:.3f}  (nx={nx} ny={ny} nz={nz})", fontsize=11)
    viz = os.path.join(run, "viz"); os.makedirs(viz, exist_ok=True)
    dst = os.path.join(viz, f"fields_{step:06d}.png")
    fig.savefig(dst, dpi=130); plt.close(fig)
    # scalars: peak velocity, and how much of it sits FAR from the drop (a box acoustic mode fills
    # the domain; an interface blow-up is localized near c=0.5).
    far = np.abs(y) > 0.75 * yb[-1]
    print(f"{os.path.basename(run)} step={step} t/t_r={step*dt/t_r:.3f}: "
          f"max|vy|={np.abs(vy).max():.4e}  max|v|={vmag.max():.4e}  "
          f"max|vy| in far-field rows={np.abs(vy[far]).max():.4e}")
    print(f"  wrote {dst}")


if __name__ == "__main__":
    main()
