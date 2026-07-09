#!/usr/bin/env python3
"""Grid-convergence overlay for the 3D thermocapillary rise -> figures/fig5_3d_overlay.png.

Draws the fixed-Ma=1 grid family (every runs/grid_ma/nx*/ma1 with restart data, coarse->fine, color =
Nx) on Samareh's digitized Fig 6 (3D VOF) axes. As dx->0 the plateau climbs toward Samareh's converged
~0.95 / the analytic zero-Ma YGB ceiling of 1.0. Reuses the 2D example's curve machinery
(color_weighted_vy + v_ygb_ratio + PLATE_STYLE) so 3D runs land on the same v/v_YGB-vs-t/t_r axes.
Unstable runs (max|v/v_YGB|>1.3) are drawn faded+dashed and flagged, not silently dropped. Skips
incomplete leaves, so it can be re-run while a grid point is still in flight.

Companion to plot_fig5_ma.py (the Marangoni-number sweep at fixed Nx). Usage: python3 plot_fig5.py
"""

import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TWO_D = os.path.abspath(os.path.join(HERE, "..", "2D_thermocapillary_migration"))
sys.path.insert(0, TWO_D)
from plot import PLATE_STYLE, color_weighted_vy, v_ygb_ratio  # noqa: E402

FIGS = os.path.join(HERE, "figures")
SWEEP = os.path.join(HERE, "runs", "grid_ma")

# Samareh Fig 6 (3D VOF), v/v_YGB vs t/t_r, digitized from the published page-5 right panel (same
# arrays as plot_fig5_ma.py): nx=128 (solid) plateaus ~0.945, nx=64 (dashed) ~0.928.
SAMAREH_FIG6_NX128 = np.array([(0.00, 0.00), (0.08, 0.55), (0.13, 0.74), (0.18, 0.83), (0.25, 0.88), (0.35, 0.91), (0.5, 0.925), (0.7, 0.935), (1.0, 0.94), (1.3, 0.945), (1.7, 0.944), (2.0, 0.943), (2.5, 0.945), (3.0, 0.947)])
SAMAREH_FIG6_NX64 = np.array([(0.00, 0.00), (0.08, 0.52), (0.13, 0.70), (0.18, 0.79), (0.25, 0.855), (0.35, 0.885), (0.5, 0.905), (0.7, 0.915), (1.0, 0.92), (1.3, 0.924), (1.7, 0.925), (2.0, 0.926), (2.5, 0.928), (3.0, 0.931)])


def curve_of(run_dir):
    """(x=t/t_r, y=v/v_YGB, plateau, unstable) or None if no data."""
    out = color_weighted_vy(run_dir)
    if out is None or len(out[0]) < 5:
        return None
    x, y = v_ygb_ratio(out)
    unstable = float(np.max(np.abs(y))) > 1.3
    win = y[(x >= 1.0) & (x <= 2.0)]
    plateau = float(np.median(win)) if win.size else float(y[-1])
    return x, y, plateau, unstable


def discover():
    """Completed (Nx, run_dir) at Ma=1, coarse->fine."""
    runs = {}
    for d in glob.glob(os.path.join(SWEEP, "nx*", "ma1")):
        m = re.match(r"nx(\d+)$", os.path.basename(os.path.dirname(d)))
        if m and os.path.isdir(os.path.join(d, "restart_data")):
            runs[int(m.group(1))] = d
    return sorted(runs.items())


def main():
    points = discover()
    if not points:
        sys.exit(f"no Ma=1 grid runs under {SWEEP}/nx*/ma1 yet")
    nxs = [nx for nx, _ in points]
    lo, hi = min(nxs), max(nxs)

    os.makedirs(FIGS, exist_ok=True)
    cmap = plt.get_cmap("viridis")
    printed = []
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
        ax.plot(SAMAREH_FIG6_NX128[:, 0], SAMAREH_FIG6_NX128[:, 1], "-", color="0.0", lw=1.8, zorder=2, label=r"Samareh Fig 6, 3D VOF $n_x=128$ (digitized)")
        ax.plot(SAMAREH_FIG6_NX64[:, 0], SAMAREH_FIG6_NX64[:, 1], "--", color="0.0", lw=1.6, zorder=2, label=r"Samareh Fig 6, 3D VOF $n_x=64$ (digitized)")
        ax.axhline(1.0, color="0.3", lw=1.1, ls=":", zorder=1, label=r"$v_{\mathrm{YGB}}$ (analytic, zero-Ma sphere)")

        for nx, d in points:
            c = curve_of(d)
            if c is None:
                print(f"  skip nx{nx} {d} (no data)")
                continue
            x, y, plateau, unstable = c
            frac = (nx - lo) / (hi - lo) if hi > lo else 0.85
            color = cmap(0.15 + 0.7 * frac)
            tag = f"Nx={nx} (cells/D={nx / 5:.1f}), Ma=1"
            if unstable:
                ax.plot(x, np.clip(y, -0.05, 1.1), "--", color=color, lw=1.0, alpha=0.4, zorder=3, label=f"{tag}  (unstable)")
            else:
                ax.plot(x, y, "-", marker="o", color=color, ms=3.2, lw=1.8, alpha=0.95, zorder=4, solid_capstyle="round", label=f"{tag}, plateau = {plateau:.2f}")
            printed.append((nx, plateau, unstable))

        ax.set_xlim(0.0, 2.3)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
        ax.set_ylabel(r"normalized rise velocity   $v / v_{\mathrm{YGB}}$")
        ax.set_title(r"3D thermocapillary rise: grid convergence vs Samareh Fig 6 (3D VOF)", fontsize=12, loc="left")
        ax.legend(loc="lower right", fontsize=7.5, frameon=False, ncol=1)
        dst = os.path.join(FIGS, "fig5_3d_overlay.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(printed)} curve(s))")
    print(f"\n  {'Nx':>4} {'plateau':>8}  note")
    for nx, pl, unstable in printed:
        print(f"  {nx:>4} {pl:>8.3f}  {'UNSTABLE' if unstable else ''}")


if __name__ == "__main__":
    main()
