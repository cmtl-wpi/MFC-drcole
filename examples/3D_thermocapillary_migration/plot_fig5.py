#!/usr/bin/env python3
"""Overlay the 3D thermocapillary rise on Samareh's digitized Fig 5(d) VOF curve.

Reuses the 2D example's digitized Samareh data (SAMAREH_VOF) AND its exact curve machinery
(color_weighted_vy + v_ygb_ratio + PLATE_STYLE) so the 3D run lands on the same v/v_YGB-vs-t/t_r
axes as the 2D validation -- the only difference is the run dir fed in. color_weighted_vy is
dimension-agnostic (reshapes to (nz, ny, nx); vy is the y-momentum index 3 in both 2D and 3D), so
the 3D restart data plots through the identical path.

Note the comparison this draws: SAMAREH_VOF is Samareh's Fig 5(d), the 2D-planar VOF result
(plateau ~0.83). The 3D run here is at Ma=1.0 (finite), NOT the Ma=0 limit, so it sits near the
same ~0.82 -- a finite-Ma coincidence with the 2D number, not a zero-Ma 3D validation (that target
is ~0.95, reached only as Ma->0 + grid converges). The dotted line is the analytic YGB limit (1.0).

Usage:  python3 plot_fig5.py [run_dir]      (default: runs/fig6_anchor/nx064)
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# Pull the digitized Samareh data + curve math straight from the 2D example (the user's reference).
TWO_D = os.path.abspath(os.path.join(HERE, "..", "2D_thermocapillary_migration"))
sys.path.insert(0, TWO_D)
from plot import PLATE_STYLE, SAMAREH_VOF, color_weighted_vy, v_ygb_ratio  # noqa: E402

FIGS = os.path.join(HERE, "figures")
run_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "runs", "fig6_anchor", "nx064")


def main():
    out = color_weighted_vy(run_dir)
    if out is None:
        sys.exit(f"no restart data in {run_dir} -- run the 3D case first")
    x, y = v_ygb_ratio(out)
    win = y[(x >= 1.0) & (x <= 2.0)]
    plateau = float(np.median(win)) if win.size else float(y[-1])

    os.makedirs(FIGS, exist_ok=True)
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(SAMAREH_VOF[:, 0], SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0, mfc="none", mew=1.3, lw=1.0, label=r"Samareh Fig 5(d), VOF (2D planar, digitized)")
        ax.axhline(1.0, color="0.3", lw=1.1, ls=":", zorder=1, label=r"$v_{\mathrm{YGB}}$ (analytic, zero-Ma sphere)")
        ax.axhline(0.95, color="#d62728", lw=1.1, ls="--", zorder=1, label=r"Samareh 3D converged $\approx 0.95$")
        ax.plot(x, y, "-o", color="#1f77b4", ms=3.5, lw=1.8, alpha=0.95, solid_capstyle="round", label=rf"MFC 3D sphere (Nx=64, Ma=1.0), plateau = {plateau:.2f}")
        ax.set_xlim(0.0, max(2.3, float(x.max()) + 0.2))
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
        ax.set_ylabel(r"normalized rise velocity   $v / v_{\mathrm{YGB}}$")
        ax.set_title(r"3D thermocapillary rise vs Samareh Fig 5(d) (digitized)", fontsize=12, loc="left")
        ax.legend(loc="lower right", fontsize=9, frameon=False)
        dst = os.path.join(FIGS, "fig5_3d_overlay.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}")
    print(f"  MFC 3D plateau v/v_YGB = {plateau:.3f}  (Samareh Fig 5(d) VOF ~ 0.83; 3D converged ~ 0.95)")


if __name__ == "__main__":
    main()
