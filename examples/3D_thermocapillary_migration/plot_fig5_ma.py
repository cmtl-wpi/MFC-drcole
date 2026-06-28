#!/usr/bin/env python3
"""Overlay the 3D thermocapillary Marangoni-number series (fixed grid) on Samareh's Fig 5(d) axes.

Companion to plot_fig5.py (which overlays the GRID family at fixed Ma=1.0). This holds the grid fixed
at the stable Nx=32 and sweeps Ma downward -- smaller Ma = stronger conduction = T held closer to the
imposed gradient = Samareh's invariant-T limit (target plateau ~0.95). Reuses the 2D example's curve
machinery (color_weighted_vy + v_ygb_ratio + PLATE_STYLE) so every run lands on the same
v/v_YGB-vs-t/t_r axes as the validation overlay.

Discovers runs/grid_ma/nx<NX>/ma<tok>/ for the chosen NX (default 32). Numerically unstable runs
(max|v/v_YGB| > 1.3 -- the conduction/interface coupling blows up as Ma drops on under-resolved grids)
are drawn faded + dashed and flagged "(unstable)", NOT silently dropped, so the fragility is visible.

Usage:  python3 plot_fig5_ma.py [NX]     (default NX=32)
"""

import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TWO_D = os.path.abspath(os.path.join(HERE, "..", "2D_thermocapillary_migration"))
sys.path.insert(0, TWO_D)
from plot import PLATE_STYLE, SAMAREH_VOF, color_weighted_vy, v_ygb_ratio  # noqa: E402
from plot_grid_ma import is_complete  # noqa: E402

FIGS = os.path.join(HERE, "figures")
SWEEP = os.path.join(HERE, "runs", "grid_ma")
NX = int(sys.argv[1]) if len(sys.argv) > 1 else 32


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
    """Completed Ma points at the chosen NX: list of (Ma, run_dir), high Ma -> low Ma."""
    runs = []
    for d in glob.glob(os.path.join(SWEEP, f"nx{NX:03d}", "ma*")):
        if is_complete(d):
            ma = float(os.path.basename(d)[2:].replace("p", "."))
            runs.append((ma, d))
    return sorted(runs, reverse=True)


def main():
    points = discover()
    if not points:
        sys.exit(f"no completed runs under {SWEEP}/nx{NX:03d}/ma* yet")

    os.makedirs(FIGS, exist_ok=True)
    cmap = plt.get_cmap("plasma")
    mas = [ma for ma, _ in points]
    # color by log(Ma): bright/yellow = low Ma (strong conduction), dark/purple = high Ma
    lg = np.log10(mas)
    lo, hi = min(lg), max(lg)
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(SAMAREH_VOF[:, 0], SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0, mfc="none", mew=1.3, lw=1.0, label=r"Samareh Fig 5(d), VOF (2D planar, digitized)")
        ax.axhline(1.0, color="0.3", lw=1.1, ls=":", zorder=1, label=r"$v_{\mathrm{YGB}}$ (analytic, zero-Ma sphere)")
        ax.axhline(0.95, color="#d62728", lw=1.1, ls="--", zorder=1, label=r"Samareh 3D converged $\approx 0.95$")

        printed = []
        for ma, d in points:
            c = curve_of(d)
            if c is None:
                print(f"  skip {d} (no data)")
                continue
            x, y, plateau, unstable = c
            frac = (np.log10(ma) - lo) / (hi - lo) if hi > lo else 0.5
            color = cmap(0.1 + 0.8 * (1.0 - frac))  # low Ma -> bright
            if unstable:
                ax.plot(x, np.clip(y, -0.05, 1.1), "--", color=color, lw=1.1, alpha=0.45, label=rf"MFC Ma={ma:g}  (unstable, $\max|v/v_{{YGB}}|>1.3$)")
            else:
                ax.plot(x, y, "-o", color=color, ms=3.0, lw=1.9, alpha=0.95, solid_capstyle="round", label=rf"MFC Ma={ma:g}, plateau = {plateau:.2f}")
            printed.append((ma, plateau, unstable))

        ax.set_xlim(0.0, 2.3)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
        ax.set_ylabel(r"normalized rise velocity   $v / v_{\mathrm{YGB}}$")
        ax.set_title(rf"3D thermocapillary rise: Marangoni-number sweep @ Nx={NX}", fontsize=12, loc="left")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        dst = os.path.join(FIGS, "fig5_3d_ma_overlay.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(printed)} curve(s))")
    print(f"\n  {'Ma':>7} {'plateau':>8}  note")
    for ma, pl, unstable in sorted(printed, reverse=True):
        print(f"  {ma:>7g} {pl:>8.3f}  {'UNSTABLE' if unstable else ''}")


if __name__ == "__main__":
    main()
