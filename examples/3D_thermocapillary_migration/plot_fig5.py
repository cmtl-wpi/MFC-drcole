#!/usr/bin/env python3
"""Overlay the 3D thermocapillary grid-convergence family on Samareh's digitized Fig 5(d) VOF curve.

Reuses the 2D example's digitized Samareh data (SAMAREH_VOF) AND its exact curve machinery
(color_weighted_vy + v_ygb_ratio + PLATE_STYLE) so every 3D run lands on the same v/v_YGB-vs-t/t_r
axes as the 2D validation. color_weighted_vy is dimension-agnostic (reshapes to (nz, ny, nx); vy is
the y-momentum index 3 in both 2D and 3D), so the 3D restart data plots through the identical path.

This overlays the grid-convergence family at one or more fixed Marangoni numbers (default Ma=1.0 and
Ma=0.5, both stable). Within each Ma family every completed runs/grid_ma/nx<NNN>/ma<tok> leaf is drawn
coarse->fine (color = Nx); the Ma families are distinguished by marker+linestyle. The runs/fig6_anchor/
nx064 is the (Nx=64, Ma=1.0) reference. As dx->0 the plateau climbs toward Samareh's converged ~0.95.
Numerically unstable runs (max|v/v_YGB|>1.3) are drawn faded+dashed and flagged, not silently dropped.
Skips not-yet-complete leaves, so it can be re-run while a sweep is still in flight.

Usage:  python3 plot_fig5.py            # overlay the Ma=1.0 and Ma=0.5 grid families
        python3 plot_fig5.py <run_dir>  # single-run mode (one curve, as before)
"""

import glob
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
from plot_grid_ma import is_complete  # noqa: E402  -- shared completeness gate (>=20 snaps OR [100%])

FIGS = os.path.join(HERE, "figures")
SWEEP = os.path.join(HERE, "runs", "grid_ma")
ANCHOR = os.path.join(HERE, "runs", "fig6_anchor", "nx064")

# (dir-token, Ma value, linestyle, marker) per grid-convergence family. Tokens match sweep_grid_ma fmt.
MA_FAMILIES = [("1", 1.0, "-", "o"), ("0p5", 0.5, "--", "^")]


def curve_of(run_dir):
    """(x=t/t_r, y=v/v_YGB, plateau, unstable) for a run, or None if no data."""
    out = color_weighted_vy(run_dir)
    if out is None or len(out[0]) < 5:
        return None
    x, y = v_ygb_ratio(out)
    unstable = float(np.max(np.abs(y))) > 1.3
    win = y[(x >= 1.0) & (x <= 2.0)]
    plateau = float(np.median(win)) if win.size else float(y[-1])
    return x, y, plateau, unstable


def discover(token):
    """Completed grid points (Nx, run_dir) for one Ma family, coarse->fine."""
    runs = {}
    for d in glob.glob(os.path.join(SWEEP, "nx*", f"ma{token}")):
        nx_tok = os.path.basename(os.path.dirname(d))[2:]
        if nx_tok.isdigit() and is_complete(d):
            runs[int(nx_tok)] = d
    if token == "1" and 64 not in runs and is_complete(ANCHOR):
        runs[64] = ANCHOR
    return sorted(runs.items())


def main():
    single = len(sys.argv) > 1
    if single:
        families = [("", None, "-", "o", [(None, sys.argv[1])])]
    else:
        families = [(tok, ma, ls, mk, discover(tok)) for tok, ma, ls, mk in MA_FAMILIES]
    all_nx = [nx for *_, pts in families for nx, _ in pts if nx is not None]
    if not all_nx and not single:
        sys.exit(f"no completed grid runs under {SWEEP} (run the sweep first)")
    lo, hi = (min(all_nx), max(all_nx)) if all_nx else (0, 1)

    os.makedirs(FIGS, exist_ok=True)
    cmap = plt.get_cmap("viridis")
    printed = []
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
        ax.plot(SAMAREH_VOF[:, 0], SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0, mfc="none", mew=1.3, lw=1.0, label=r"Samareh Fig 5(d), VOF (2D planar, digitized)")
        ax.axhline(1.0, color="0.3", lw=1.1, ls=":", zorder=1, label=r"$v_{\mathrm{YGB}}$ (analytic, zero-Ma sphere)")
        ax.axhline(0.95, color="#d62728", lw=1.1, ls="--", zorder=1, label=r"Samareh 3D converged $\approx 0.95$")

        for tok, ma, ls, mk, points in families:
            for nx, d in points:
                c = curve_of(d)
                if c is None:
                    print(f"  skip {d} (no data)")
                    continue
                x, y, plateau, unstable = c
                frac = (nx - lo) / (hi - lo) if (nx is not None and hi > lo) else 0.85
                color = cmap(0.15 + 0.7 * frac)
                tag = f"Nx={nx} (cells/D={nx / 5:.1f})" if nx is not None else "MFC 3D"
                tag += f", Ma={ma:g}" if ma is not None else ""
                if unstable:
                    ax.plot(x, np.clip(y, -0.05, 1.1), ls, color=color, lw=1.0, alpha=0.4, label=f"{tag}  (unstable)")
                else:
                    ax.plot(x, y, ls, marker=mk, color=color, ms=3.2, lw=1.8, alpha=0.95, solid_capstyle="round", label=f"{tag}, plateau = {plateau:.2f}")
                printed.append((ma, nx, plateau, unstable))

        ax.set_xlim(0.0, 2.3)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
        ax.set_ylabel(r"normalized rise velocity   $v / v_{\mathrm{YGB}}$")
        ax.set_title(r"3D thermocapillary rise: grid convergence vs Samareh Fig 5(d)", fontsize=12, loc="left")
        ax.legend(loc="lower right", fontsize=7.5, frameon=False, ncol=1)
        dst = os.path.join(FIGS, "fig5_3d_overlay.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(printed)} curve(s))")
    for ma, nx, pl, unstable in sorted(printed, key=lambda t: (t[0] is None, -(t[0] or 0), t[1] or 0)):
        print(f"    Ma={ma}  Nx={nx}  plateau = {pl:.3f}  {'UNSTABLE' if unstable else ''}")


if __name__ == "__main__":
    main()
