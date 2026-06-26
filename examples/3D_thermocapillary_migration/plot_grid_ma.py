#!/usr/bin/env python3
"""Convergence figure for the 3D grid+Ma sweep: separate the grid effect from the Ma effect.

Left panel  -- GRID convergence at Ma=0.1: plateau v/v_YGB vs 1/Nx. If it extrapolates toward
               ~0.95 as 1/Nx->0, the 0.82 deficit is resolution (the expected cause).
Right panel -- Ma INDEPENDENCE at Nx=32: plateau v/v_YGB vs Ma (log axis). If flat, the
               temperature is already invariant and Marangoni number is NOT the lever.

Reuses the 2D example's curve machinery (color_weighted_vy + v_ygb_ratio); it's dimension-agnostic.
Reads runs/grid_ma/nx<NNN>/ma<tok>/ (written by sweep_grid_ma.py) plus the existing
runs/fig6_anchor/nx064 as the (Nx=64, Ma=1.0) reference. Skips not-yet-complete leaves, so it can
run while the sweep is in flight. Writes figures/grid_ma_convergence.png.
"""

import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "2D_thermocapillary_migration")))
from plot import PLATE_STYLE, color_weighted_vy, v_ygb_ratio  # noqa: E402

FIGS = os.path.join(HERE, "figures")
SWEEP = os.path.join(HERE, "runs", "grid_ma")


def is_complete(d):
    """Accept finished runs ([100%]) OR partials that already reached the plateau (>=20 snapshots).
    The plateau settles by t/t_r~0.25, so a run stopped at t/t_r~1 is fully usable; this still
    rejects the stunted early-partial garbage (<20 snaps)."""
    n = len(glob.glob(os.path.join(d, "restart_data", "lustre_[0-9]*.dat")))
    if n >= 20:
        return True
    mfc = os.path.join(d, "MFC.out")
    return os.path.isfile(mfc) and "[100%]" in open(mfc, errors="ignore").read()


def plateau_of(run_dir):
    """Quasi-steady plateau v/v_YGB = median over the flat region (t/t_r >= 0.5), or None if not
    ready. The 0.5 floor skips the rise/overshoot and works for both full (0-2) and stopped (0-1) runs."""
    if not is_complete(run_dir):
        return None
    out = color_weighted_vy(run_dir)
    if out is None or len(out[0]) < 5:
        return None
    x, y = v_ygb_ratio(out)
    if float(np.max(np.abs(y))) > 1.3:  # numerically unstable (Ma=0.1 runs oscillate to +-1.6)
        return None
    win = y[x >= 0.5]
    return float(np.median(win)) if win.size else float(y[-1])


def discover():
    """Map every completed leaf -> (Nx, Ma, plateau). Includes the fig6_anchor (64,1.0) reference."""
    pts = []
    for d in glob.glob(os.path.join(SWEEP, "nx*", "ma*")):
        nx = int(os.path.basename(os.path.dirname(d))[2:])
        ma = float(os.path.basename(d)[2:].replace("p", "."))
        pl = plateau_of(d)
        if pl is not None:
            pts.append((nx, ma, pl))
    ref = os.path.join(HERE, "runs", "fig6_anchor", "nx064")
    if os.path.isdir(os.path.join(ref, "restart_data")):
        pl = plateau_of(ref)
        if pl is not None:
            pts.append((64, 1.0, pl))
    return pts


def plot_curves():
    """Line plot of every completed run's rise-velocity history v/v_YGB vs t/t_r -> grid_ma_curves.png."""
    runs = []
    for d in sorted(glob.glob(os.path.join(SWEEP, "nx*", "ma*"))):
        if not is_complete(d):
            continue
        nx = int(os.path.basename(os.path.dirname(d))[2:])
        ma = float(os.path.basename(d)[2:].replace("p", "."))
        runs.append((nx, ma, d))
    ref = os.path.join(HERE, "runs", "fig6_anchor", "nx064")
    if os.path.isdir(os.path.join(ref, "restart_data")):
        runs.append((64, 1.0, ref))
    if not runs:
        return
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
        ax.axhline(0.95, color="#d62728", ls="--", lw=1.1, label=r"Samareh 3D $\approx 0.95$")
        ax.axhline(1.0, color="0.3", ls=":", lw=1.0, label=r"$v_{\mathrm{YGB}}$")
        for nx, ma, d in sorted(runs):
            out = color_weighted_vy(d)
            if out is None:
                continue
            x, y = v_ygb_ratio(out)
            unstable = float(np.max(np.abs(y))) > 1.3  # Ma=0.1 runs oscillate to +-1.6 (numerically unstable)
            if unstable:
                ax.plot(x, y, "-", lw=0.8, alpha=0.3, label=f"Nx={nx}, Ma={ma:g}  (UNSTABLE)")
            else:
                ax.plot(x, y, "-o", ms=2.5, lw=1.8, label=f"Nx={nx}, Ma={ma:g}  (cells/D={nx / 5:.1f})")
        ax.set_xlim(0, 2.0)
        ax.set_ylim(0, 1.1)  # clip the unstable oscillations so the clean Ma=1.0 plateaus read clearly
        ax.set_xlabel(r"$t / t_r$")
        ax.set_ylabel(r"$v / v_{\mathrm{YGB}}$")
        ax.set_title("3D rise-velocity curves -- grid+Ma sweep (Ma=0.1 unstable, faded)", fontsize=12, loc="left")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        fig.savefig(os.path.join(FIGS, "grid_ma_curves.png"), dpi=200)
        plt.close(fig)


def main():
    pts = discover()
    if not pts:
        sys.exit(f"no completed runs under {SWEEP} yet (run sweep_grid_ma.py first)")
    os.makedirs(FIGS, exist_ok=True)
    plot_curves()  # the rise-velocity line plot (primary), refreshed alongside the convergence scatter

    # Grid convergence at Ma=1.0 -- the Ma=0.1 runs are numerically unstable (excluded by plateau_of),
    # so Ma=1.0 is the usable stable series for the convergence trend.
    grid = sorted([(nx, pl) for nx, ma, pl in pts if abs(ma - 1.0) < 1e-9])
    maser = sorted([(ma, pl) for nx, ma, pl in pts if nx == 32])

    with plt.rc_context(PLATE_STYLE):
        fig, (axg, axm) = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)

        axg.axhline(0.95, color="#d62728", lw=1.1, ls="--", label=r"Samareh 3D $\approx 0.95$")
        axg.axhline(1.0, color="0.3", lw=1.0, ls=":", label=r"$v_{\mathrm{YGB}}$")
        if grid:
            nxs = np.array([nx for nx, _ in grid])
            pls = np.array([pl for _, pl in grid])
            axg.plot(1.0 / nxs, pls, "o-", color="#1f77b4", ms=7, lw=1.8)
            for nx, pl in grid:
                axg.annotate(f"  Nx={nx}\n  {pl:.2f}", (1.0 / nx, pl), fontsize=8, va="center")
        axg.set_xlabel(r"$1/N_x$  (grid spacing $\to 0$)")
        axg.set_ylabel(r"plateau  $v / v_{\mathrm{YGB}}$")
        # Reference: the reliable (64, 1.0) anchor (different Ma, but the trustworthy data point so far).
        ref = [pl for nx, ma, pl in pts if nx == 64 and abs(ma - 1.0) < 1e-9]
        if ref:
            axg.plot(1.0 / 64, ref[0], "D", color="0.4", ms=8, label=rf"Nx=64, Ma=1.0 ref = {ref[0]:.2f}")
            axg.annotate(f"  {ref[0]:.2f}", (1.0 / 64, ref[0]), fontsize=8, va="center")
        axg.set_xlim(left=0.0)
        axg.set_ylim(0.3, 1.02)
        axg.set_title("Grid convergence @ Ma=1.0 (Ma=0.1 unstable)", fontsize=11, loc="left")
        axg.legend(loc="lower right", fontsize=8, frameon=False)

        axm.axhline(0.95, color="#d62728", lw=1.1, ls="--", label=r"Samareh 3D $\approx 0.95$")
        if maser:
            mas = np.array([ma for ma, _ in maser])
            pls = np.array([pl for _, pl in maser])
            axm.plot(mas, pls, "s-", color="#2ca02c", ms=7, lw=1.8)
            for ma, pl in maser:
                axm.annotate(f"  {pl:.2f}", (ma, pl), fontsize=8, va="center")
        axm.set_xscale("log")
        axm.set_xlabel(r"Marangoni number $Ma$  (invariant-$T$ limit $\to 0$)")
        axm.set_ylabel(r"plateau  $v / v_{\mathrm{YGB}}$")
        axm.set_ylim(0.3, 1.02)
        axm.set_title("Ma independence @ Nx=32", fontsize=11, loc="left")
        axm.legend(loc="lower right", fontsize=8, frameon=False)

        dst = os.path.join(FIGS, "grid_ma_convergence.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(pts)} runs: grid@Ma0.1={len(grid)}, Ma@Nx32={len(maser)})")

    print(f"\n  {'Nx':>4} {'Ma':>7} {'plateau':>8}")
    for nx, ma, pl in sorted(pts):
        print(f"  {nx:>4} {ma:>7} {pl:>8.3f}")


if __name__ == "__main__":
    main()
