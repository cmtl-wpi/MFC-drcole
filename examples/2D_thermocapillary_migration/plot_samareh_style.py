#!/usr/bin/env python3
"""Direct, Samareh-style comparison figures for TC1 (Fig 5) and TC2 (Fig 7).

Unlike plot_curves.py (which fades the acoustic ring with a running mean), these
figures plot each trace's raw per-snapshot points as MARKERS joined by a THIN
DASHED line, with NO smoothing. MFC is compressible, so the closed box rings
acoustically (an aliased acoustic standing wave, not migration -- Samareh's
incompressible solver has no acoustics); here that ring shows up directly as the
line's small zig-zag. Everything is drawn in Samareh's plain published style
(white plate, full box frame, no grid) over his digitized reference curves.

  figures/case1_fig5_samareh_style.png  (TC1)
      v/v_YGB vs t/t_r on [0,10]. MFC with bulk conduction (Ma=0.1, which holds T
      near-linear so the curve sits at a FLAT plateau like Samareh's Ma=0
      invariant-T limit) vs Samareh Fig 5(d) VOF.
  figures/case2_fig7_samareh_style.png  (TC2)
      U*=U/U_r vs t*=t/t_r on [0,20]. MFC (64, 128 cells/width) vs the digitized
      Nas & Tryggvason transient.

Run-dependent constants are read from each run's simulation.inp via
plot_curves.color_weighted_vy, so the figures can't silently disagree with the
data. This script reuses color_weighted_vy / _v_ygb_ratio / NAS_TRYGGVASON from
plot_curves rather than duplicating them.

Usage:  python3 plot_samareh_style.py
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_curves as pc  # reuses color_weighted_vy / _v_ygb_ratio / NAS_TRYGGVASON

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figures")
RUNS = os.path.join(HERE, "runs")
R = 0.5  # drop radius (D = 1)

# Samareh Fig 5(d) VOF curve (sharp-interface analogue of MFC), digitized by eye from the published
# raster (~ +/-0.02 in v/v_YGB); his invariant-T plateau holds flat ~0.82-0.84 out to t/t_r = 10.
SAMAREH_VOF = np.array([
    (0.0, 0.0), (0.12, 0.42), (0.28, 0.70), (0.45, 0.80), (0.7, 0.815), (1.0, 0.82), (2.0, 0.825),
    (3.0, 0.83), (4.0, 0.83), (5.0, 0.835), (6.0, 0.83), (7.0, 0.835), (8.0, 0.838), (9.0, 0.835),
    (10.0, 0.84)])

# Samareh's plain plate style: white background, full box frame, outward ticks, no gridlines.
SAMAREH_STYLE = {
    "axes.grid": False, "axes.facecolor": "white", "figure.facecolor": "white",
    "axes.edgecolor": "0.0", "axes.linewidth": 1.0, "font.size": 13,
    "xtick.direction": "out", "ytick.direction": "out",
}
# RINGNOTE removed per request: do not annotate plots with the ring note.


def fig5_tc1():
    """TC1, Fig 5 (Ma=0 limit): MFC bulk-conduction runs at two grids vs Samareh's VOF curve.

    Both grids overshoot to a peak that climbs with refinement (toward/over Samareh's
    ~0.80) and then drift back down -- they do NOT hold Samareh's flat plateau, so the
    refinement trend in the peak is the honest comparison.
    """
    # (run dir, color). box width = 5D, so cells/D = (m+1)/5 -- computed per run, not hardcoded.
    runs = [
        ("tc1_cond_ma01_w64", "#C44E52"),    # ~12.8 cells/D
        ("fig5_cond_w128_tr10", "#4C72B0"),  # ~25.6 cells/D
    ]
    series = []
    for name, color in runs:
        run = os.path.join(RUNS, name)
        if not os.path.isdir(os.path.join(run, "restart_data")):
            print(f"  fig5: {name} not found, skipping")
            continue
        out = pc.color_weighted_vy(run)
        if out is None:
            print(f"  fig5: {name} unreadable, skipping")
            continue
        x, y = pc._v_ygb_ratio(out)
        cells_per_D = (int(out[2]["m"]) + 1) / 5.0
        series.append((x, y, color, cells_per_D))
    if not series:
        print("  fig5: no conduction runs found")
        return
    with plt.rc_context(SAMAREH_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(SAMAREH_VOF[:, 0], SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0, mfc="none", mew=1.3,
                lw=1.0, label=r"Samareh Fig 5(d), VOF ($Ma=0$, digitized)")
        for x, y, color, cpd in series:
            ax.plot(x, y, "o--", color=color, ms=4.0, mew=0, lw=0.9, alpha=0.8,
                    label=rf"MFC bulk conduction ($Ma=0.1$), {cpd:.1f}/$D$ — to $t/t_r={x.max():.1f}$")
        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel(r"Time   $t/t_r$")
        ax.set_ylabel(r"Normalized Rise Velocity   $v/v_{\mathrm{YGB}}$")
        ax.set_title("Fig 5 — 2D thermocapillary rise in the $Ma=0$ limit", fontsize=12)
        ax.legend(loc="lower right", fontsize=10, frameon=False)
        # note removed
        dst = os.path.join(FIGS, "case1_fig5_samareh_style.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(series)} grids: " + ", ".join(f"{c:.1f}/D" for *_, c in series) + ")")


def fig7_tc2():
    """TC2, Fig 7 (Re=5, Ma=20, Ca=0.01666): MFC migration vs the digitized Nas & Tryggvason transient."""
    nt = pc.NAS_TRYGGVASON
    with plt.rc_context(SAMAREH_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(nt[:, 0], nt[:, 1], "^--", color="0.0", ms=6.5, mfc="none", mew=1.3, lw=1.0,
                zorder=5, label="Nas & Tryggvason (digitized)")
        plotted = False
        for name, nx, color in [("fig7_w064", 64, "#4C72B0"), ("fig7_w128", 128, "#DD8452")]:
            out = pc.color_weighted_vy(os.path.join(RUNS, name))
            if out is None or len(out[0]) < 10:  # skip absent or depleted runs (need a real curve)
                if out is not None:
                    print(f"  fig7: skipping {name} -- only {len(out[0])} snapshots on disk (data depleted)")
                continue
            t, u_lab, params = out

            # Build the Marangoni reference scales (velocity U_r, time t_r) from this run's constants.
            mu_b = 1.0 / float(params["fluid_pp(1)%re(1)"])
            sigma_T = float(params["sigma_dtdt"])
            Ly = float(params["y_domain%end"]) - float(params["y_domain%beg"])
            if "bc_y%twall_out" in params:
                gradT = abs(float(params["bc_y%twall_out"]) - float(params["bc_y%twall_in"])) / Ly
            else:
                gradT = 1.0 / Ly
            marangoni_stress = abs(sigma_T * gradT)
            U_r = marangoni_stress * R / mu_b
            t_r = mu_b / marangoni_stress
            ts, us = t / t_r, u_lab / U_r

            ax.plot(ts, us, "o--", color=color, ms=3.5, mew=0, lw=0.9, alpha=0.6, zorder=3,
                    label=f"MFC {nx} cells/width ({nx // 2}/$D$)")
            plotted = True
        if not plotted:
            print("  fig7: no runs found")
            plt.close(fig)
            return
        ax.axhline(0.0, color="0.75", lw=0.8, zorder=1)  # rest baseline (raw scatter dips below it)
        ax.set_xlim(0.0, 20.0)
        ax.set_ylim(-0.025, 0.15)  # slightly past Samareh's 0-0.14 so the raw acoustic scatter is visible
        ax.set_xlabel(r"$t^* = t/t_r$")
        ax.set_ylabel(r"$U^* = U/U_r$")
        ax.set_title("Fig 7 — 2D migration at finite $Ma$ (Re=5, Ma=20, Ca=0.0167)", fontsize=12)
        ax.legend(loc="upper right", fontsize=10, frameon=False)
        # note removed
        dst = os.path.join(FIGS, "case2_fig7_samareh_style.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig5_tc1()
    fig7_tc2()
