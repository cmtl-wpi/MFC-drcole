#!/usr/bin/env python3
"""The validation overlay figures for TC1 (Fig 5) and TC2 (Fig 7), in Samareh's published style.

Each trace is drawn as raw per-snapshot MARKERS joined by a THIN DASHED line, with NO smoothing.
MFC is compressible, so the closed box rings acoustically (an aliased acoustic standing wave, not
migration -- Samareh's incompressible solver has no acoustics); that ring shows up directly as the
line's small zig-zag, and we leave it visible. Everything is drawn in Samareh's plain published
style (white plate, full box frame, outward ticks, no grid) over his digitized reference curves.

  figures/case1_fig5_samareh_style.png  (TC1)
      v/v_YGB vs t/t_r on [0,10]. MFC with bulk conduction (Ma=0.1, which holds T near-linear so the
      curve sits at a FLAT plateau like Samareh's Ma=0 invariant-T limit) vs Samareh Fig 5(d) VOF.
  figures/case2_fig7_samareh_style.png  (TC2)
      U*=U/U_r vs t*=t/t_r on [0,20]. MFC (64, 128 cells/width) vs the digitized Nas & Tryggvason
      transient.

Every run-dependent constant is read from each run's simulation.inp (via color_weighted_vy), so the
figures can't silently disagree with the data. Conserved layout (model_eqns=3, num_fluids=2): 0,1 =
partial densities, 3 = y-momentum, color c last (second-to-last when a thermal_scalar T_s is appended).

Usage:  python3 plot_samareh_style.py
"""

import glob
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figures")
RUNS = os.path.join(HERE, "runs")
R = 0.5  # drop radius (D = 1)
GRADT = 2.0 / 15.0  # imposed |dT/dy|, common to TC1/TC2


def read_namelist(path):
    """Parse a Fortran namelist file's plain "name = value" lines into a dict (lowercase keys)."""
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def color_weighted_vy(run_dir):
    """Per-snapshot color-weighted lab-frame y-velocity history of a slip-wall run.
    Returns (t, u_lab, params) or None. The drop migrates in +y, so u_lab IS the rise velocity."""
    inp = os.path.join(run_dir, "simulation.inp")
    rd = os.path.join(run_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    P = read_namelist(inp)
    f = lambda k: float(P[k.lower()])  # noqa: E731
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    ts = str(P.get("thermal_scalar", "F")).strip(". ").upper().startswith("T")
    cells = nx * ny * nz
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    c_idx = nvars - 2 if ts else nvars - 1  # color function (T_s appended after it in ts mode)
    t, u_lab = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        vy = fld(3) / (fld(0) + fld(1))
        c = np.clip(fld(c_idx), 0.0, None)
        t.append(s * f("dt"))
        u_lab.append((c * vy).sum() / c.sum())
    return np.array(t), np.array(u_lab), P


def v_ygb_ratio(out):
    """(t/t_r, v/v_YGB) from a color_weighted_vy result, using each run's own constants."""
    t, u_lab, P = out
    mu = 1.0 / float(P["fluid_pp(1)%re(1)"])
    dsdt = float(P["sigma_dtdt"])
    t_r = mu / abs(dsdt * GRADT)
    v_YGB = (2.0 / 15.0) * (-dsdt) * GRADT * R / mu
    return t / t_r, u_lab / v_YGB


# Nas & Tryggvason U*(t*) transient, digitized BY EYE from Samareh Fig 7 (the red open triangles; the
# two Samareh grids nearly coincide with it). Accuracy ~ +/-0.005 in U*. Anchors match the paper text:
# broad peak ~0.131 at t*~4-5, terminal ~0.10 at t*=20 (the fine grid is within 1.7% of N&T).
NAS_TRYGGVASON = np.array([
    (0.0, 0.0), (1.0, 0.055), (2.0, 0.100), (3.0, 0.122), (4.0, 0.130), (5.0, 0.131), (6.0, 0.128),
    (7.0, 0.124), (8.0, 0.120), (10.0, 0.114), (12.0, 0.110), (14.0, 0.106), (16.0, 0.103),
    (18.0, 0.101), (20.0, 0.0995)])

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
        out = color_weighted_vy(run)
        if out is None:
            print(f"  fig5: {name} unreadable, skipping")
            continue
        x, y = v_ygb_ratio(out)
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
        dst = os.path.join(FIGS, "case1_fig5_samareh_style.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(series)} grids: " + ", ".join(f"{c:.1f}/D" for *_, c in series) + ")")


def fig7_tc2():
    """TC2, Fig 7 (Re=5, Ma=20, Ca=0.01666): MFC migration vs the digitized Nas & Tryggvason transient."""
    with plt.rc_context(SAMAREH_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(NAS_TRYGGVASON[:, 0], NAS_TRYGGVASON[:, 1], "^--", color="0.0", ms=6.5, mfc="none", mew=1.3, lw=1.0,
                zorder=5, label="Nas & Tryggvason (digitized)")
        plotted = False
        for name, nx, color in [("fig7_w064", 64, "#4C72B0"), ("fig7_w128", 128, "#DD8452")]:
            out = color_weighted_vy(os.path.join(RUNS, name))
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
        dst = os.path.join(FIGS, "case2_fig7_samareh_style.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig5_tc1()
    fig7_tc2()
