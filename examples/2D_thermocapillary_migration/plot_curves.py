#!/usr/bin/env python3
"""Build the two curated VALIDATION curves into figures/, straight from the current runs:

  figures/tc1_fig5_rise_velocity_2D.png  -- TC1 / Samareh Fig 5(d): 2D zero-Marangoni (Ma=0) rise
                                            velocity v/v_YGB vs t/t_r, grid convergence (64/128/256), anchor 0.80.
  figures/tc2_fig7_migration_2D.png      -- TC2 / Samareh Fig 7: finite-Ma (Re=5, Ma=20, Ca=0.01666)
                                            migration U* vs t*, grids 64/128, vs the Nas & Tryggvason peak ~0.13.

Both curves are drawn as UNCONNECTED DOTS, one marker per saved snapshot. The runs are compressible
and the closed slip-wall box rings acoustically (the initial unbalanced Laplace jump sigma/r launches
standing waves), and we save only ~80-100 snapshots -- about 2.6 per acoustic period, right at the
Nyquist limit. Connecting those samples with lines turns the aliased acoustic oscillation into a
spurious sawtooth; plotting them as a dot cloud shows the data honestly (the eye averages it) and the
migration signal -- the mean of the cloud -- is what compares to Samareh. The acoustic ripple is set
by sound speed and box size, not by dx, so it does not shrink with grid refinement.

Reads everything run-dependent from each run's simulation.inp, so it can't silently disagree with the
data. Conserved layout (model_eqns=3, num_fluids=2): 0,1 = partial densities, 2 = x-mom, 3 = y-mom,
color c last (second-to-last when a thermal_scalar T_s is appended).

Usage:  python3 plot_curves.py
"""

import glob
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIGS = os.path.join(HERE, "figures")
GRADT = 2.0 / 15.0  # imposed |dT/dy|, common to both cases
R = 0.5  # drop radius (D = 1)

# seaborn theme + a colorblind-safe categorical palette. Each grid gets one distinct hue (no reuse);
# reference lines are drawn in neutral black/gray so they never collide with a data color.
sns.set_theme(style="whitegrid", palette="colorblind", font_scale=1.15)
PAL = sns.color_palette("colorblind")
INK = "0.15"  # near-black for reference lines


def read_namelist(path):
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
    c_idx = nvars - 2 if ts else nvars - 1
    t, u_lab = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        vy = fld(3) / (fld(0) + fld(1))
        c = np.clip(fld(c_idx), 0.0, None)
        t.append(s * f("dt"))
        u_lab.append((c * vy).sum() / c.sum())
    return np.array(t), np.array(u_lab), P


def fig5_rise_velocity():
    """Samareh Fig 5(d): v/v_YGB vs t/t_r, 2D grid convergence in the slip-wall box."""
    grids = [("fig5_2D_w064", 64, PAL[0], "o"), ("fig5_2D_w128", 128, PAL[1], "s"), ("fig5_2D_w256", 256, PAL[2], "^")]
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    plotted = False
    for name, nx, col, mk in grids:
        out = color_weighted_vy(os.path.join(RUNS, name))
        if out is None:
            continue
        t, u_lab, P = out
        mu = 1.0 / float(P["fluid_pp(1)%re(1)"])
        dsdt = float(P["sigma_dtdt"])
        t_r = mu / abs(dsdt * GRADT)
        v_YGB = (2.0 / 15.0) * (-dsdt) * GRADT * R / mu
        ax.plot(t / t_r, u_lab / v_YGB, mk, color=col, ms=6.5, mew=0, alpha=0.85, linestyle="none",
                label=f"MFC {nx} cells/width ({nx / 5:.0f}/$D$)")
        plotted = True
    if not plotted:
        print("  fig5: no runs found")
        return
    ax.axhline(1.0, ls=":", color=INK, lw=1.5, label=r"$v_{\mathrm{YGB}}$ (Samareh Eq. 29, ratio = 1)")
    ax.axhline(0.80, ls="--", color=INK, lw=1.5, label=r"Samareh 2D $\approx$ 0.80 (Fig 5, slip-wall box)")
    ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"$v / v_{\mathrm{YGB}}$")
    ax.set_title("Fig 5: 2D thermocapillary rise, grid convergence (Ma = 0)")
    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, 1.4)
    ax.legend(loc="lower right", fontsize=11, frameon=True, framealpha=0.92)
    sns.despine(ax=ax)
    out = os.path.join(FIGS, "tc1_fig5_rise_velocity_2D.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  wrote {out}")


def fig7_migration():
    """Samareh Fig 7: U* = U/U_r vs t* = t/t_r, finite-Ma (Re=5, Ma=20, Ca=0.01666), grids 64/128."""
    grids = [("fig7_w064", 64, PAL[0], "o"), ("fig7_w128", 128, PAL[3], "s")]
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    plotted = False
    for name, nx, col, mk in grids:
        out = color_weighted_vy(os.path.join(RUNS, name))
        if out is None:
            continue
        t, u_lab, P = out
        mu_b = 1.0 / float(P["fluid_pp(1)%re(1)"])
        sigma_T = float(P["sigma_dtdt"])
        Ly = float(P["y_domain%end"]) - float(P["y_domain%beg"])
        gradT = abs(float(P["bc_y%twall_out"]) - float(P["bc_y%twall_in"])) / Ly if "bc_y%twall_out" in P else 1.0 / Ly
        G = abs(sigma_T * gradT)  # Marangoni stress scale
        U_r, t_r = G * R / mu_b, mu_b / G
        ax.plot(t / t_r, u_lab / U_r, mk, color=col, ms=6.5, mew=0, alpha=0.85, linestyle="none",
                label=f"MFC {nx} cells/width ({nx / 2:.0f}/$D$)")
        plotted = True
    if not plotted:
        print("  fig7: no runs found")
        return
    ax.axhline(0.13, ls="--", color=INK, lw=1.5, label=r"Nas & Tryggvason peak $\approx$ 0.13")
    ax.axhline(0.0, ls="-", color="0.6", lw=1.0)
    ax.set_xlabel(r"$t^* = t / t_r$  ($t_r = \mu_b / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"$U^* = U / U_r$")
    ax.set_title("Fig 7: 2D finite-Ma migration (Re=5, Ma=20, Ca=0.0167)")
    ax.set_xlim(left=0.0)
    ax.set_ylim(-0.04, 0.185)
    ax.legend(loc="upper right", fontsize=11, frameon=True, framealpha=0.92)
    sns.despine(ax=ax)
    out = os.path.join(FIGS, "tc2_fig7_migration_2D.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  wrote {out}")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig5_rise_velocity()
    fig7_migration()
