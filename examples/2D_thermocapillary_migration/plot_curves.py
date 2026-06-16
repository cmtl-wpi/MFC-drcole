#!/usr/bin/env python3
"""Build the two curated VALIDATION curves into figures/, straight from the current runs:

  figures/case1_zero_marangoni_2D_fig5_rise_velocity.png -- TC1 / Samareh Fig 5(d): 2D zero-Marangoni
                                            (Ma=0) rise velocity v/v_YGB vs t/t_r, plotted on Samareh's
                                            FULL 0-10 window. The 64-cell grid is run to t/t_r=10 so the
                                            curve spans the paper's axis; the 128/256 grids are the
                                            shorter (t/t_r<=2) convergence runs that show the frozen-T
                                            upward drift growing with refinement. Anchors: 0.80 and 1.0.
  figures/case2_low_marangoni_nas_tryggvason_fig7.png    -- TC2 / Samareh Fig 7: finite-Ma (Re=5, Ma=20,
                                            Ca=0.01666) migration U* vs t* on the paper's 0-20 axis, with
                                            the Nas & Tryggvason transient (digitized from Fig 7) overlaid
                                            as the reference -- not just its peak -- so ramp, peak timing,
                                            and decline can all be compared.

The raw data is drawn as one marker per saved snapshot. The runs are compressible and the closed
slip-wall box rings acoustically (the initial unbalanced Laplace jump sigma/r launches standing
waves), and we save only ~80-100 snapshots -- about 2.6 per acoustic period, right at the Nyquist
limit. Connecting those RAW samples with lines turns the aliased acoustic oscillation into a spurious
sawtooth, so we never do that. Instead fig7 fades the raw markers and overlays a centered running
mean: the migration signal -- the mean of the cloud -- is what compares to Samareh, and a moving
average is a legitimate low-pass below the acoustic band (unlike connecting individual samples). The
acoustic ripple is set by sound speed and box size, not by dx, so it does not shrink with grid
refinement. fig5 still shows the bare dot cloud (the eye averages it).

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


def running_mean(y, w):
    """Centered moving average with a shrinking window at the ends (min one sample), so the early
    overshoot is not damped by edge padding. Window w spans several acoustic periods (~2.6 snapshots
    per period), so it averages out the aliased acoustic ripple while preserving the migration signal."""
    y = np.asarray(y, float)
    half = w // 2
    return np.array([y[max(0, i - half):min(len(y), i + half + 1)].mean() for i in range(len(y))])


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


def _v_ygb_ratio(out):
    """(t/t_r, v/v_YGB) from a color_weighted_vy result, using each run's own constants."""
    t, u_lab, P = out
    mu = 1.0 / float(P["fluid_pp(1)%re(1)"])
    dsdt = float(P["sigma_dtdt"])
    t_r = mu / abs(dsdt * GRADT)
    v_YGB = (2.0 / 15.0) * (-dsdt) * GRADT * R / mu
    return t / t_r, u_lab / v_YGB


def fig5_rise_velocity():
    """Samareh Fig 5(d): v/v_YGB vs t/t_r in the slip-wall box, on the paper's FULL 0-10 axis.

    The 64-cell grid is run to t/t_r=10 (fig5_2D_w064_tr10) so the MFC curve spans Samareh's window
    rather than stopping at the t/t_r~2 of the convergence runs. The 128/256 grids are those shorter
    runs (t/t_r<=2); they sit ABOVE the 64-grid and drift up toward/over v_YGB, the frozen-T (Ma=0)
    advection drift that grows with refinement -- so the coarse-grid landing on 0.80 is partly
    fortuitous, not a converged limit. Plotted as marker clouds (the compressible box rings; connecting
    samples would alias into a sawtooth), with a running-mean line on the long 64-grid run."""
    fig, ax = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    # 64-grid LONG run: spans 0-10. Faint raw cloud (acoustic ring) + bold centered running mean.
    long_out = color_weighted_vy(os.path.join(RUNS, "fig5_2D_w064_tr10"))
    plotted = False
    if long_out is not None:
        x, y = _v_ygb_ratio(long_out)
        ax.plot(x, y, "o", color=PAL[0], ms=4.0, mew=0, alpha=0.20, linestyle="none", zorder=2)
        ax.plot(x, running_mean(y, 7), "-", color=PAL[0], lw=2.6, zorder=5,
                label=rf"MFC 64 cells/width (13/$D$), to $t/t_r\approx{x.max():.1f}$")
        plotted = True
    # Finer grids: the shorter (t/t_r<=2) convergence runs -- show the upward drift with refinement.
    for name, nx, col, mk in [("fig5_2D_w128", 128, PAL[1], "s"), ("fig5_2D_w256", 256, PAL[2], "^")]:
        out = color_weighted_vy(os.path.join(RUNS, name))
        if out is None:
            continue
        x, y = _v_ygb_ratio(out)
        ax.plot(x, y, mk, color=col, ms=5.0, mew=0, alpha=0.65, linestyle="none",
                label=rf"MFC {nx} cells/width ({nx / 5:.0f}/$D$), $t/t_r\leq2$")
        plotted = True
    if not plotted:
        print("  fig5: no runs found")
        plt.close(fig)
        return
    # Quasi-steady window where the frozen-T IC still matches Samareh's invariant-T field; beyond it the
    # advected density-proxy T distorts and MFC drifts away (Samareh holds T invariant, so theirs stays flat).
    ax.axvspan(0.5, 3.0, color="0.85", alpha=0.35, zorder=0)
    ax.text(1.75, 0.06, "quasi-steady\ncomparison window", ha="center", va="bottom", fontsize=8.5, color="0.45")
    ax.annotate("frozen-$T$ drift (Ma = 0;\nSamareh holds $T$ invariant)", xy=(6.2, 0.66), xytext=(3.2, 0.40),
                fontsize=8.5, color="0.45", arrowprops=dict(arrowstyle="->", color="0.55", lw=1.0))
    ax.axhline(1.0, ls=":", color=INK, lw=1.5, label=r"$v_{\mathrm{YGB}}$ (Samareh Eq. 29, ratio = 1)")
    ax.axhline(0.80, ls="--", color=INK, lw=1.5, label=r"Samareh 2D $\approx$ 0.80 (Fig 5)")
    ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"$v / v_{\mathrm{YGB}}$")
    ax.set_title("Fig 5: 2D thermocapillary rise (Ma = 0), Samareh's $0$–$10$ window")
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, 1.25)
    ax.legend(loc="lower right", fontsize=10, frameon=True, framealpha=0.92)
    sns.despine(ax=ax)
    out = os.path.join(FIGS, "case1_zero_marangoni_2D_fig5_rise_velocity.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  wrote {out}")


# Nas & Tryggvason U*(t*) transient, digitized BY EYE from Samareh Fig 7 (the red open triangles; the
# two Samareh grids nearly coincide with it). Accuracy ~ +/-0.005 in U*. Anchors match the paper text:
# broad peak ~0.131 at t*~4-5, terminal ~0.10 at t*=20 (the fine grid is within 1.7% of N&T).
NAS_TRYGGVASON = np.array([
    (0.0, 0.0), (1.0, 0.055), (2.0, 0.100), (3.0, 0.122), (4.0, 0.130), (5.0, 0.131), (6.0, 0.128),
    (7.0, 0.124), (8.0, 0.120), (10.0, 0.114), (12.0, 0.110), (14.0, 0.106), (16.0, 0.103),
    (18.0, 0.101), (20.0, 0.0995)])


def fig7_migration():
    """Samareh Fig 7: U* = U/U_r vs t* = t/t_r, finite-Ma (Re=5, Ma=20, Ca=0.01666), grids 64/128,
    on the paper's 0-20 axis with the digitized Nas & Tryggvason transient overlaid as the reference."""
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
        ts, us = t / t_r, u_lab / U_r
        # faint raw cloud (one marker per snapshot) + a bold centered running mean = the migration signal
        ax.plot(ts, us, mk, color=col, ms=4.0, mew=0, alpha=0.22, linestyle="none", zorder=2)
        ax.plot(ts, running_mean(us, 11), "-", color=col, lw=2.6, zorder=4,
                label=f"MFC {nx} cells/width ({nx / 2:.0f}/$D$)")
        plotted = True
    if not plotted:
        print("  fig7: no runs found")
        plt.close(fig)
        return
    # Reference: the full Nas & Tryggvason transient (open triangles + faint connector), not just a peak.
    ax.plot(NAS_TRYGGVASON[:, 0], NAS_TRYGGVASON[:, 1], "-", color=INK, lw=1.0, alpha=0.45, zorder=5)
    ax.plot(NAS_TRYGGVASON[:, 0], NAS_TRYGGVASON[:, 1], "^", color=INK, ms=6.5, mfc="none", mew=1.4,
            zorder=6, label="Nas & Tryggvason (digitized, Samareh Fig 7)")
    ax.axhline(0.0, ls="-", color="0.6", lw=0.8)
    ax.text(0.99, 0.02, "MFC markers: raw snapshots · lines: running mean (w=11)   ·   N&T digitized $\\pm$0.005",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color="0.4")
    ax.set_xlabel(r"$t^* = t / t_r$  ($t_r = \mu_b / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"$U^* = U / U_r$")
    ax.set_title("Fig 7: 2D finite-Ma migration (Re=5, Ma=20, Ca=0.0167)")
    ax.set_xlim(0.0, 20.0)
    ax.set_ylim(-0.02, 0.16)
    ax.legend(loc="upper right", fontsize=10, frameon=True, framealpha=0.92)
    sns.despine(ax=ax)
    out = os.path.join(FIGS, "case2_low_marangoni_nas_tryggvason_fig7.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  wrote {out}")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig5_rise_velocity()
    fig7_migration()
