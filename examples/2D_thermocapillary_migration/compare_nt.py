#!/usr/bin/env python3
"""Overlay MFC against Nas & Tryggvason (2003) Figs 2 and 3 -- the head-to-head comparison the
validation actually calls for. Two figures:

  figures/nt_fig2_comparison.png  creeping single drop (Re=Ma=2.5e-3, Ca=1e-3) vs their Fig 2
  figures/nt_fig3_comparison.png  finite-Re drop (Re=5, Ma=20, Ca=0.0167)   vs their Fig 3

The paper curves are digitized by eye from the published plots (pages 1123-1124) -- approximate, for
visual comparison only, not exact data. MFC's Fig 3 curves are measured live from the case_Ma_20.py
runs in runs/tc2/{w064,w128}; MFC's Fig 2 curve is read from results/nt_fig2_mfc.json (the creeping
Nx=32 run; its restart_data was deleted, see that file's _provenance).

Usage:  python3 compare_nt.py
"""

import glob
import json
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figures")
R = 0.5

# --- digitized-by-eye reference curves from the paper (approximate) --------------------------------
# Fig 2: V* vs t*, converged ~64/D curve rises monotonically to a ~0.133-0.14 plateau (16/D lower).
NT_FIG2_TSTAR = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
NT_FIG2_VSTAR = [0.0, 0.025, 0.055, 0.085, 0.105, 0.118, 0.130, 0.137, 0.140, 0.141]
NT_FIG2_PLATEAU = (0.133, 0.141)  # 16/D .. 64/D converged band
# Fig 3: V* vs t*, converged 128x256 curve overshoots to ~0.129 at t*~5 then settles to ~0.10.
NT_FIG3_TSTAR = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0, 27.0]
NT_FIG3_VSTAR = [0.0, 0.045, 0.095, 0.118, 0.127, 0.129, 0.128, 0.122, 0.116, 0.109, 0.104, 0.101, 0.099, 0.098]
# ---------------------------------------------------------------------------------------------------


def read_namelist(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                name, value = line.split("=", 1)
                out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def measure_vstar(wd):
    """Color-weighted drop rise velocity U* = U/U_r vs t* = t/t_r from a run's restart_data."""
    p = read_namelist(os.path.join(wd, "simulation.inp"))

    def P(n):
        return float(p[n.lower()])

    nx, ny = int(P("m")) + 1, int(P("n")) + 1
    dt = P("dt")
    Ly = P("y_domain%end") - P("y_domain%beg")
    mu = 1.0 / P("fluid_pp(1)%re(1)")
    sigma_T = P("sigma_dtdt")
    gradT = abs(P("bc_y%twall_out") - P("bc_y%twall_in")) / Ly
    G = abs(sigma_T * gradT)
    U_r, t_r = G * R / mu, mu / G

    rd = os.path.join(wd, "restart_data")
    steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
    cells = nx * ny
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    c_idx = nvars - 1

    def fld(snap, i):
        return snap[i * cells : (i + 1) * cells].reshape(ny, nx)

    t_star, v_star = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        rho = fld(snap, 0) + fld(snap, 1)
        c = np.clip(fld(snap, c_idx), 0.0, 1.0)
        t_star.append(s * dt / t_r)
        v_star.append((c * fld(snap, 3) / rho).sum() / c.sum() / U_r)
    return np.array(t_star), np.array(v_star)


def fig2():
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.fill_between([0, 2.05], NT_FIG2_PLATEAU[0], NT_FIG2_PLATEAU[1], color="C3", alpha=0.12, zorder=0)
    ax.plot(NT_FIG2_TSTAR, NT_FIG2_VSTAR, "s--", color="C3", lw=1.6, ms=4, label="Nas & Tryggvason Fig 2 (digitized, ~64/D)")
    mfc = json.load(open(os.path.join(HERE, "results", "nt_fig2_mfc.json")))
    ax.plot(mfc["t_star"], mfc["V_star"], "o-", color="C0", lw=1.2, ms=4, label="MFC creeping Nx=32 (16/D)")
    ax.axhline(0.0, color="0.7", lw=0.8)
    ax.set_xlim(0, 2.05)
    ax.set_xlabel(r"$t^* = t/t_r$")
    ax.set_ylabel(r"$V^* = V/U_r$")
    ax.set_title("Fig 2: creeping drop ($Re=Ma=2.5\\times10^{-3}$, $Ca=10^{-3}$) — MFC vs Nas & Tryggvason")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5)
    ax.text(
        0.98,
        0.04,
        "MFC does NOT match: noisy mean $V^*\\approx0.4$ vs paper $\\approx0.13$\n(compressible/acoustic at the low $c_\\mathrm{ref}$ used to afford the creeping run)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.3",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9),
    )
    fig.tight_layout()
    out = os.path.join(FIGS, "nt_fig2_comparison.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def fig3():
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.plot(NT_FIG3_TSTAR, NT_FIG3_VSTAR, "s--", color="C3", lw=1.6, ms=4, label="Nas & Tryggvason Fig 3 (digitized, converged)")
    for name, lbl, col in [("w064", "MFC 32/D", "C0"), ("w128", "MFC 64/D", "C1")]:
        wd = os.path.join(HERE, "runs", "tc2", name)
        if not os.path.isfile(os.path.join(wd, "simulation.inp")):
            print(f"  (skip {name}: no run)")
            continue
        ts, vs = measure_vstar(wd)
        ax.plot(ts, vs, "o-", color=col, lw=1.0, ms=2.5, alpha=0.85, label=lbl)
    ax.axhline(0.0, color="0.7", lw=0.8)
    ax.set_xlim(left=0)
    ax.set_xlabel(r"$t^* = t/t_r$")
    ax.set_ylabel(r"$V^* = V/U_r$")
    ax.set_title("Fig 3: finite-$Re$ drop ($Re=5$, $Ma=20$, $Ca=0.0167$) — MFC vs Nas & Tryggvason")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8.5)
    ax.text(
        0.98,
        0.04,
        "MFC tracks the rise + overshoot (peak ~0.13-0.15) but over-declines\nafter the peak (compressible relaxation); shape matches, late tail does not",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.3",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9),
    )
    fig.tight_layout()
    out = os.path.join(FIGS, "nt_fig3_comparison.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig2()
    fig3()
