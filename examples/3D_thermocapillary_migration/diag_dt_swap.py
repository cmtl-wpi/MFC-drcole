#!/usr/bin/env python3
"""Diagnosis summary: the Ma-instability is set by the TIMESTEP, not the Marangoni number.

Plots v/v_YGB(t) for the Nx=32 controlled experiments. The dt-swap pair is the proof: at FIXED Ma
(fixed conduction strength k_therm), the big acoustic-limited dt is unstable while the small
diffusion-limited dt is stable -- in BOTH directions:
  Ma=0.1 : big dt (ma0p1) unstable  vs  small dt (ma0p1_dtlo) stable
  Ma=0.03: small dt (ma0p03) stable vs  big dt (ma0p03_dthi) -> ICFL blow-up
So dt is the lever; the apparent non-monotonicity in Ma is just that dt is acoustic-pinned for
Ma>=0.0625 then diffusion-limited (smaller) below.

Usage:  python3 diag_dt_swap.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "2D_thermocapillary_migration"))
from plot import PLATE_STYLE, color_weighted_vy, v_ygb_ratio  # noqa: E402

BASE = os.path.join(HERE, "runs", "diag_nx32")
# (label, legend, color, linestyle)
SERIES = [
    ("ma0p3",      "Ma=0.3,  acoustic dt  (stable ref)",            "#1b9e77", "-"),
    ("ma0p1",      "Ma=0.1,  big (acoustic) dt  -> UNSTABLE",       "#d62728", "-"),
    ("ma0p1_dtlo", "Ma=0.1,  small dt (= Ma=0.03's)  -> stable",    "#d62728", "--"),
    ("ma0p05",     "Ma=0.05, dt at crossover  -> UNSTABLE",         "#ff7f0e", "-"),
    ("ma0p03",     "Ma=0.03, small (diffusion) dt  -> stable",      "#1f77b4", "-"),
]


def main():
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
        ax.axhline(1.3, color="0.6", ls=":", lw=1.0)
        ax.text(0.02, 1.32, "unstable threshold (max|v/v_YGB|>1.3)", color="0.5", fontsize=8)
        for label, leg, c, ls in SERIES:
            out = color_weighted_vy(os.path.join(BASE, label))
            if out is None or len(out[0]) < 3:
                print(f"  {label}: few snaps, skip"); continue
            x, y = v_ygb_ratio(out)
            ax.plot(x, np.clip(y, -0.3, 2.5), ls, color=c, lw=1.8, alpha=0.9, label=leg)
        # ma0p03_dthi blew up (ICFL) at step 13 -- mark it
        ax.scatter([13 * 2.668e-3 / 7.5], [2.4], marker="X", s=120, color="k", zorder=6)
        ax.annotate("Ma=0.03 at big dt:\nICFL blow-up, step 13", xy=(13 * 2.668e-3 / 7.5, 2.4),
                    xytext=(0.25, 2.15), fontsize=8.5, color="k",
                    arrowprops=dict(arrowstyle="->", color="0.4"))
        ax.set_xlim(0, 0.65); ax.set_ylim(-0.3, 2.6)
        ax.set_xlabel(r"$t / t_r$")
        ax.set_ylabel(r"rise velocity   $v / v_{\mathrm{YGB}}$")
        ax.set_title("3D thermocapillary @ Nx=32: dt sets stability, not Ma\n"
                     "(at fixed Ma, big dt unstable / small dt stable -- both directions)",
                     fontsize=11, loc="left")
        ax.legend(loc="upper right", fontsize=8.5, frameon=False)
        dst = os.path.join(HERE, "figures", "diag_dt_swap.png")
        fig.savefig(dst, dpi=200); plt.close(fig)
        print(f"  wrote {dst}")


if __name__ == "__main__":
    main()
