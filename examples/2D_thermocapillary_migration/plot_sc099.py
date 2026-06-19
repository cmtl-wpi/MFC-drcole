#!/usr/bin/env python3
"""Interface-thickness effect on the TC1/Fig 5 rise velocity: smooth_coeff 0.5 vs 0.99.

Overlays the conduction Ma-sweep (Ma = 0.1, 0.01, 0.001) at the default interface
(smooth_coeff = 0.5, half-width w = dx/0.5 = 2*dx) against the sharper interface
(smooth_coeff = 0.99, w ~ dx, ~2x sharper) rerun by run_sc099.py, plus Samareh's VOF
reference. Reuses plot.py's color-weighted rise-velocity extraction verbatim, so the
curves are read exactly as in case1_fig5.png.

Writes figures/case1_fig5_smooth_coeff.png.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot as P  # color_weighted_vy, v_ygb_ratio, SAMAREH_VOF, PLATE_STYLE, RUNS, FIGS

# (Ma label, baseline dir [sc=0.5], sharp dir [sc=0.99], color). Both interfaces are regenerated
# from the current committed case (same 2*t_r window, same dt) -- the only difference is smooth_coeff.
PAIRS = [
    (r"$Ma=0.1$",   "tc1/ma0p1/w064/sc050",   "tc1/ma0p1/w064/sc099",   "#9ecae1"),
    (r"$Ma=0.01$",  "tc1/ma0p01/w064/sc050",  "tc1/ma0p01/w064/sc099",  "#4292c6"),
    (r"$Ma=0.001$", "tc1/ma0p001/w064/sc050", "tc1/ma0p001/w064/sc099", "#084594"),
]


def curve(name):
    out = P.color_weighted_vy(os.path.join(P.RUNS, name))
    if out is None or len(out[0]) < 5:
        return None
    x, y = P.v_ygb_ratio(out)
    return x, y


def main():
    os.makedirs(P.FIGS, exist_ok=True)
    with plt.rc_context(P.PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
        ax.plot(P.SAMAREH_VOF[:, 0], P.SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0,
                mfc="none", mew=1.3, lw=1.0, label=r"Samareh Fig 5(d), VOF (digitized)")
        ax.axhline(1.0, color="0.3", lw=1.0, ls="--", zorder=1,
                   label=r"$u_{\mathrm{YGB}}$ (analytic terminal)")
        n = 0
        for ma_txt, base, sharp, color in PAIRS:
            cb, cs = curve(base), curve(sharp)
            if cb is not None:
                ax.plot(cb[0], cb[1], "o--", color=color, ms=3.5, mew=0, lw=0.9, alpha=0.55,
                        label=rf"{ma_txt}, $w=2\,\Delta x$ (smooth_coeff 0.5)")
                n += 1
            if cs is not None:
                ax.plot(cs[0], cs[1], "o-", color=color, ms=4.0, mew=0, lw=1.6,
                        label=rf"{ma_txt}, $w\approx\Delta x$ (smooth_coeff 0.99)")
                n += 1
        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 1.1)
        ax.set_xlabel(r"Time   $t/t_r$")
        ax.set_ylabel(r"Normalized Rise Velocity   $u/u_{\mathrm{YGB}}$")
        ax.set_title(r"Fig 5 — interface-thickness effect on rise velocity ($Ma\to0$)",
                     fontsize=12)
        ax.legend(loc="lower right", fontsize=8.5, frameon=False, ncol=1)
        dst = os.path.join(P.FIGS, "case1_fig5_smooth_coeff.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({n} curves)")


if __name__ == "__main__":
    main()
