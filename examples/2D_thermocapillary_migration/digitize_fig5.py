#!/usr/bin/env python3
"""Overlay the four numerical methods of Samareh, Mostaghimi & Moreau (2014) Fig 5 on a single
axes -- the published figure shows them in four separate panels, so they cannot be compared
directly. Each method's converged (n_x = 256, solid) "normalized rise velocity vs time" trace is
digitized BY EYE from the published raster (accuracy ~ +/- 0.02 in v/v_YGB), then MFC's own n_x = 256
curve is overlaid for a direct head-to-head against the paper's VOF method (MFC is a compressible
finite-volume VOF-type solver, so panel (d) is the closest comparison).

The four methods (Samareh Fig 5 a-d):
  (a) conservative level-set        -- overshoots, settles ~0.93-0.97
  (b) combined level-set / VOF      -- plateau ~0.80-0.86, oscillatory
  (c) refined level-set on a grid   -- cleanest, flat plateau ~0.82
  (d) volume of fluid (VOF)         -- plateau ~0.82-0.84   <- MFC's closest analogue

Samareh's Fig 5 x-axis is labeled "Time"; we take it to be the capillary-thermal time t/t_r
(t_r = mu/|sigma_T gradT|, the paper's declared time scale -- the caption does not state the unit, so
this is an assumption) and plot MFC on the same normalization. MFC's 64-grid run extends to t/t_r ~ 6:
it matches the VOF/RLSG plateau (~0.80) in the quasi-steady window t/t_r ~ [0.5, 3], then drifts DOWN as
the frozen (Ma=0) density-proxy T is advected and distorts -- Samareh holds T invariant (infinite
diffusivity), so their curves stay flat to 10.

Usage:  python3 digitize_fig5.py
Writes: figures/case1_fig5_samareh_methods_overlay.png
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
FIGS = os.path.join(HERE, "figures")
GRADT = 2.0 / 15.0  # imposed |dT/dy|
R = 0.5  # drop radius (D = 1)

sns.set_theme(style="whitegrid", font_scale=1.15)

# -- Digitized control points (t/t_r, v/v_YGB) of each method's n_x = 256 (solid) curve in Fig 5.
#    Read by eye from the published raster; dense enough that linear interpolation looks smooth. --
SAMAREH = {
    "(a) conservative level-set": dict(
        color="#000000",
        ls="-",
        pts=[
            (0.0, 0.0),
            (0.15, 0.45),
            (0.3, 0.70),
            (0.5, 0.80),
            (0.8, 0.83),
            (1.2, 0.85),
            (1.8, 0.89),
            (2.4, 0.94),
            (2.8, 0.97),
            (3.4, 0.92),
            (4.2, 0.905),
            (5.0, 0.95),
            (5.8, 0.915),
            (6.5, 0.92),
            (7.4, 0.955),
            (8.2, 0.945),
            (9.0, 0.935),
            (10.0, 0.975),
        ],
    ),
    "(b) combined level-set / VOF": dict(
        color="#009E73",
        ls="-",
        pts=[
            (0.0, 0.0),
            (0.15, 0.45),
            (0.3, 0.72),
            (0.5, 0.80),
            (0.8, 0.81),
            (1.2, 0.83),
            (1.8, 0.86),
            (2.2, 0.865),
            (2.8, 0.83),
            (3.5, 0.80),
            (4.2, 0.82),
            (5.0, 0.85),
            (5.8, 0.83),
            (6.5, 0.80),
            (7.2, 0.82),
            (8.0, 0.845),
            (8.8, 0.83),
            (10.0, 0.835),
        ],
    ),
    "(c) refined level-set (RLSG)": dict(
        color="#56B4E9", ls="-", pts=[(0.0, 0.0), (0.1, 0.40), (0.25, 0.70), (0.4, 0.79), (0.6, 0.81), (1.0, 0.815), (2.0, 0.818), (4.0, 0.82), (6.0, 0.82), (8.0, 0.821), (10.0, 0.822)]
    ),
    "(d) volume of fluid (VOF)": dict(
        color="#D55E00",
        ls="-",
        pts=[
            (0.0, 0.0),
            (0.12, 0.42),
            (0.28, 0.70),
            (0.45, 0.80),
            (0.7, 0.815),
            (1.0, 0.82),
            (2.0, 0.825),
            (3.0, 0.83),
            (4.0, 0.83),
            (5.0, 0.835),
            (6.0, 0.83),
            (7.0, 0.835),
            (8.0, 0.838),
            (9.0, 0.835),
            (10.0, 0.84),
        ],
    ),
}


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def running_mean(y, w):
    y = np.asarray(y, float)
    half = w // 2
    return np.array([y[max(0, i - half) : min(len(y), i + half + 1)].mean() for i in range(len(y))])


def mfc_curve(run_dir):
    """MFC color-weighted lab-frame rise velocity -> (t/t_r, v/v_YGB). None if the run is absent."""
    inp = os.path.join(run_dir, "simulation.inp")
    rd = os.path.join(run_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    P = read_namelist(inp)
    f = lambda k: float(P[k.lower()])  # noqa: E731
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    cells = nx * ny * nz
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    ts = str(P.get("thermal_scalar", "F")).strip(". ").upper().startswith("T")
    c_idx = nvars - 2 if ts else nvars - 1
    mu = 1.0 / f("fluid_pp(1)%re(1)")
    dsdt = f("sigma_dtdt")
    t_r = mu / abs(dsdt * GRADT)
    v_YGB = (2.0 / 15.0) * (-dsdt) * GRADT * R / mu
    t, u = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        vy = fld(3) / (fld(0) + fld(1))
        c = np.clip(fld(c_idx), 0.0, None)
        t.append(s * f("dt") / t_r)
        u.append((c * vy).sum() / c.sum() / v_YGB)
    return np.array(t), np.array(u)


fig, ax = plt.subplots(figsize=(9.2, 5.6), constrained_layout=True)
tt = np.linspace(0, 10, 400)
for label, d in SAMAREH.items():
    px, py = zip(*d["pts"])
    ax.plot(tt, np.interp(tt, px, py), d["ls"], color=d["color"], lw=2.4, label=f"Samareh {label}")


# MFC overlay (closest analogue is the VOF method, panel d). Use the FINEST full-duration (tr=10) run
# available so the trace spans the same t/t_r range as the published curves; fall back to the short
# tr=2 256-grid run only if no full-duration run has finished yet.
def pick_mfc_run():
    cands = []
    for dd in glob.glob(os.path.join(HERE, "runs", "fig5_2D_w*_tr10")):
        mm = re.search(r"fig5_2D_w(\d+)_tr10", dd)
        if mm and glob.glob(os.path.join(dd, "restart_data", "lustre_[0-9]*.dat")):
            cands.append((int(mm.group(1)), dd))
    return max(cands) if cands else (256, os.path.join(HERE, "runs", "fig5_2D_w256"))


mfc_nx, mfc_dir = pick_mfc_run()
mfc = mfc_curve(mfc_dir)
if mfc is not None:
    t, u = mfc
    w = max(5, len(t) // 14)
    ax.plot(t, u, ".", color="#CC79A7", ms=3, alpha=0.12, zorder=3)  # faint raw cloud (the acoustic ring)
    ax.plot(t, running_mean(u, w), "-", color="#CC79A7", lw=3.0, zorder=6, label=f"MFC (this work, VOF, $n_x={mfc_nx}$)")

ax.axhline(0.80, ls=":", color="0.5", lw=1.3)
ax.text(9.9, 0.805, "0.80", ha="right", va="bottom", fontsize=10, color="0.4")
ax.set_xlabel(r"Time  $t/t_r$  ($t_r = \mu/|\sigma_T \nabla T|$)")
ax.set_ylabel(r"Normalized rise velocity  $v/v_{\mathrm{YGB}}$")
ax.set_title("Samareh (2014) Fig 5 — four methods digitized & overlaid (converged $n_x=256$), vs MFC")
ax.set_xlim(0, 10)
ax.set_ylim(0, 1.05)
ax.legend(loc="lower right", fontsize=10, frameon=True, framealpha=0.93)
sns.despine(ax=ax)
os.makedirs(FIGS, exist_ok=True)
out = os.path.join(FIGS, "case1_fig5_samareh_methods_overlay.png")
fig.savefig(out, dpi=200)
print(f"wrote {out}")
