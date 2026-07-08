#!/usr/bin/env python3
"""
Phase 1 seam analysis for one static-drop run (uniform or AMR).

Measures, at each saved time:
  * domain-wide max|u|            (parasitic-current level)
  * seam-band max|u|  (AMR only)  (spurious coarse/fine seam current)
  * Laplace pressure jump vs sigma/R
  * drop-center drift
  * containment audit: how far cf / alpha depart from {0,1} in the seam band

Raw |u| fields are rendered (t=0 and final) BEFORE any scalar is reported
(RESEARCH_WORKFLOWS s7). Emits one RESULT_JSON line and writes results/<label>.json.

Usage: seam_analysis.py <run_dir> [--label NAME]
"""
import argparse
import json
import math
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mfc_read as mr

EXP = "/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st"

ap = argparse.ArgumentParser()
ap.add_argument("run_dir")
ap.add_argument("--label", default=None)
ap.add_argument("--band-block", default=None,
                help="bx0,bx1,by0,by1 cell indices to force a seam-band mask on a "
                     "non-AMR run (control: coarse max|u| in the AMR seam region)")
args = ap.parse_args()

run_dir = os.path.abspath(args.run_dir)
label = args.label or os.path.basename(run_dir.rstrip("/"))
figdir = os.path.join(EXP, "results", "figures")
os.makedirs(figdir, exist_ok=True)

inp = mr.read_inp(run_dir)                     # grid/EOS/algorithm/AMR
pp = mr.read_inp(run_dir, "pre_process.inp")   # patch geometry/density
x, y = mr.grid(run_dir)
steps = mr.list_steps(run_dir)
t_save = inp["t_save"]
X, Y = np.meshgrid(x, y)                       # (ny, nx)

# physics scales
sigma = inp["sigma"]
R = pp["patch_icpp(2)%radius"]
rho_l = pp["patch_icpp(2)%alpha_rho(1)"] / (1.0 - 1e-9)       # (1-eps)*rho_l / (1-eps)
U_sigma = math.sqrt(sigma / (rho_l * R))       # capillary velocity scale
dP_laplace = sigma / R

# geometry: interface radius from origin
Rc = np.sqrt(X**2 + Y**2)

# seam band mask (AMR run, or forced on a control run via --band-block)
blk = mr.amr_block(inp)
if args.band_block is not None:
    blk = tuple(int(v) for v in args.band_block.split(","))
w = 2 * mr.buff_size(inp)
seam_mask = None
if blk is not None:
    bx0, bx1, by0, by1 = blk
    ix = np.arange(len(x))[None, :] * np.ones((len(y), 1))
    iy = np.arange(len(y))[:, None] * np.ones((1, len(x)))
    near_x = (np.abs(ix - bx0) <= w) | (np.abs(ix - bx1) <= w)
    near_y = (np.abs(iy - by0) <= w) | (np.abs(iy - by1) <= w)
    inside_ext = (ix >= bx0 - w) & (ix <= bx1 + w) & (iy >= by0 - w) & (iy <= by1 + w)
    seam_mask = inside_ext & (near_x | near_y)   # a frame of width w around the block edge


def scalars(fld):
    speed = fld["speed"]
    a1, cf, pres = fld["alpha1"], fld["cf"], fld["pres"]
    drop = a1 > 0.9
    far = (a1 < 0.1) & (Rc > 2.0 * R)
    p_jump = (pres[drop].mean() - pres[far].mean()) if drop.any() and far.any() else np.nan
    wgt = np.clip(a1, 0.0, 1.0)
    xc = (wgt * X).sum() / wgt.sum()
    yc = (wgt * Y).sum() / wgt.sum()
    drift = math.hypot(xc, yc) / R
    out = dict(
        max_speed=float(speed.max()),
        p_jump=float(p_jump),
        drift_over_R=float(drift),
    )
    if seam_mask is not None:
        out["max_speed_seam"] = float(speed[seam_mask].max())
        out["cf_dist_seam"] = float(np.minimum(cf[seam_mask], 1.0 - cf[seam_mask]).max())
        out["alpha_dist_seam"] = float(np.minimum(a1[seam_mask], 1.0 - a1[seam_mask]).max())
    return out


# ---- time series --------------------------------------------------------------
series = {k: [] for k in
          ("t", "max_speed", "max_speed_seam", "p_jump", "drift_over_R",
           "cf_dist_seam", "alpha_dist_seam")}
first_fld = last_fld = None
for s in steps:
    fld = mr.read_step(run_dir, s, inp, x, y)
    sc = scalars(fld)
    series["t"].append(s * t_save)
    for k in ("max_speed", "p_jump", "drift_over_R"):
        series[k].append(sc[k])
    for k in ("max_speed_seam", "cf_dist_seam", "alpha_dist_seam"):
        series[k].append(sc.get(k, np.nan))
    if s == steps[0]:
        first_fld = fld
    last_fld = fld

t = np.array(series["t"])
tau = math.sqrt(rho_l * R**3 / sigma)


# ---- RAW FIELD FIRST: |u| pcolormesh at t0 and final --------------------------
def plot_speed(fld, tag, tval):
    fig, ax = plt.subplots(figsize=(6, 5.2))
    sp = fld["speed"] / U_sigma
    pcm = ax.pcolormesh(X / R, Y / R, sp, shading="auto", cmap="inferno")
    ax.contour(X / R, Y / R, fld["alpha1"], levels=[0.5], colors="cyan", linewidths=0.8)
    if blk is not None:
        ax.add_patch(Rectangle((x[bx0] / R, y[by0] / R),
                               (x[bx1] - x[bx0]) / R, (y[by1] - y[by0]) / R,
                               fill=False, ec="lime", lw=1.2, ls="--"))
    ax.set_aspect("equal")
    ax.set_xlabel("x / R"); ax.set_ylabel("y / R")
    ax.set_title(f"{label}: |u|/U_sigma  {tag} (t={tval:.2e}s = {tval/tau:.2f} tau)")
    fig.colorbar(pcm, ax=ax, label="|u| / U_sigma", fraction=0.046)
    fp = os.path.join(figdir, f"{label}_speed_{tag}.png")
    fig.savefig(fp, dpi=150, bbox_inches="tight"); plt.close(fig)
    return fp


f_t0 = plot_speed(first_fld, "t0", t[0])
f_tf = plot_speed(last_fld, "final", t[-1])

# time-series of parasitic current
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t / tau, np.array(series["max_speed"]) / U_sigma, "-o", ms=3, label="domain max|u|")
if seam_mask is not None:
    ax.plot(t / tau, np.array(series["max_speed_seam"]) / U_sigma, "-s", ms=3,
            color="crimson", label="seam-band max|u|")
ax.set_xlabel("t / tau"); ax.set_ylabel("max|u| / U_sigma")
ax.set_title(f"{label}: parasitic current vs time")
ax.legend(); ax.grid(alpha=0.3)
f_ts = os.path.join(figdir, f"{label}_maxu_vs_t.png")
fig.savefig(f_ts, dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- summary ------------------------------------------------------------------
ms = np.array(series["max_speed"])
half = len(ms) // 2


def growth(a):
    """2nd-half mean / 1st-half mean; >~1.5 => growing. None if 1st half ~0."""
    if half == 0:
        return None
    d = a[:half].mean()
    return float(a[half:].mean() / d) if abs(d) > 1e-30 else None
summary = dict(
    label=label,
    run_dir=run_dir,
    amr=blk is not None,
    n_steps=len(steps),
    t_final=float(t[-1]),
    t_final_over_tau=float(t[-1] / tau),
    U_sigma=U_sigma,
    dP_laplace=dP_laplace,
    tau=tau,
    seam_band_w=int(w) if blk else None,
    # parasitic current
    max_speed_final=float(ms[-1]),
    max_speed_peak=float(ms.max()),
    max_speed_final_over_Usigma=float(ms[-1] / U_sigma),
    # growth: ratio of 2nd-half mean to 1st-half mean (>~1.5 => growing)
    growth_ratio=growth(ms),
    # pressure jump / drift at final
    p_jump_final=float(series["p_jump"][-1]),
    p_jump_rel_err=float((series["p_jump"][-1] - dP_laplace) / dP_laplace),
    drift_final_over_R=float(series["drift_over_R"][-1]),
)
if seam_mask is not None:
    sm = np.array(series["max_speed_seam"])
    summary.update(
        max_speed_seam_final=float(sm[-1]),
        max_speed_seam_peak=float(sm.max()),
        seam_growth_ratio=growth(sm),
        containment_cf_max=float(np.nanmax(series["cf_dist_seam"])),
        containment_alpha_max=float(np.nanmax(series["alpha_dist_seam"])),
        containment_pass=bool(np.nanmax(series["cf_dist_seam"]) < 1e-9),
    )

os.makedirs(os.path.join(EXP, "results"), exist_ok=True)
with open(os.path.join(EXP, "results", f"{label}.json"), "w") as f:
    json.dump({"summary": summary, "series": series}, f, indent=2)

print(f"[seam_analysis] figures: {f_t0}, {f_tf}, {f_ts}", file=sys.stderr)
print("RESULT_JSON " + json.dumps(summary))
