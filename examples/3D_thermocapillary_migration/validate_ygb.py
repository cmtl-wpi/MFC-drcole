#!/usr/bin/env python3
"""Analyze the 3D thermocapillary u_YGB sweep: does v_t/v_YGB -> 1.0 as the three deficits vanish?

This does NOT run MFC -- run run_ygb.py first. It walks runs/ygb/, measures each leaf with
measure.py, parses (geom, W, Nx, Ma) back out of the leaf path, and produces three convergence
reductions, each holding the other two axes fixed at the converged corner (cube, W=10, Nx=80, Ma=0.5):

  confinement  ratio vs 1/W  -- the headline. Linear fit, extrapolate 1/W -> 0 (W -> infinity).
  grid         ratio vs dx   -- confirm the intercept is grid-independent.
  ma           ratio vs Ma   -- confirm it survives Ma -> 0 (perfectly invariant T).

Each writes a figure to figures/ and a reduction block to results/ygb_summary.json.

Usage:  python3 validate_ygb.py [confinement|grid|ma|all]   (default: all)
"""

import glob
import json
import os
import subprocess
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs", "ygb")
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "results", "ygb_summary.json")

# Converged corner: the two axes each reduction holds fixed while sweeping the third.
CORNER = {"W": 10.0, "Nx": 80, "Ma": 0.5}


def decode(token, prefix):
    """'w7p5' -> 7.5, 'nx080' -> 80.0, 'ma0p25' -> 0.25 (strip prefix, 'p' -> '.')."""
    return float(token[len(prefix) :].replace("p", "."))


def collect():
    """Measure every leaf under runs/ygb/ and tag it with (geom, W, Nx, Ma) from its path."""
    rows = []
    for inp in glob.glob(os.path.join(RUNS, "*", "*", "*", "*", "simulation.inp")):
        wd = os.path.dirname(inp)
        if not os.path.isdir(os.path.join(wd, "restart_data")):
            continue
        geom, wtok, nxtok, matok = os.path.relpath(wd, RUNS).split(os.sep)
        m = subprocess.run([sys.executable, os.path.join(HERE, "measure.py"), wd], capture_output=True, text=True, check=False)
        res = next((json.loads(l[len("RESULT_JSON ") :]) for l in m.stdout.splitlines() if l.startswith("RESULT_JSON ")), None)
        if res is None:
            print(f"  measure failed: {os.path.relpath(wd, HERE)}")
            continue
        rows.append({"geom": geom, "W": decode(wtok, "w"), "Nx": int(decode(nxtok, "nx")), "Ma": decode(matok, "ma"), **res})
    return rows


def save_block(key, entry):
    summ = json.load(open(SUMMARY)) if os.path.isfile(SUMMARY) else {}
    summ[key] = entry
    os.makedirs(os.path.dirname(SUMMARY), exist_ok=True)
    json.dump(summ, open(SUMMARY, "w"), indent=2)


def near(a, b):
    return abs(a - b) < 1e-6


def converge_confinement(rows):
    """ratio_plateau vs 1/W at fixed (cube, Nx, Ma); linear-extrapolate to 1/W -> 0 (unbounded)."""
    # The confinement line holds drop resolution fixed (cells_per_D = Nx/W ~ 8), so vary only W.
    pts = sorted((r["W"], r["ratio_plateau"]) for r in rows if r["geom"] == "cube" and near(r["Ma"], CORNER["Ma"]) and abs(r["Nx"] / r["W"] - 8.0) < 0.5)
    if len(pts) < 2:
        print(f"  confinement: need >=2 cube/cells_per_D~8/Ma{CORNER['Ma']} runs, have {len(pts)} -- skipping")
        return
    W = np.array([p[0] for p in pts])
    ratio = np.array([p[1] for p in pts])
    inv = 1.0 / W
    slope, intercept = np.polyfit(inv, ratio, 1)  # ratio ~ intercept + slope*(1/W); W->inf is the intercept
    print(f"  confinement: W={W.tolist()}  ratio={[round(x, 3) for x in ratio]}  -> W->inf intercept = {intercept:+.3f}")

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    xf = np.linspace(0.0, inv.max() * 1.05, 100)
    ax.plot(inv, ratio, "o", color="#084594", ms=7, label="MFC cube runs")
    ax.plot(xf, intercept + slope * xf, "-", color="#084594", lw=1.3, label=f"linear fit (W$\\to\\infty$: {intercept:.3f})")
    ax.plot(0.0, intercept, "*", color="k", ms=16, zorder=5, label=f"unbounded extrapolation = {intercept:.3f}")
    ax.axhline(1.0, ls="--", color="0.4", lw=1.1, label=r"$u_{\mathrm{YGB}}$ (analytic, ratio $=1$)")
    ax.set_xlabel(r"$1/W$  (inverse box width in $D$;  $0$ = unbounded)")
    ax.set_ylabel(r"plateau  $v_t / u_{\mathrm{YGB}}$")
    ax.set_xlim(left=-0.005)
    ax.set_title(r"Confinement convergence: $v_t/u_{\mathrm{YGB}} \to 1$ as $W\to\infty$")
    ax.legend(loc="best", fontsize=9, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "ygb_vs_confinement.png"), dpi=160)
    plt.close(fig)
    save_block("confinement", {"W": W.tolist(), "ratio_plateau": ratio.tolist(), "slope": float(slope), "ratio_Winf": float(intercept), "Nx": CORNER["Nx"], "Ma": CORNER["Ma"]})


def converge_grid(rows):
    """ratio_plateau vs dx at fixed (cube, W, Ma); dx = W/Nx."""
    pts = sorted((r["Nx"], r["ratio_plateau"]) for r in rows if r["geom"] == "cube" and near(r["W"], CORNER["W"]) and near(r["Ma"], CORNER["Ma"]))
    if len(pts) < 2:
        print(f"  grid: need >=2 cube/W{CORNER['W']}/Ma{CORNER['Ma']} runs, have {len(pts)} -- skipping")
        return
    Nx = np.array([p[0] for p in pts])
    ratio = np.array([p[1] for p in pts])
    dx = CORNER["W"] / Nx
    print(f"  grid: Nx={Nx.tolist()}  ratio={[round(x, 3) for x in ratio]}")

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(dx, ratio, "s-", color="#C44E52", ms=7, lw=1.2, label=f"MFC cube W={CORNER['W']:g}")
    ax.axhline(1.0, ls="--", color="0.4", lw=1.1, label=r"$u_{\mathrm{YGB}}$ (ratio $=1$)")
    ax.set_xlabel(r"cell size  $\Delta x = W/N_x$")
    ax.set_ylabel(r"plateau  $v_t / u_{\mathrm{YGB}}$")
    ax.set_title(r"Grid convergence (fixed box $W={}$, Ma$={}$)".format(int(CORNER["W"]), CORNER["Ma"]))
    ax.legend(loc="best", fontsize=9, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "ygb_vs_dx.png"), dpi=160)
    plt.close(fig)
    save_block("grid", {"Nx": Nx.tolist(), "dx": dx.tolist(), "ratio_plateau": ratio.tolist(), "W": CORNER["W"], "Ma": CORNER["Ma"]})


def converge_ma(rows):
    """ratio_plateau vs Ma at fixed (cube, W, Nx); extrapolate Ma -> 0 (perfectly invariant T)."""
    pts = sorted((r["Ma"], r["ratio_plateau"]) for r in rows if r["geom"] == "cube" and near(r["W"], CORNER["W"]) and r["Nx"] == CORNER["Nx"])
    if len(pts) < 2:
        print(f"  ma: need >=2 cube/W{CORNER['W']}/Nx{CORNER['Nx']} runs, have {len(pts)} -- skipping")
        return
    Ma = np.array([p[0] for p in pts])
    ratio = np.array([p[1] for p in pts])
    slope, intercept = np.polyfit(Ma, ratio, 1)  # ratio ~ intercept + slope*Ma; Ma->0 is the intercept
    print(f"  ma: Ma={Ma.tolist()}  ratio={[round(x, 3) for x in ratio]}  -> Ma->0 intercept = {intercept:+.3f}")

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    xf = np.linspace(0.0, Ma.max() * 1.05, 100)
    ax.plot(Ma, ratio, "^", color="#2ca25f", ms=8, label=f"MFC cube W={CORNER['W']:g}, Nx={CORNER['Nx']}")
    ax.plot(xf, intercept + slope * xf, "-", color="#2ca25f", lw=1.3, label=f"linear fit (Ma$\\to0$: {intercept:.3f})")
    ax.plot(0.0, intercept, "*", color="k", ms=16, zorder=5)
    ax.axhline(1.0, ls="--", color="0.4", lw=1.1, label=r"$u_{\mathrm{YGB}}$ (ratio $=1$)")
    ax.set_xlabel(r"Marangoni number  Ma  ($0$ = perfectly invariant $T$)")
    ax.set_ylabel(r"plateau  $v_t / u_{\mathrm{YGB}}$")
    ax.set_xlim(left=-0.02)
    ax.set_title(r"Finite-Ma convergence (fixed box $W={}$, $N_x={}$)".format(int(CORNER["W"]), CORNER["Nx"]))
    ax.legend(loc="best", fontsize=9, frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "ygb_vs_Ma.png"), dpi=160)
    plt.close(fig)
    save_block("ma", {"Ma": Ma.tolist(), "ratio_plateau": ratio.tolist(), "slope": float(slope), "ratio_Ma0": float(intercept), "W": CORNER["W"], "Nx": CORNER["Nx"]})


COMMANDS = {"confinement": converge_confinement, "grid": converge_grid, "ma": converge_ma}


def main():
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    if sel not in COMMANDS and sel != "all":
        sys.exit(f"usage: validate_ygb.py {{{'|'.join(COMMANDS)}|all}}")
    rows = collect()
    print(f"measured {len(rows)} runs under {os.path.relpath(RUNS, HERE)}/")
    # Record the raw per-run table too, so the summary stands alone.
    save_block(
        "runs", {os.path.join(r["geom"], f"w{r['W']:g}", f"nx{r['Nx']:03d}", f"ma{r['Ma']:g}"): {k: r[k] for k in ("ratio_plateau", "W", "Nx", "Ma", "cells_per_D", "t_end_tr", "rises")} for r in rows}
    )
    for name in COMMANDS if sel == "all" else [sel]:
        COMMANDS[name](rows)


if __name__ == "__main__":
    main()
