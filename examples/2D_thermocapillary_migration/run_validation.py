#!/usr/bin/env python3
"""Reproduce Samareh-2014 Fig 5 (2D zero-Marangoni thermocapillary RISE velocity) with MFC.

A drop starts 1.5D above the cold floor of a slip-walled box and rises toward the hot top in the
Ma=0 (frozen linear-T) limit; Fig 5 is the grid-convergence of the normalized rise velocity
v/v_YGB at Samareh's box-width resolutions 64/128/256 -> quasi-steady plateau ~ 0.80.

Each variant runs a COPY of case.py from its own runs/<name>/ directory (MFC writes output next to
the case file) via `./mfc.sh run`, then this script reads the restart data directly and overlays the
v/v_YGB(t/t_r) curves -- mirroring measure.py's rise-velocity protocol (lab-frame in the slip-wall
box). Results are aggregated to results/summary.json and the comparison figure is written to
results/. CAVEAT (no bulk conduction, frozen-T IC): the curve ramps, overshoots, settles to a
quasi-steady plateau, then slowly drifts -- the plateau (not the endpoint) is the Samareh comparison.

This example is 2D ONLY. Samareh's 3D companion (Fig 6) and the finite-Ma figures (7/8/12/13/16) are
NOT here: Fig 6 belongs in a separate 3D example and is not validatable on this no-conduction branch
(the 3D rise drifts unboundedly, no plateau); 7/8/12/13/16 need bulk conduction + mu(T) MFC lacks.

Usage (invokes mpirun, so run from a normal shell):
    python3 run_validation.py [run|remeasure]   (default: run)

`remeasure` rebuilds results/summary.json and the figure from existing runs/<name>/ directories
WITHOUT running any simulation -- use it after editing the measurement or plotting code.
"""

import glob
import json
import os
import re
import shutil
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")
RESULTS = os.path.join(HERE, "results")

SAMAREH_2D = 0.80  # Samareh Fig 5: converged 2D ratio in their slip-wall box

# 2D grid convergence (Samareh slip-wall box). Rank counts respect MFC's decomposition rule: each
# local block must hold >= num_stcls_min * weno_order = 5*5 = 25 cells per split dimension. tr = run
# length in capillary-thermal times t_r; long enough to show ramp + overshoot + the quasi-steady
# plateau (the Samareh comparison) + the onset of the frozen-T drift.
# ranks must factor into blocks of >= 25 cells per split dim on the nx x (1.5*nx) grid:
#   w064 (64x96):   6 = 2x3   -> 32x32
#   w128 (128x192): 16 = 4x4  -> 32x48
#   w256 (256x384): 64 = 8x8  -> 32x48
# Frozen-T grid convergence (Samareh Fig 5) plus the finite-Ma modes at w128: bulk conduction
# (SAMAREH_MA, with the wall-only isothermal BC fix) and the independent temperature scalar
# (SAMAREH_TS, density-decoupled). The mode variants share the slip-wall geometry so all curves are
# lab-frame-comparable. `env` overrides the per-variant case.py knobs.
VARIANTS = [
    dict(name="fig5_2D_w064", nx=64, ranks=6, tr=2.0, env={}),
    dict(name="fig5_2D_w128", nx=128, ranks=16, tr=2.0, env={}),
    dict(name="fig5_2D_w256", nx=256, ranks=64, tr=2.0, env={}),
    dict(name="modes_conduction_w128", nx=128, ranks=16, tr=2.0, env={"SAMAREH_MA": "0.3"}),
    dict(name="modes_thermal_scalar_w128", nx=128, ranks=16, tr=2.0, env={"SAMAREH_TS": "1", "SAMAREH_MA": "0.3"}),
]

# Curves overlaid in the finite-Ma modes-comparison figure (reuses the frozen-T w128 run as baseline).
MODE_CURVES = [
    ("fig5_2D_w128", "C7", "frozen-T (Ma=0)"),
    ("modes_conduction_w128", "C0", "conduction (Ma=0.3, wall isothermal BC)"),
    ("modes_thermal_scalar_w128", "C3", "thermal_scalar (Ma=0.3, T decoupled)"),
]


def run_variant(v):
    """Run one MFC case (2D, Samareh slip-wall geometry). Returns True on success."""
    workdir = os.path.join(RUNS, v["name"])
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    shutil.copy(os.path.join(HERE, "case.py"), os.path.join(workdir, "case.py"))

    env = {**os.environ, "SAMAREH_NX": str(v["nx"]), "SAMAREH_WALL": "1", "SAMAREH_TR": str(v["tr"]), **v.get("env", {})}
    rel_case = os.path.relpath(os.path.join(workdir, "case.py"), REPO)
    print(f"\n>>> {v['name']}: width={v['nx']} ranks={v['ranks']} tr={v['tr']}", flush=True)
    proc = subprocess.run(["./mfc.sh", "run", rel_case, "-n", str(v["ranks"])], cwd=REPO, env=env, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(f"  RUN FAILED (exit {proc.returncode}). Last stderr/stdout:")
        print("\n".join((proc.stdout + proc.stderr).splitlines()[-25:]))
        return False
    return True


def measure(workdir):
    """Run measure.py on a case directory and return its RESULT_JSON dict (or None)."""
    meas = subprocess.run([sys.executable, os.path.join(HERE, "measure.py"), workdir], capture_output=True, text=True, check=False)
    print("\n".join(meas.stdout.splitlines()[-5:]))
    if meas.returncode != 0:
        print(f"  MEASURE FAILED: {meas.stderr[-500:]}")
        return None
    for line in meas.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


def ratio_curve(case_dir):
    """Full v(t)/v_YGB rise history of a run, read straight from restart data. Returns
    (t/t_r, ratio) or None. Mirrors measure.py's rise protocol (lab-frame in the slip-wall box)."""
    inp = os.path.join(case_dir, "simulation.inp")
    rd = os.path.join(case_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    p = {}
    with open(inp) as fh:
        for line in fh:
            if "=" in line:
                k, val = line.split("=", 1)
                p[k.strip().lower()] = val.strip().rstrip(",")
    f = lambda k: float(p[k.lower()])  # noqa: E731
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    dt, mu, dsdt = f("dt"), 1.0 / f("fluid_pp(1)%re(1)"), f("sigma_dtdt")
    wall = int(f("bc_y%beg")) == -2
    ts = str(p.get("thermal_scalar", "F")).upper().strip(". ").startswith("T")  # T_s appended after color c
    gradT = 2.0 / 15.0
    t_r = mu / abs(dsdt * gradT)
    v_YGB = (2.0 / 15.0) * (-dsdt) * gradT * 0.5 / mu
    cells = nx * ny * nz
    Ly = f("y_domain%end") - f("y_domain%beg")
    far = 0.75 * Ly / 2.0
    yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
    y = 0.5 * (yb[:-1] + yb[1:])
    is_far = np.abs(y) > far
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    t_tr, ratio = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        vy = fld(3) / (fld(0) + fld(1))
        c = np.clip(fld(nvars - 2 if ts else nvars - 1), 0.0, None)  # color c; T_s is last in ts mode
        u = (c * vy).sum() / c.sum() - (0.0 if wall else vy[:, is_far, :].mean())
        t_tr.append(s * dt / t_r)
        ratio.append(u / v_YGB)
    return np.array(t_tr), np.array(ratio)


def make_figure():
    """Overlay the v/v_YGB(t/t_r) rise curves for every grid against the Samareh 2D plateau."""
    colors = ["C0", "C1", "C2"]
    curves = []
    for i, v in enumerate(VARIANTS):
        cur = ratio_curve(os.path.join(RUNS, v["name"]))
        if cur is not None:
            curves.append((f"{v['nx']} cells/box width ({v['nx'] / 5:.0f}/$D$)", colors[i % len(colors)], cur))
    if not curves:
        print("  (no runs found; nothing to plot)")
        return
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for lab, col, (t, rr) in curves:
        ax.plot(t, rr, "-", color=col, lw=1.7, label="MFC " + lab)
    ax.axhline(1.0, ls=":", color="0.45", lw=1.2, label=r"$v_{\mathrm{YGB}}$ (Samareh Eq. 29)")
    ax.axhline(SAMAREH_2D, ls="--", color="0.2", lw=1.4, label=r"Samareh 2D $\approx$ 0.80 (Fig 5)")
    ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$, Samareh time scale)")
    ax.set_ylabel(r"normalized rise velocity  $v / v_{\mathrm{YGB}}$")
    ax.set_title("Fig 5: 2D thermocapillary rise, grid convergence (Samareh slip-wall box)")
    ax.set_xlim(left=0.0)
    ax.set_ylim(0.0, 1.15)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, "fig5_rise_velocity_2D.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def make_modes_figure():
    """Overlay frozen-T vs conduction vs thermal_scalar rise curves at w128 (the finite-Ma modes).
    All in the slip-wall box (lab frame). y-axis is left unclipped so any reversal would be visible."""
    curves = []
    for name, col, lab in MODE_CURVES:
        cur = ratio_curve(os.path.join(RUNS, name))
        if cur is not None:
            curves.append((lab, col, cur))
    if not curves:
        return
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for lab, col, (t, rr) in curves:
        ax.plot(t, rr, "-", color=col, lw=1.7, label="MFC " + lab)
    ax.axhline(1.0, ls=":", color="0.45", lw=1.2, label=r"$v_{\mathrm{YGB}}$ (Samareh Eq. 29)")
    ax.axhline(SAMAREH_2D, ls="--", color="0.2", lw=1.4, label=r"Samareh 2D $\approx$ 0.80")
    ax.axhline(0.0, ls="-", color="0.6", lw=0.8)
    ax.set_xlabel(r"$t / t_r$  ($t_r = \mu / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"normalized rise velocity  $v / v_{\mathrm{YGB}}$")
    ax.set_title("Finite-Ma modes at 128 cells/width (slip-wall box): forward migration restored")
    ax.set_xlim(left=0.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, "finite_ma_modes_2D.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def print_table(summary):
    # Only rows produced by the current measure.py have these keys; skip stale entries from older
    # example versions (e.g. a pre-merge summary.json) rather than crashing on a missing key.
    cols = ("cells_per_D", "rho_drop", "ratio_plateau", "t_plateau_tr", "overshoot", "ratio_final", "samareh_ratio")
    rows = {k: v for k, v in summary.items() if all(c in v for c in cols)}
    skipped = [k for k in summary if k not in rows]
    print("\n" + "=" * 84)
    print(f"{'variant':>16} {'cells/D':>8} {'rho_drop':>8} {'plateau':>8} {'@t/t_r':>7} {'oversh':>7} {'endpt':>7} {'ref':>5}")
    for name, r in sorted(rows.items(), key=lambda kv: kv[1]["cells_per_D"]):
        print(
            f"{name:>16} {r['cells_per_D']:>8.1f} {r['rho_drop']:>8.3f} {r['ratio_plateau']:>+8.3f} "
            f"{r['t_plateau_tr']:>7.2f} {r['overshoot']:>+7.3f} {r['ratio_final']:>+7.3f} {r['samareh_ratio']:>5.2f}"
        )
    print("=" * 84)
    if skipped:
        print(f"(skipped {len(skipped)} stale entr{'y' if len(skipped) == 1 else 'ies'} missing current fields: {', '.join(sorted(skipped))})")


def load_summary():
    path = os.path.join(RESULTS, "summary.json")
    return json.load(open(path)) if os.path.isfile(path) else {}


def save_summary(summary):
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(summary, open(os.path.join(RESULTS, "summary.json"), "w"), indent=2)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"

    if mode == "remeasure":
        summary = {}
        for name in sorted(os.listdir(RUNS)) if os.path.isdir(RUNS) else []:
            workdir = os.path.join(RUNS, name)
            if not os.path.isfile(os.path.join(workdir, "simulation.inp")):
                continue
            print(f"\n>>> remeasure {name}", flush=True)
            res = measure(workdir)
            if res is not None:
                summary[name] = res
        save_summary(summary)
    else:
        summary = load_summary()
        for v in VARIANTS:
            if run_variant(v):
                res = measure(os.path.join(RUNS, v["name"]))
                if res is not None:
                    summary[v["name"]] = res
                    save_summary(summary)  # checkpoint after every successful variant
                    make_figure()  # rebuild the figures incrementally
                    make_modes_figure()

    make_figure()
    make_modes_figure()
    print_table(summary)
    print(f"\nwrote {os.path.join(RESULTS, 'summary.json')} and figure in {RESULTS}/")
