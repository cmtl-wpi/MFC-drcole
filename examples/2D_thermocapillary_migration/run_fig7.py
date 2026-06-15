#!/usr/bin/env python3
"""Reproduce Samareh-2014 Fig 7 (2D finite-Marangoni thermocapillary migration) with MFC.

Runs case_fig7.py (the Nas & Tryggvason Re=5/Ma=20/Ca=0.01666 test: a real two-fluid drop of
property ratio 0.5 in a closed slip-wall box, with bulk conduction of an independent temperature
scalar) at Samareh's two box-width resolutions (n_x = 64, 128), measures U*(t*) for each via
measure_fig7.py, and overlays them against the Nas & Tryggvason overshoot peak (~0.13).

The walls use the faithful isothermal Dirichlet BC (the case default), which sinks the drop's thermal
wake so the migration reaches a clean plateau (peak ~0.14, terminal ~0.10, matching Samareh). This
requires the MPI rank-guard fix in s_apply_thermal_conduction_bc: before it, the isothermal BC
overwrote interior ranks' halo cells, scrambling T into rank-boundary bands and reversing the drop
(the CONDUCTION_REVERSAL_SAGA root cause -- not a density proxy or advective throughflow). FIG7_ADIABATIC=1
runs adiabatic walls instead (forward, but the un-sunk wake makes the plateau over-decline). See
case_fig7.py and CONDUCTION_REVERSAL_SAGA.md.

Unlike run_validation.py (Fig 5, the Ma=0 frozen-T case normalized by v_YGB), this driver normalizes
by the Marangoni scales U_r, t_r and compares to the published finite-Ma curve -- the validation of
MFC's heat-transfer (thermal_conduction + thermal_scalar) module.

Usage (invokes mpirun, so run from a normal shell):
    python3 run_fig7.py [run|remeasure]   (default: run)
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

NAS_TRYGGVASON_PEAK = 0.13  # Samareh Fig 7 overshoot peak of U* = U/U_r

# Samareh's two grids. Rank counts respect MFC's decomposition rule (>= 5*5 = 25 cells per split dim
# on the n_x x 2*n_x grid):  w064 (64x128): 8 = 2x4 -> 32x32;  w128 (128x256): 16 = 4x4 -> 32x64.
VARIANTS = [
    dict(name="fig7_w064", nx=64, ranks=8, color="C0"),
    dict(name="fig7_w128", nx=128, ranks=16, color="C3"),
]


def run_variant(v):
    """Run one case_fig7.py variant in its own runs/<name>/ directory. Returns True on success."""
    workdir = os.path.join(RUNS, v["name"])
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    shutil.copy(os.path.join(HERE, "case_fig7.py"), os.path.join(workdir, "case_fig7.py"))

    env = {**os.environ, "FIG7_NX": str(v["nx"])}
    rel_case = os.path.relpath(os.path.join(workdir, "case_fig7.py"), REPO)
    print(f"\n>>> {v['name']}: width={v['nx']} ranks={v['ranks']}", flush=True)
    proc = subprocess.run(["./mfc.sh", "run", rel_case, "-n", str(v["ranks"])], cwd=REPO, env=env, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(f"  RUN FAILED (exit {proc.returncode}). Last output:")
        print("\n".join((proc.stdout + proc.stderr).splitlines()[-25:]))
        return False
    return True


def measure(workdir):
    """Run measure_fig7.py on a case directory; return its RESULT_JSON dict (or None)."""
    meas = subprocess.run([sys.executable, os.path.join(HERE, "measure_fig7.py"), workdir], capture_output=True, text=True, check=False)
    print("\n".join(meas.stdout.splitlines()[-4:]))
    if meas.returncode != 0:
        print(f"  MEASURE FAILED: {meas.stderr[-500:]}")
        return None
    for line in meas.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


def curve(case_dir):
    """Lightweight U*(t*) history of a run, read straight from restart data, for the overlay figure.
    Mirrors measure_fig7.py (lab-frame rise in the closed slip-wall box). Returns (t*, U*) or None."""
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
    dt, mu_b, sigma_T = f("dt"), 1.0 / f("fluid_pp(1)%re(1)"), f("sigma_dtdt")
    Ly = f("y_domain%end") - f("y_domain%beg")
    gradT = abs(f("bc_y%twall_out") - f("bc_y%twall_in")) / Ly if "bc_y%twall_out" in p else 1.0 / Ly
    G = abs(sigma_T * gradT)
    U_r, t_r = G * 0.5 / mu_b, mu_b / G
    cells = nx * ny * nz
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    t_star, U_star = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        vy = fld(3) / (fld(0) + fld(1))
        c = np.clip(fld(nvars - 2), 0.0, 1.0)  # color c; T_s is the last conserved variable
        t_star.append(s * dt / t_r)
        U_star.append(((c * vy).sum() / c.sum()) / U_r)
    return np.array(t_star), np.array(U_star)


def make_figure():
    """Overlay U*(t*) for both grids against the Nas & Tryggvason peak."""
    curves = []
    for v in VARIANTS:
        cur = curve(os.path.join(RUNS, v["name"]))
        if cur is not None:
            curves.append((f"{v['nx']} cells/width ({v['nx'] / 2:.0f}/$D$)", v["color"], cur))
    if not curves:
        print("  (no runs found; nothing to plot)")
        return
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for lab, col, (t, uu) in curves:
        ax.plot(t, uu, "-", color=col, lw=1.7, label="MFC " + lab)
    ax.axhline(NAS_TRYGGVASON_PEAK, ls="--", color="0.2", lw=1.4, label=rf"Nas & Tryggvason peak $\approx$ {NAS_TRYGGVASON_PEAK:.2f}")
    ax.axhline(0.0, ls="-", color="0.6", lw=0.8)
    ax.set_xlabel(r"$t^* = t / t_r$  ($t_r = \mu_b / |\sigma_T \nabla T|$)")
    ax.set_ylabel(r"$U^* = U / U_r$")
    ax.set_title("Fig 7: finite-Ma thermocapillary migration, grid convergence (Re=5, Ma=20, Ca=0.01666)")
    ax.set_xlim(left=0.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, "fig7_migration_2D.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def load_summary():
    path = os.path.join(RESULTS, "fig7_summary.json")
    return json.load(open(path)) if os.path.isfile(path) else {}


def save_summary(summary):
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(summary, open(os.path.join(RESULTS, "fig7_summary.json"), "w"), indent=2)


def print_table(summary):
    print("\n" + "=" * 72)
    print(f"{'variant':>10} {'cells/D':>8} {'peak U*':>9} {'@t*':>6} {'terminal U*':>12} {'N&T peak':>9}")
    for name, r in sorted(summary.items(), key=lambda kv: kv[1]["cells_per_D"]):
        print(f"{name:>10} {r['cells_per_D']:>8.1f} {r['peak']:>9.4f} {r['t_peak_tr']:>6.2f} {r['terminal']:>12.4f} {r['nas_tryggvason_peak']:>9.2f}")
    if len(summary) == 2:
        peaks = [r["peak"] for r in summary.values()]
        print(f"\n  grid-to-grid peak spread: {100 * abs(peaks[0] - peaks[1]) / max(peaks):.1f}%   [Samareh: 1.2% (grids), 1.7% vs N&T]")
    print("=" * 72)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "remeasure":
        summary = {}
        for v in VARIANTS:
            workdir = os.path.join(RUNS, v["name"])
            if os.path.isfile(os.path.join(workdir, "simulation.inp")):
                print(f"\n>>> remeasure {v['name']}", flush=True)
                res = measure(workdir)
                if res is not None:
                    summary[v["name"]] = res
        save_summary(summary)
    else:
        summary = load_summary()
        for v in VARIANTS:
            if run_variant(v):
                res = measure(os.path.join(RUNS, v["name"]))
                if res is not None:
                    summary[v["name"]] = res
                    save_summary(summary)  # checkpoint after every successful variant
                    make_figure()

    make_figure()
    print_table(summary)
    print(f"\nwrote {os.path.join(RESULTS, 'fig7_summary.json')} and figure in {RESULTS}/")
