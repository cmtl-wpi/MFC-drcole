#!/usr/bin/env python3
"""Run the thermocapillary validation cases, measure them, and refresh the figures.

Each grid variant runs in its own runs/<name>/ dir (MFC writes output next to the case
file), is measured by measure.py, and its result is collected into results/<summary>.json.
The overlay figures are regenerated at the end.

Targets:
  fig5  case_Ma_0.py     -- 2D zero-Ma rise, grid convergence (v/v_YGB -> ~0.80)
  fig7  case_Ma_20.py    -- 2D low-Ma migration (U*/U_r, peak ~0.13)
  tc3   3D case_Ma_1723  -- 3D large-Ma + mu(T)
  all   fig5 + fig7

Usage (invokes mpirun -- run from a normal shell):
    python3 run.py <fig5|fig7|tc3|all> [run|remeasure]   (default: fig5 run)

  run        run each variant, then measure and plot
  remeasure  re-measure existing runs/ and replot, WITHOUT running any simulation

The canonical case files have no grid knob, so a variant's grid is set by copying the case
into its run dir and rewriting its `Nx = <n>` line; the committed case is never touched.
Runs are pinned off cores 0-15 with MPI binding disabled (safe on the shared server).
"""

import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")
RESULTS = os.path.join(HERE, "results")

PIN = ["taskset", "-c", "16-255"]  # keep off cores 0-15 (a neighbour's job may live there)
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}  # don't let prterun pin to cores

# Per target: the case file, the measure.py mode, the summary file, and the grid variants to run.
# Each variant is (run-dir name, MPI ranks, Nx override or None for the case default). Summary
# filenames keep their historical names so the plot scripts and notes that read them resolve.
TARGETS = {
    "fig5": dict(case="case_Ma_0.py", mode="fig5", summary="summary.json", variants=[
        ("fig5_2D_w064", 6, 64),
        ("fig5_2D_w128", 16, 128),
        ("fig5_2D_w256", 64, 256),
    ]),
    "fig7": dict(case="case_Ma_20.py", mode="fig7", summary="fig7_summary.json", variants=[
        ("fig7_w064", 8, 64),
        ("fig7_w128", 16, 128),
    ]),
    "tc3": dict(case="../3D_thermocapillary_migration/case_Ma_1723.py", mode="tc3", summary="tc3_summary.json", variants=[
        ("tc3_run", 16, None),  # set ranks/Nx to match your grid before running
    ]),
}

# Headline numbers to show per mode in the end-of-run table.
TABLE_KEYS = {
    "fig5": ("cells_per_D", "ratio_plateau", "overshoot", "ratio_final"),
    "fig7": ("cells_per_D", "peak", "t_peak_tr", "terminal"),
    "tc3": ("cells", "dist_end_mm", "peak_rise_velocity_mms"),
}


def run_variant(case_file, name, ranks, nx):
    """Run one variant in a fresh runs/<name>/ dir. Returns True on success.

    The case is copied in (basename only, so a ../ source lands flat); if nx is given, the
    copy's `Nx = <n>` line is rewritten. The committed case file is never touched.
    """
    wd = os.path.join(RUNS, name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)

    dst = os.path.join(wd, os.path.basename(case_file))
    shutil.copy(os.path.join(HERE, case_file), dst)
    if nx is not None:
        text, hits = re.subn(r"(?m)^Nx = \d+", f"Nx = {nx}", open(dst).read(), count=1)
        if hits != 1:
            print(f"  WARNING: no `Nx = <int>` line in {name} -- using case default")
        open(dst, "w").write(text)

    print(f"\n>>> {name}: {os.path.basename(case_file)} ranks={ranks} Nx={nx or 'default'}", flush=True)
    rel = os.path.relpath(dst, REPO)
    p = subprocess.run(
        PIN + ["./mfc.sh", "run", rel, "-n", str(ranks)],
        cwd=REPO, env={**os.environ, **NOBIND}, capture_output=True, text=True, check=False,
    )
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}):")
        print("\n".join((p.stdout + p.stderr).splitlines()[-25:]))
        return False
    return True


def measure(wd, mode):
    """Run measure.py <wd> <mode> and return its RESULT_JSON dict (or None on failure)."""
    m = subprocess.run([sys.executable, os.path.join(HERE, "measure.py"), wd, mode], capture_output=True, text=True, check=False)
    print("\n".join(m.stdout.splitlines()[-4:]))
    if m.returncode != 0:
        print(f"  MEASURE FAILED: {m.stderr[-500:]}")
        return None
    for line in m.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


def regenerate_figures():
    """Rebuild the curated overlay figures from whatever runs are on disk."""
    for script in ("plot_samareh_style.py",):
        path = os.path.join(HERE, script)
        if not os.path.isfile(path):
            continue
        r = subprocess.run([sys.executable, path], cwd=HERE, capture_output=True, text=True, check=False)
        print(f"  {script}: {'ok' if r.returncode == 0 else 'FAILED -- ' + r.stderr[-300:]}")


def process_target(target, run_mode):
    """Run (or just re-measure) every variant of one target and update its summary file."""
    cfg = TARGETS[target]
    spath = os.path.join(RESULTS, cfg["summary"])
    summary = json.load(open(spath)) if os.path.isfile(spath) else {}

    for name, ranks, nx in cfg["variants"]:
        wd = os.path.join(RUNS, name)
        if run_mode == "run":
            if not run_variant(cfg["case"], name, ranks, nx):
                continue
        elif not os.path.isfile(os.path.join(wd, "simulation.inp")):
            continue  # remeasure: skip variants that were never run
        res = measure(wd, cfg["mode"])
        if res is not None:
            summary[name] = res
            json.dump(summary, open(spath, "w"), indent=2)  # checkpoint after each variant

    keys = TABLE_KEYS[target]
    print(f"\n=== {target} summary ({len(summary)} runs) -> {spath} ===")
    for name, res in sorted(summary.items()):
        print(f"  {name:>16}: " + "  ".join(f"{k}={res[k]}" for k in keys if k in res))


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "fig5"
    run_mode = sys.argv[2] if len(sys.argv) > 2 else "run"
    targets = ["fig5", "fig7"] if target == "all" else [target]
    if any(t not in TARGETS for t in targets):
        sys.exit(f"unknown target {target!r}; choose from {', '.join(TARGETS)}, all")

    os.makedirs(RESULTS, exist_ok=True)
    for t in targets:
        process_target(t, run_mode)
    regenerate_figures()


if __name__ == "__main__":
    main()
