#!/usr/bin/env python3
"""Run the thermocapillary validation cases and aggregate the results -- one driver for all three.

Each variant runs a COPY of its case in its own runs/<name>/ directory (MFC writes output next to the
case file) via `./mfc.sh run`, is measured with measure.py (mode passed explicitly), and its
RESULT_JSON is aggregated into results/<target>_summary.json. The curated overlay figures are
regenerated at the end via plot_curves.py and plot_samareh_style.py.

  fig5  case_Ma_0.py                          Sec 4.1.1 / Fig 5   -- 2D zero-Ma rise, grid convergence (v/v_YGB -> ~0.80)
  fig7  case_Ma_20.py                         Sec 4.1.2 / Fig 7   -- 2D low-Ma migration (U*/U_r, peak ~0.13)
  tc3   ../3D_thermocapillary_migration/case_Ma_1723.py   Sec 4.2 / Figs 8,13 -- 3D large-Ma + mu(T)
  all   fig5 + fig7

The case files are hardcoded canonical examples (no env knobs), so a grid variant is produced by
rewriting the `Nx = <n>` line in that run's COPY of the case before launching; the committed case
files are never touched. The target/mode/summary names stay fig5/fig7/tc3 to match results/*.json and
the plot scripts. Runs launch pinned to cores 16-255 with MPI binding off (the safe pattern on a
shared box; harmless on an idle one), sequentially, so the per-case pre_process builds don't collide
and same-case variants reuse the build. Rank counts respect MFC's decomposition rule (>= 25 cells per
split dim).

Usage (invokes mpirun, so run from a normal shell):
    python3 run.py <fig5|fig7|tc3|all> [run|remeasure]   (default: fig5 run)

`remeasure` rebuilds the summary from existing runs/ WITHOUT running any simulation.
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
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}  # don't let prterun pin to specific cores

# Per-target config. variants are (run-dir name, MPI ranks, Nx override or None for the case default).
# The summary filenames preserve the historical names so plot scripts / notes that read them resolve.
TARGETS = {
    "fig5": dict(
        case="case_Ma_0.py",
        mode="fig5",
        summary="summary.json",
        variants=[
            ("fig5_2D_w064", 6, 64),
            ("fig5_2D_w128", 16, 128),
            ("fig5_2D_w256", 64, 256),
        ],
    ),
    "fig7": dict(
        case="case_Ma_20.py",
        mode="fig7",
        summary="fig7_summary.json",
        variants=[
            ("fig7_w064", 8, 64),
            ("fig7_w128", 16, 128),
        ],
    ),
    "tc3": dict(
        case="../3D_thermocapillary_migration/case_Ma_1723.py",
        mode="tc3",
        summary="tc3_summary.json",
        variants=[
            ("tc3_run", 16, None),  # 3D large-Ma + mu(T); set ranks/Nx to match your grid before running
        ],
    ),
}

# Headline-number keys to show per mode in the end-of-run table.
TABLE_KEYS = {
    "fig5": ("cells_per_D", "ratio_plateau", "overshoot", "ratio_final"),
    "fig7": ("cells_per_D", "peak", "t_peak_tr", "terminal"),
    "tc3": ("cells", "dist_end_mm", "peak_rise_velocity_mms"),
}


def run_variant(case_file, name, ranks, nx):
    """Run one case variant in its own runs/<name>/ directory. Returns True on success.

    The case is copied in (basename only, so a ../ source path lands flat in runs/<name>/); if nx is
    given, the copy's `Nx = <n>` line is rewritten to set the grid (the committed case is untouched).
    """
    wd = os.path.join(RUNS, name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    dst = os.path.join(wd, os.path.basename(case_file))
    shutil.copy(os.path.join(HERE, case_file), dst)
    if nx is not None:
        text = open(dst).read()
        text, n = re.subn(r"(?m)^Nx = \d+", f"Nx = {nx}", text, count=1)
        if n != 1:
            print(f"  WARNING: could not set Nx={nx} in {name} (no `Nx = <int>` line) -- using case default")
        open(dst, "w").write(text)
    rel = os.path.relpath(dst, REPO)
    e = {**os.environ, **NOBIND}
    print(f"\n>>> {name}: {os.path.basename(case_file)} ranks={ranks} Nx={nx or 'default'}", flush=True)
    p = subprocess.run(PIN + ["./mfc.sh", "run", rel, "-n", str(ranks)], cwd=REPO, env=e, capture_output=True, text=True, check=False)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}):")
        print("\n".join((p.stdout + p.stderr).splitlines()[-25:]))
        return False
    return True


def measure(wd, mode):
    """Run measure.py <wd> <mode> and return its RESULT_JSON dict (or None)."""
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
    for script in ("plot_curves.py", "plot_samareh_style.py"):
        path = os.path.join(HERE, script)
        if not os.path.isfile(path):
            continue
        r = subprocess.run([sys.executable, path], cwd=HERE, capture_output=True, text=True, check=False)
        print(f"  {script}: {'ok' if r.returncode == 0 else 'FAILED -- ' + r.stderr[-300:]}")


def main():
    target = sys.argv[1] if len(sys.argv) > 1 else "fig5"
    run_mode = sys.argv[2] if len(sys.argv) > 2 else "run"
    targets = ["fig5", "fig7"] if target == "all" else [target]
    if any(t not in TARGETS for t in targets):
        sys.exit(f"unknown target {target!r}; choose from {', '.join(TARGETS)}, all")
    os.makedirs(RESULTS, exist_ok=True)

    for t in targets:
        cfg = TARGETS[t]
        spath = os.path.join(RESULTS, cfg["summary"])
        summary = json.load(open(spath)) if os.path.isfile(spath) else {}
        for name, ranks, nx in cfg["variants"]:
            wd = os.path.join(RUNS, name)
            if run_mode == "run":
                if not run_variant(cfg["case"], name, ranks, nx):
                    continue
            elif not os.path.isfile(os.path.join(wd, "simulation.inp")):
                continue  # remeasure: skip variants that haven't been run
            res = measure(wd, cfg["mode"])
            if res is not None:
                summary[name] = res
                json.dump(summary, open(spath, "w"), indent=2)  # checkpoint after each variant
        keys = TABLE_KEYS[t]
        print(f"\n=== {t} summary ({len(summary)} runs) -> {spath} ===")
        for nm, r in sorted(summary.items()):
            print(f"  {nm:>16}: " + "  ".join(f"{k}={r[k]}" for k in keys if k in r))

    regenerate_figures()


if __name__ == "__main__":
    main()
