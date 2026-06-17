#!/usr/bin/env python3
"""Run the thermocapillary validation cases and aggregate the results -- one driver for all three.

Each variant runs a COPY of its case in its own runs/<name>/ directory (MFC writes output next to the
case file) via `./mfc.sh run`, is measured with measure.py (auto-detected mode passed explicitly), and
its RESULT_JSON is aggregated into results/<target>_summary.json. The curated overlay figures are
regenerated at the end via plot_curves.py and plot_samareh_style.py.

  fig5  case.py       Sec 4.1.1 / Fig 5   -- 2D zero-Ma rise, grid convergence (v/v_YGB -> ~0.80)
  fig7  case_fig7.py  Sec 4.1.2 / Fig 7   -- 2D finite-Ma migration (U*/U_r, peak ~0.13)
  tc3   case_tc3.py   Sec 4.2 / Figs 8,13 -- 3D large-Ma + mu(T) (rise mm/s vs distance from wall)
  all   fig5 + fig7

Runs launch pinned to cores 16-255 with MPI binding off (the safe pattern on a shared box; harmless
on an idle one), sequentially, so the per-case pre_process builds don't collide and same-case variants
reuse the build. Rank counts respect MFC's decomposition rule (>= 25 cells per split dim).

Usage (invokes mpirun, so run from a normal shell):
    python3 run.py <fig5|fig7|tc3|all> [run|remeasure]   (default: fig5 run)

`remeasure` rebuilds the summary from existing runs/ WITHOUT running any simulation.
"""

import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")
RESULTS = os.path.join(HERE, "results")
PIN = ["taskset", "-c", "16-255"]  # keep off cores 0-15 (a neighbour's job may live there)
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}  # don't let prterun pin to specific cores

# Per-target config. variants are (run-dir name, MPI ranks, env overrides for the case). The summary
# filenames preserve the historical names so plot scripts / notes that read them still resolve.
TARGETS = {
    "fig5": dict(case="case.py", mode="fig5", summary="summary.json", variants=[
        ("fig5_2D_w064", 6, {"SAMAREH_NX": "64", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
        ("fig5_2D_w128", 16, {"SAMAREH_NX": "128", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
        ("fig5_2D_w256", 64, {"SAMAREH_NX": "256", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
    ]),
    "fig7": dict(case="case_fig7.py", mode="fig7", summary="fig7_summary.json", variants=[
        ("fig7_w064", 8, {"FIG7_NX": "64"}),
        ("fig7_w128", 16, {"FIG7_NX": "128"}),
    ]),
    "tc3": dict(case="case_tc3.py", mode="tc3", summary="tc3_summary.json", variants=[
        ("tc3_run", 16, {}),  # 3D large-Ma + mu(T); set ranks/env to match your grid before running
    ]),
}

# Headline-number keys to show per mode in the end-of-run table.
TABLE_KEYS = {
    "fig5": ("cells_per_D", "ratio_plateau", "overshoot", "ratio_final"),
    "fig7": ("cells_per_D", "peak", "t_peak_tr", "terminal"),
    "tc3": ("cells", "dist_end_mm", "peak_rise_velocity_mms"),
}


def run_variant(case_file, name, ranks, env):
    """Run one case variant in its own runs/<name>/ directory. Returns True on success."""
    wd = os.path.join(RUNS, name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    shutil.copy(os.path.join(HERE, case_file), os.path.join(wd, case_file))
    rel = os.path.relpath(os.path.join(wd, case_file), REPO)
    e = {**os.environ, **NOBIND, **env}
    print(f"\n>>> {name}: {case_file} ranks={ranks} {env}", flush=True)
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
        for name, ranks, env in cfg["variants"]:
            wd = os.path.join(RUNS, name)
            if run_mode == "run":
                if not run_variant(cfg["case"], name, ranks, env):
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
