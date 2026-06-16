#!/usr/bin/env python3
"""Re-run the 5 production variants behind the two curated validation figures, with the balanced-IC
fix now default in both case files, then regenerate the figures via plot_curves.py.

Multi-rank (the tc3 babysitter's `pkill -x prterun` is dormant while tc3 runs steadily). Pinned to
cores 16+ and with MPI binding off, to stay clear of tc3 (taskset on cores 0-11). Runs sequentially so
the two per-case pre_process builds don't collide; variants of the same case reuse the build.

Usage:  python3 run_reproduce.py
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")

# (case file, run dir name, MPI ranks, env) -- matches run_fig7.py / run_validation.py production configs.
VARIANTS = [
    ("case_fig7.py", "fig7_w064", 8, {"FIG7_NX": "64"}),
    ("case_fig7.py", "fig7_w128", 16, {"FIG7_NX": "128"}),
    ("case.py", "fig5_2D_w064", 6, {"SAMAREH_NX": "64", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
    ("case.py", "fig5_2D_w128", 16, {"SAMAREH_NX": "128", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
    ("case.py", "fig5_2D_w256", 64, {"SAMAREH_NX": "256", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
]

PIN = ["taskset", "-c", "16-255"]  # keep off tc3's cores 0-11
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}  # don't let prterun pin to specific cores


def main():
    for case, name, ranks, env in VARIANTS:
        wd = os.path.join(RUNS, name)
        if os.path.isdir(wd):
            shutil.rmtree(wd)
        os.makedirs(wd)
        shutil.copy(os.path.join(HERE, case), os.path.join(wd, case))
        rel = os.path.relpath(os.path.join(wd, case), REPO)
        e = {**os.environ, **NOBIND, **env}
        print(f"\n>>> {name}: {case} ranks={ranks}", flush=True)
        p = subprocess.run(PIN + ["./mfc.sh", "run", rel, "-n", str(ranks)], cwd=REPO, env=e, capture_output=True, text=True, check=False)
        if p.returncode != 0:
            print(f"FAILED {name} (exit {p.returncode}):")
            print("\n".join((p.stdout + p.stderr).splitlines()[-25:]))
            sys.exit(1)
        print(f"  {name} done")

    print("\n>>> regenerating curated figures", flush=True)
    r = subprocess.run([sys.executable, "plot_curves.py"], cwd=HERE, capture_output=True, text=True, check=False)
    print(r.stdout + r.stderr)
    print("REPRODUCE_DONE" if r.returncode == 0 else "PLOT_FAILED")


if __name__ == "__main__":
    main()
