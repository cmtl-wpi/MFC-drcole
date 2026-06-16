#!/usr/bin/env python3
"""Re-run the fig5 runs that plot_curves.fig5 actually reads, now with the balanced IC default, then
regenerate the curated figures. The main curve is the long-window 64-grid (fig5_2D_w064_tr10, tr=10);
w256 is the tr=2 convergence overlay that the babysitter killed mid-run. Multi-rank, pinned off tc3.

Usage:  python3 run_fig5fix.py
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")

VARIANTS = [
    ("case.py", "fig5_2D_w064_tr10", 6, {"SAMAREH_NX": "64", "SAMAREH_WALL": "1", "SAMAREH_TR": "10"}),
    ("case.py", "fig5_2D_w256", 64, {"SAMAREH_NX": "256", "SAMAREH_WALL": "1", "SAMAREH_TR": "2.0"}),
]
PIN = ["taskset", "-c", "16-255"]
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}


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
        p = subprocess.run(PIN + ["./mfc.sh", "run", rel, "-n", str(ranks)],
                           cwd=REPO, env=e, capture_output=True, text=True, check=False)
        if p.returncode != 0:
            print(f"FAILED {name} (exit {p.returncode}):")
            print("\n".join((p.stdout + p.stderr).splitlines()[-25:]))
            continue  # keep going so the main curve still regenerates even if w256 is killed again
        print(f"  {name} done")

    print("\n>>> regenerating curated figures", flush=True)
    r = subprocess.run([sys.executable, "plot_curves.py"], cwd=HERE, capture_output=True, text=True, check=False)
    print(r.stdout + r.stderr)
    print("FIG5FIX_DONE" if r.returncode == 0 else "PLOT_FAILED")


if __name__ == "__main__":
    main()
