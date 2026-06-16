#!/usr/bin/env python3
"""Turns the saved heat-conduction runs into MP4 animations of the temperature field.

The validation harness (validate.py) keeps each run's checkpoints under runs/ but
never post-processes them, so there is no field data for the viz tool to read.
This script fills that gap: for each run it restores the checkpoints, runs MFC's
post_process to write a Silo database, and calls `./mfc.sh viz --mp4` on the
temperature scalar. Output lands in animations/.

Every case is energy-coupled: temperature is encoded in the density at uniform
pressure and recovered from the stiffened-gas EOS. The case files set `T_wrt`, so
post_process writes that recovered field as `temperature` -- the field viz plots here.

Run from the repo root (post_process needs mpirun, so disable the command sandbox):
  python3 examples/Thermal_Conduction_Validation/animate.py

Needs the runs/ directory populated first (run validate.py for whichever benchmarks
you want). Missing runs are skipped with a note.
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASES_REL = os.path.join("examples", "Thermal_Conduction_Validation", "cases")
CASES = os.path.join(HERE, "cases")
RUNS = os.path.join(HERE, "runs")
ANIM = os.path.join(HERE, "animations")

# label (run dir) , case file , output name , viz variable , viz styling.
# vmin/vmax are fixed so the animation shows the real decay instead of rescaling
# every frame. The 3D mode is sin(kx)sin(ky)sin(kz): its z=L/2 midplane is identically
# flat, so it is sliced at z=L/4 where the structure is strongest. Every run plots the
# EOS-recovered `temperature` field.
ANIMS = [
    ("1d", "case_1d.py", "heat_1d_sine_decay", "temperature", ["--vmin", "9.5", "--vmax", "13.2"]),
    ("2d_mode", "case_2d_mode.py", "heat_2d_mode", "temperature", ["--cmap", "inferno", "--vmin", "7", "--vmax", "13"]),
    ("3d_hotspot", "case_3d_hotspot.py", "heat_3d_hotspot_midplane", "temperature", ["--cmap", "inferno", "--vmin", "10", "--vmax", "15", "--slice-axis", "z"]),
    ("3d_mode", "case_3d_mode.py", "heat_3d_mode_slice", "temperature", ["--cmap", "inferno", "--vmin", "7", "--vmax", "13", "--slice-axis", "z", "--slice-value", "0.25"]),
]


def mfc(*args):
    res = subprocess.run(["./mfc.sh", *args], cwd=REPO, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        print(res.stdout[-2500:], res.stderr[-1500:])
        raise RuntimeError("mfc.sh " + " ".join(args) + " failed")
    return res.stdout


def animate(label, case, name, var, style):
    rundir = os.path.join(RUNS, label)
    if not os.path.isdir(os.path.join(rundir, "restart_data")):
        print(f"  skip {label}: no runs/{label}/restart_data (run validate.py first)")
        return
    # Stage this run's checkpoints where post_process looks for them.
    shutil.rmtree(os.path.join(CASES, "restart_data"), ignore_errors=True)
    shutil.rmtree(os.path.join(CASES, "silo_hdf5"), ignore_errors=True)
    shutil.copytree(os.path.join(rundir, "restart_data"), os.path.join(CASES, "restart_data"))
    shutil.copy(os.path.join(rundir, "simulation.inp"), CASES)

    print(f"  {label}: post_process -> silo")
    mfc("run", os.path.join(CASES_REL, case), "-t", "post_process", "--no-build", "-n", "1")
    print(f"  {label}: viz -> animations/{name}.mp4")
    mfc("viz", CASES, "--var", var, "--step", "all", "--mp4", "--fps", "8", *style, "-o", ANIM)
    shutil.move(os.path.join(ANIM, f"{var}.mp4"), os.path.join(ANIM, f"{name}.mp4"))


if __name__ == "__main__":
    os.makedirs(ANIM, exist_ok=True)
    wanted = [a for a in sys.argv[1:] if not a.startswith("--")]
    for label, case, name, var, style in ANIMS:
        if wanted and label not in wanted:
            continue
        animate(label, case, name, var, style)
    print(f"done -> {ANIM}")
