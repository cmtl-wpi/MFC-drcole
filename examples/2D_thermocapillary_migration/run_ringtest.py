#!/usr/bin/env python3
"""Short before/after diagnostic for the acoustic-ring mitigation (FIG7_RINGFIX in case_fig7.py).

Runs case_fig7.py at a coarse grid (n_x=64) over a short window (3 t_r) for three variants:
  v0  before                 : Samareh's uniform-pressure IC (rings)
  v1  balanced IC            : p_in = p_out + sigma/r at t=0
  v2  balanced IC + bulk visc: also fluid_pp%Re(2) (mu_bulk small enough to keep dt acoustic-limited)
then overlays U*(t*) and reports the ring amplitude (std of the detrended U* residual) per variant.

This machine runs other MFC jobs whose babysitter does `pkill -9 -x prterun`, which would also kill
our MPI launcher. So we run each binary as a DETACHED SINGLETON (no prterun): inputs are generated
with `./mfc.sh run --dry-run` (pure Python, no launcher), then the case-optimized pre_process and
simulation binaries are exec'd directly, single-rank, in their own session.

Usage:  python3 run_ringtest.py
"""

import os
import re
import shutil
import subprocess
import time

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_curves import color_weighted_vy  # reuse the exact reader the figures use

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "out", "ringtest")
NX, TR, BV = 64, 1.5, 0.05  # 1.5 t_r covers the ring + overshoot onset; keeps each singleton ~8 min
MAX_WALL = 1200  # s
STAGGER = 22  # s between singleton launches: OMPI falls back to prterun if several MPI_Init race

VARIANTS = [
    dict(k=0, ringfix=0, label="before (uniform p, rings)", color="C3"),
    dict(k=1, ringfix=1, label="balanced IC", color="C0"),
    dict(k=2, ringfix=2, label="balanced IC + bulk visc", color="C2"),
]


def binaries():
    """pre_process: case-built today from case_fig7.py (handles the analytic-T IC, writes the IC).
    simulation: a GENERIC (non-case-optimized) build -- it only reads restart_data, so it needs no
    case-specific compile, and it accepts the FULL namelist that --dry-run emits (the case-optimized
    binary bakes weno_order etc. out of the namelist and would reject the full .inp). Verified both."""
    pp = os.path.join(REPO, "build/install/18cc893151/bin/pre_process")
    sm = os.path.join(REPO, "build/install/655c1b4249/bin/simulation")
    for b in (pp, sm):
        assert os.path.isfile(b), f"missing binary {b}"
    return pp, sm


def gen_inputs(wd, ringfix):
    env = {**os.environ, "FIG7_NX": str(NX), "FIG7_TR": str(TR), "FIG7_RINGFIX": str(ringfix), "FIG7_BULKVISC": str(BV)}
    rel = os.path.relpath(os.path.join(wd, "case_fig7.py"), REPO)
    p = subprocess.run(["./mfc.sh", "run", rel, "--dry-run", "--no-build"], cwd=REPO, env=env, capture_output=True, text=True, check=False)
    if p.returncode != 0:
        print("  dry-run FAILED:\n" + "\n".join((p.stdout + p.stderr).splitlines()[-15:]))
        return False
    return True


def main():
    pp_bin, sm_bin = binaries()
    print(f"pre_process: {pp_bin}\nsimulation : {sm_bin}\n")

    procs = []
    for v in VARIANTS:
        wd = os.path.join(OUT, f"v{v['k']}")
        if os.path.isdir(wd):
            shutil.rmtree(wd)
        os.makedirs(wd)
        shutil.copy(os.path.join(HERE, "case_fig7.py"), os.path.join(wd, "case_fig7.py"))
        if not gen_inputs(wd, v["ringfix"]):
            continue
        # singleton pre_process (blocking; fast) -> grid + step-0 IC
        with open(os.path.join(wd, "pre.log"), "w") as f:
            subprocess.run([pp_bin], cwd=wd, stdout=f, stderr=subprocess.STDOUT, check=False)
        # singleton simulation (detached; no prterun for the babysitter to kill)
        log = open(os.path.join(wd, "sim.log"), "w")
        subprocess.Popen([sm_bin], cwd=wd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        procs.append((v, wd))
        print(f">>> v{v['k']} ({v['label']}) launched", flush=True)
        if v is not VARIANTS[-1]:
            time.sleep(STAGGER)  # let this singleton finish MPI_Init before the next launches

    # poll the sim logs until each reaches its final step
    t_stop = {}
    for v, wd in procs:
        t_stop[wd] = int(re.search(r"t_step_stop\s*=\s*(\d+)", open(os.path.join(wd, "simulation.inp")).read()).group(1))
    t0 = time.time()
    done = set()
    while len(done) < len(procs) and time.time() - t0 < MAX_WALL:
        time.sleep(20)
        msg = []
        for v, wd in procs:
            steps = re.findall(r"Time step\s+(\d+)", open(os.path.join(wd, "sim.log")).read() or "")
            last = int(steps[-1]) if steps else 0
            if last >= t_stop[wd] - 80:
                done.add(wd)
            msg.append(f"v{v['k']}:{last}/{t_stop[wd]}")
        print(f"  [{int(time.time() - t0)}s] " + "  ".join(msg), flush=True)

    # analyze
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    summary = []
    for v, wd in procs:
        out = color_weighted_vy(wd)
        if out is None:
            print(f"  v{v['k']}: no snapshots")
            continue
        t, u_lab, P = out
        mu_b = 1.0 / float(P["fluid_pp(1)%re(1)"])
        G = abs(float(P["sigma_dtdt"]) * (1.0 / (float(P["y_domain%end"]) - float(P["y_domain%beg"]))))
        U_r, t_r = G * 0.5 / mu_b, mu_b / G
        ts, us = t / t_r, u_lab / U_r
        ring = float((us - np.polyval(np.polyfit(ts, us, 5), ts)).std())
        ax.plot(ts, us, ".-", ms=4, lw=1.0, color=v["color"], label=f"{v['label']}  (ring std={ring:.4f})")
        summary.append((v["label"], ring, len(ts)))
        print(f"  v{v['k']}: {len(ts)} snaps, ring std = {ring:.4f}")

    ax.axhline(0.13, ls="--", color="0.4", lw=1.2, label="Nas & Tryggvason peak ~0.13")
    ax.set(xlabel=r"$t^* = t/t_r$", ylabel=r"$U^* = U/U_r$", title=f"Acoustic-ring mitigation, n_x={NX}, {TR} $t_r$ (before vs balanced IC vs +bulk visc)")
    ax.legend(loc="lower right", fontsize=9)
    png = os.path.join(HERE, "results", "ringtest.png")
    fig.savefig(png, dpi=150)
    print(f"\nwrote {png}")
    if summary:
        base = summary[0][1]
        print("\n=== ring amplitude (std of U* residual) ===")
        for label, ring, n in summary:
            print(f"  {label:30s} {ring:.4f}  ({ring / base * 100:5.1f}% of before)")


if __name__ == "__main__":
    main()
