#!/usr/bin/env python3
"""Resume the paused muscl_thinc TC2 run from its last on-disk checkpoint.

The muscl_thinc variant of the THINC experiment (see sweep_recon.py and
[[recon-tc2-thinc-experiment-state]]) was stopped cleanly partway through; its restart_data
checkpoints (lustre_<step>.dat, parallel-IO, one per t_step_save) are intact on disk. MFC
resumes a fixed-dt run by setting `t_step_start` to a saved step: simulation then reads
restart_data/lustre_<t_step_start>.dat and integrates to t_step_stop (m_start_up.fpp:404-406,
926 -- restart-read branch). This script finds the latest checkpoint, regenerates the
muscl_thinc case into the run dir (idempotent -- same MUSCL/int_comp patch as sweep_recon),
injects t_step_start, and launches ONLY the simulation target (no pre_process -> restart_data
is left untouched; no rmtree -> checkpoints are preserved).

t_step_start is a runtime namelist param, NOT part of the compiled analytic IC, so the build is
unchanged -- ./mfc.sh run detects no source change and skips recompilation, then continues from
the checkpoint. To resume the muscl baseline instead you would not need this (it already finished).

Prereq: source the NVHPC toolchain first (GPU is the only backend here):
    source ../../.nighthawk_gpu_env.sh && python3 resume_thinc.py

After it completes: verify simulation.inp still has recon_type=2 + int_comp=T, then run
measure.py (fig7) and band_recon.py for the comparison.
"""

import glob
import os
import re
import subprocess
import sys

import sweep_recon as S  # reuse make_case (the WENO->MUSCL + int_comp=T patch)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNDIR = os.path.join(HERE, "runs", "recon_tc2", "muscl_thinc")
CASE = "case_Ma_20.py"
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}


def latest_checkpoint(rd):
    """Largest numeric restart step in rd (ignores the lustre_*_cb.dat grid files)."""
    steps = []
    for f in glob.glob(os.path.join(rd, "lustre_*.dat")):
        m = re.fullmatch(r"lustre_(\d+)\.dat", os.path.basename(f))
        if m:
            steps.append(int(m.group(1)))
    return max(steps) if steps else None


def make_resume_case(dst, n_start):
    """Regenerate the muscl_thinc case (same patch as the sweep) and inject t_step_start=n_start."""
    S.make_case(dst, "muscl_thinc")
    text = open(dst).read()
    inject = f'data["t_step_start"] = {n_start}\nprint(json.dumps(data))'
    text, n = re.subn(r"(?m)^print\(json\.dumps\(data\)\)\s*$", inject, text, count=1)
    if n != 1:
        sys.exit(f"  inject FAILED for {dst}: replaced {n} print() lines (expected 1)")
    open(dst, "w").write(text)


def main():
    rd = os.path.join(RUNDIR, "restart_data")
    if not os.path.isdir(rd):
        sys.exit(f"no restart_data in {RUNDIR}; nothing to resume (run sweep_recon.py muscl_thinc instead)")
    n_start = latest_checkpoint(rd)
    if not n_start:
        sys.exit(f"no usable checkpoint (>0) in {rd}; rerun fresh with sweep_recon.py muscl_thinc")

    dst = os.path.join(RUNDIR, CASE)
    make_resume_case(dst, n_start)
    # Confirm t_step_stop from the generated case so we don't "resume" an already-finished run.
    js = subprocess.run([sys.executable, dst], capture_output=True, text=True, check=True).stdout
    import json

    cfg = json.loads(js)
    t_stop = int(cfg["t_step_stop"])
    if n_start >= t_stop:
        sys.exit(f"latest checkpoint {n_start} >= t_step_stop {t_stop}: run is already complete, nothing to resume")
    print(f"resuming muscl_thinc from checkpoint {n_start} -> {t_stop} ({100 * n_start / t_stop:.0f}% done)  recon_type={cfg.get('recon_type')} int_comp={cfg.get('int_comp')}")

    rel = os.path.relpath(dst, REPO)
    # -t simulation only: reads restart_data/lustre_<n_start>.dat, no pre_process, no rmtree.
    cmd = ["./mfc.sh", "run", rel, "-t", "simulation", "--no-debug", "-j", "8", "-n", "1", "--gpu", "acc"]
    print(f">>> {' '.join(cmd)}", flush=True)
    p = subprocess.run(cmd, cwd=REPO, env={**os.environ, **NOBIND}, check=False)
    if p.returncode != 0:
        sys.exit(f"resume run FAILED (exit {p.returncode}); inspect {RUNDIR}/MFC.out")
    # Guard (sweep_recon's gotcha #2): confirm the regenerated inp is still MUSCL+THINC.
    if not S.verify_inp(RUNDIR, "muscl_thinc"):
        sys.exit("resume produced a non-MUSCL/THINC simulation.inp -- do not trust results")
    print(f"  OK -> {RUNDIR} (resumed to t_step_stop)")


if __name__ == "__main__":
    main()
