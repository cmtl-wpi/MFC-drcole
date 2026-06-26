#!/usr/bin/env python3
"""Grid + Marangoni-number sweep on the restyled 3D Samareh Fig-6 case (case.py).

Diagnoses the 0.82 (Nx=64, Ma=1.0) vs Samareh-converged 0.95 gap by separating the two candidate
causes:
  - GRID sweep at fixed Ma=0.1:  Nx in {32, 48, 64}  -> does v/v_YGB climb toward 0.95 as dx->0?
  - Ma  sweep at fixed Nx=32:    Ma in {1.0, 0.1, 0.001} -> is the result Ma-independent (i.e. is T
    already invariant, so Marangoni number is NOT the lever)?

WHY Ma=0.001 is pinned to the coarsest grid: conduction is explicit, so dt <= 0.35*dx^2/(2*d*alpha)
with alpha ~ 1/Ma. At Ma=0.001 the diffusion-limited dt is ~100x below the acoustic limit; that is
only affordable at Nx=32 (~21 h). At Nx=64 the same point is ~28 days. So Ma=0.001 runs at Nx=32.

Nx and Ma are BOTH runtime-only here (the analytic IC rho=rho_coeff/T(y) is Nx/Ma-independent), so
ALL variants share ONE build -- the first point compiles, the rest reuse it. Runs on the single
Tesla V100 (--gpu acc --no-debug, -n 1), sequential, ordered cheap/diagnostic first so the 21 h
Ma=0.001 tail can be cancelled once the earlier points settle the Ma question.

Prereq: source the NVHPC toolchain first (GPU is the only backend here):
    source ../../.nighthawk_gpu_env.sh && python3 sweep_grid_ma.py

The existing Nx=64/Ma=1.0 run (runs/fig6_anchor/nx064, v/v_YGB=0.82) is the (64, 1.0) reference --
not re-run here.

Usage:  python3 sweep_grid_ma.py            # the information-first campaign below
        python3 sweep_grid_ma.py 32:0.1 48:0.1   # explicit Nx:Ma points only
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASE = "case.py"

# (Nx, Ma), ordered cheap/diagnostic first. grid series @ Ma=0.1: 32,48,64.  Ma series @ Nx=32: 1.0,0.1,0.001.
CAMPAIGN = [
    (32, 1.0),  # ~0.3 h  Ma-independence baseline at the coarse grid
    (32, 0.1),  # ~0.3 h  Ma-independence check vs (32,1.0)
    (48, 0.1),  # ~1.7 h  grid point
    (64, 0.1),  # ~6.6 h  grid point (+ Ma check vs the done (64,1.0)=0.82)
    (32, 0.001),  # ~21 h   deep Ma confirmation (coarse grid is the only feasible one)
]
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}


def fmt(v):
    """Number -> dir token: 1.0->1p0, 0.1->0p1, 0.001->0p001."""
    return ("%g" % v).replace(".", "p")


def make_case(dst, nx, ma):
    """Copy case.py into dst and patch the hardcoded Nx and Ma literals (the only two knobs)."""
    text = open(os.path.join(HERE, CASE)).read()
    text, n_nx = re.subn(r"(?m)^Nx = \d+", f"Nx = {nx}", text, count=1)
    text, n_ma = re.subn(r"(?m)^Ma = [\d.]+", f"Ma = {ma}", text, count=1)
    if not (n_nx == 1 and n_ma == 1):
        sys.exit(f"  PATCH FAILED for {dst}: Nx={n_nx} Ma={n_ma} (expected 1 each)")
    open(dst, "w").write(text)


def run_point(nx, ma):
    wd = os.path.join(HERE, "runs", "grid_ma", f"nx{nx:03d}", f"ma{fmt(ma)}")
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    make_case(os.path.join(wd, CASE), nx, ma)
    rel = os.path.relpath(os.path.join(wd, CASE), REPO)
    cmd = ["./mfc.sh", "run", rel, "--no-debug", "-j", "8", "-n", "1", "--gpu", "acc"]
    print(f"\n>>> Nx={nx} Ma={ma}  -> {rel}", flush=True)
    p = subprocess.run(cmd, cwd=REPO, env={**os.environ, **NOBIND}, capture_output=True, text=True, check=False)
    open(os.path.join(wd, "sweep.log"), "w").write(p.stdout + "\n===STDERR===\n" + p.stderr)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}); tail:")
        print("\n".join((p.stdout + p.stderr).splitlines()[-20:]))
        return False
    for line in reversed(p.stdout.splitlines()):
        if "Total-time" in line:
            print(f"  {line.strip()}")
            break
    print(f"  OK -> {wd}")
    return True


def main():
    args = sys.argv[1:]
    if args:
        pts = [(int(a.split(":")[0]), float(a.split(":")[1])) for a in args]
    else:
        pts = CAMPAIGN
    print(f"3D grid+Ma sweep: {pts}  (1 rank GPU each, Samareh 5D^2x7.5D box)")
    results = {(nx, ma): run_point(nx, ma) for nx, ma in pts}
    print("\n=== sweep run summary ===")
    for (nx, ma), ok in results.items():
        print(f"  Nx={nx:>3} Ma={ma:<6}: {'OK' if ok else 'FAILED'}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
