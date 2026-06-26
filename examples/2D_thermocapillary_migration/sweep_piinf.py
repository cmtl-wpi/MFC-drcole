#!/usr/bin/env python3
"""pi_inf (stiffened-gas stiffness) sweep on the TC1 Ma=0.1 w064 conduction case.

The committed case1_fig5.png sweeps GRID at fixed EOS. This sweeps the EOS stiffness p_inf
instead, at fixed grid (w064) and fixed smoothing (smooth_coeff default), so the ONLY thing
that changes is the liquid sound speed c = sqrt(gam*(p0+p_inf)/rho) -- i.e. the artificial
compressibility / effective Mach number. The case file is fully parameterised off p_inf:
cv, k_b, dt and the analytic density profile all rederive, and the thermal diffusivity
alpha_b = G*r^2/(mu*Ma) is p_inf-INDEPENDENT, so the Marangoni forcing and conduction physics
are held fixed. The sweep therefore isolates compressibility sensitivity of the rise velocity.

Each p_inf gives a different compiled analytic IC (rho_expr embeds (p0+p_inf)), so every variant
gets its own build -- `./mfc.sh run` handles this; do NOT use --no-build (stale hash for analytic
ICs). Runs on the single Tesla V100 (--gpu acc --no-debug, -n 1), variants sequential.

The committed case computes round(2.0*t_r/...) (a short conduction-study window); this patches it
back to 10*t_r so the curves span the full case1_fig5.png / Samareh range and show the late droop.

Prereq: source the NVHPC toolchain first (this box has no modules) -- it provides nvfortran AND
the hpcx MPI wrapper, which the CPU (--no-gpu) build also needs on PATH:
    source ../../.nighthawk_gpu_env.sh && python3 sweep_piinf.py

GPU is ~7x slower than CPU on this w064 grid (overhead/occupancy-bound -- see memory); use --cpu
for the cheap points and reserve the GPU for whichever points you specifically want there.

Usage:  python3 sweep_piinf.py                    # full ladder on GPU, 10 t_r
        python3 sweep_piinf.py 512 2048 4096      # high-stiffness end on GPU
        python3 sweep_piinf.py --tr 2 24375 219375    # 2 t_r window for the very stiff (c~500/1500) end
        python3 sweep_piinf.py --cpu 0.5 2 8 32 128   # low/mid end on CPU (6 ranks, ~7x faster)
        python3 sweep_piinf.py --cpu --ranks 16 128   # CPU with a custom rank count

`--tr N` sets the simulation window in capillary-thermal times (default 10). The acoustic dt
shrinks as 1/c, so steps ~ c*N; the very stiff points (c~500-1500) are only affordable at the
short 2 t_r window, which still captures the early overshoot/plateau (peaks by ~1-2 t_r) but
truncates the late droop. Such a curve sits on the same figure as the 10 t_r ladder, just ending
at t/t_r=2 -- expected, not a bug.
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASE = "case_Ma_0p1.py"
TR_WINDOW = 10.0  # capillary-thermal times (matches case1_fig5.png / the droop window)
# Don't let prterun pin to cores; CPU runs additionally stay off cores 0-15 (taskset) and
# oversubscribe is harmless. NOBIND is applied to GPU runs too (single V100, no core pinning wanted).
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}
# No taskset pinning: nighthawk is a single-user 8-core/16-thread Ryzen (CPUs 0-15), so run.py's
# "16-255" range is invalid here and there is no neighbour to avoid. Let the OS place the ranks
# (NOBIND keeps prterun from pinning them all to one core).

# Geometric ladder centred on the committed baseline (p_inf=32, c~20). x4 per step spans ~3.6
# decades in p_inf and ~15x in sound speed/Mach. p0+p_inf stays positive everywhere (p0=8), so the
# stiffened-gas EOS is valid throughout; the low end is the most compressible (highest Mach).
LADDER = [0.5, 2.0, 8.0, 32.0, 128.0, 512.0, 2048.0]


def tag_of(p_inf):
    """p_inf -> run-dir leaf tag: 0.5->p0p5, 2.0->p2, 2048.0->p2048."""
    s = str(int(p_inf)) if float(p_inf).is_integer() else str(p_inf).replace(".", "p")
    return "p" + s


def make_case(dst, p_inf):
    """Copy the committed case into dst and patch p_inf and the t_r window. Nx stays 64 (w064)."""
    text = open(os.path.join(HERE, CASE)).read()
    text, n_win = re.subn(r"round\(2\.0 \* t_r / mydt\)", f"round({TR_WINDOW} * t_r / mydt)", text, count=1)
    text, n_pi = re.subn(r"(?m)^p_inf, p0 = [^,]+,", f"p_inf, p0 = {p_inf},", text, count=1)
    if not (n_win == n_pi == 1):
        sys.exit(f"  PATCH FAILED for {dst}: window={n_win} p_inf={n_pi} (expected 1 each)")
    open(dst, "w").write(text)


def run_variant(p_inf, cpu, ranks):
    wd = os.path.join(HERE, "runs", "tc1", "ma0p1", "w064", "piinf", tag_of(p_inf))
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    dst = os.path.join(wd, CASE)
    make_case(dst, p_inf)

    rel = os.path.relpath(dst, REPO)
    # CPU: taskset off cores 0-15 + --no-gpu. GPU: single V100, --gpu acc --no-debug (the memory's
    # required Release-GPU selector). Both pass explicit backend flags so neither inherits the
    # other's build/lock.yaml default -- safe to run a CPU batch concurrently with the GPU sweep.
    backend = "CPU" if cpu else "GPU"
    # -j 8: each p_inf is a distinct compiled IC -> a fresh per-point build. Builds run while the GPU
    # is idle (the sweep is sequential: build N, run N, build N+1), so parallel compile is free here.
    base = ["./mfc.sh", "run", rel, "--no-debug", "-j", "8", "-n", str(ranks)]
    cmd = base + (["--no-gpu"] if cpu else ["--gpu", "acc"])
    print(f"\n>>> p_inf={p_inf}: ranks={ranks} ({backend}, {TR_WINDOW:g} t_r)  -> {rel}", flush=True)
    p = subprocess.run(
        cmd,
        cwd=REPO,
        env={**os.environ, **NOBIND},
        capture_output=True,
        text=True,
        check=False,
    )
    log = os.path.join(wd, "sweep_piinf.log")
    open(log, "w").write(p.stdout + "\n===STDERR===\n" + p.stderr)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}); tail of output:")
        print("\n".join((p.stdout + p.stderr).splitlines()[-20:]))
        return False
    # Echo the final progress line so the step rate / wall time is visible.
    for line in reversed(p.stdout.splitlines()):
        if "Total-time" in line or "Time/step" in line:
            print(f"  {line.strip()}")
            break
    print(f"  OK -> {wd}")
    return True


def main():
    global TR_WINDOW
    args = sys.argv[1:]
    cpu = "--cpu" in args
    args = [a for a in args if a != "--cpu"]
    ranks = 6 if cpu else 1  # CPU: 6 ranks (matches the fig5 CPU runs); GPU: 1 rank on the single V100
    if "--ranks" in args:
        i = args.index("--ranks")
        ranks = int(args[i + 1])
        del args[i : i + 2]
    if "--tr" in args:
        i = args.index("--tr")
        TR_WINDOW = float(args[i + 1])
        del args[i : i + 2]
    values = [float(a) for a in args] or LADDER
    backend = "CPU" if cpu else "GPU"
    print(f"pi_inf sweep: {values} ({ranks} rank {backend} each, {TR_WINDOW:g} t_r, w064 Ma=0.1)")
    results = {p: run_variant(p, cpu, ranks) for p in values}
    print("\n=== pi_inf sweep run summary ===")
    for p, ok in results.items():
        print(f"  p_inf={p:>8}: {'OK' if ok else 'FAILED'}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
