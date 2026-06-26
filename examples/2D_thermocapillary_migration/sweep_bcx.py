#!/usr/bin/env python3
"""Side-wall boundary-condition sweep on the TC1 Ma=0.1 w064 conduction case.

The "side" walls are the x walls (bc_x%beg/end): the drop migrates in +y, and the isothermal
gradient walls (cold floor / hot ceiling) are on y, so the x walls are the lateral boundaries
the box closes (or opens) on. The committed case uses bc_x = -2 (reflective). This sweeps the
side walls over the four wall/open options MFC supports, holding the y walls (isothermal
gradient) and everything else fixed, so the ONLY change is how the lateral boundary treats the
flow:

    -2   reflective              (the committed baseline; mirror-symmetric closed box)
    -3   ghost-cell extrapolation (OPEN side walls -- zero-gradient outflow, lets return flow leave)
    -15  slip wall               (dedicated free-slip wall)
    -16  no-slip wall            (stationary wall, vb/ve = 0 -> lateral viscous drag near the walls)

bc_x is a runtime namelist parameter, NOT part of the compiled analytic IC (rho/pres/cf strings
are bc-independent), so all four variants share ONE build -- the first point compiles, the rest
reuse it. `./mfc.sh run` handles this; do NOT pass --no-build (the analytic IC still needs the
first build). Runs on the single Tesla V100 (--gpu acc --no-debug, -n 1), variants sequential.

The committed case computes round(2.0*t_r/...) (a short conduction-study window); this patches it
to 10*t_r so the curves span the full case1_fig5.png / Samareh range. The side-wall effect is a
SLOW one (closed-box return flow vs. open-box drift grows over many t_r), so the long window is
the diagnostic one here -- at 2*t_r the variants would barely separate.

x is uniform in the imposed field (T = T(y) only), so an open (-3) or adiabatic wall on x is
physically consistent -- the isothermal-wall caveat in the docs is about the y walls, which we do
not touch. -16/-15 need viscous=T and weno_avg=T, both already set in the case.

Prereq: source the NVHPC toolchain first (this box has no modules; GPU is the only backend here):
    source ../../.nighthawk_gpu_env.sh && python3 sweep_bcx.py

Usage:  python3 sweep_bcx.py                 # full sweep (-2 -3 -15 -16) on the GPU, 10 t_r
        python3 sweep_bcx.py -3 -16           # just the open and no-slip variants
        python3 sweep_bcx.py --tr 2           # short 2 t_r window (quick look)
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASE = "case_Ma_0p1.py"
TR_WINDOW = 10.0  # capillary-thermal times (matches case1_fig5.png / the piinf sweep window)

# Side-wall BCs to compare, with the short human label used by plot_bcx.py.
BCS = [-2, -3, -15, -16]
BC_LABEL = {
    -2: "reflective",
    -3: "ghost extrap (open)",
    -15: "slip wall",
    -16: "no-slip wall",
}

# Don't let prterun pin to cores (single V100, single-user 8-core/16-thread Ryzen; run.py's
# "16-255" taskset range is invalid here). NOBIND keeps prterun from pinning all ranks to one core.
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}


def tag_of(bc):
    """bc int -> run-dir leaf tag: -2 -> 'bcm2', -16 -> 'bcm16' ('m' = minus)."""
    return "bcm" + str(abs(int(bc)))


def make_case(dst, bc):
    """Copy the committed case into dst, patch both x side walls to `bc`, and stretch the t_r
    window. Nx stays 64 (w064); the y walls (isothermal gradient) and everything else are
    untouched."""
    text = open(os.path.join(HERE, CASE)).read()
    text, n_win = re.subn(r"round\(2\.0 \* t_r / mydt\)", f"round({TR_WINDOW} * t_r / mydt)", text, count=1)
    text, n_b = re.subn(r'(?m)^(\s*"bc_x%(?:beg|end)":\s*)-?\d+,', rf"\g<1>{bc},", text)
    if not (n_win == 1 and n_b == 2):
        sys.exit(f"  PATCH FAILED for {dst}: window={n_win} (want 1) bc_x={n_b} (want 2)")
    open(dst, "w").write(text)


def run_variant(bc):
    wd = os.path.join(HERE, "runs", "tc1", "ma0p1", "w064", "bcx", tag_of(bc))
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    dst = os.path.join(wd, CASE)
    make_case(dst, bc)

    rel = os.path.relpath(dst, REPO)
    # -j 8: the first point's analytic IC build runs while the GPU is idle; later points reuse it
    # (bc_x is runtime-only), so the build is a near-instant no-op. Explicit --gpu acc so the run
    # does not inherit a stale CPU lock.yaml default.
    cmd = ["./mfc.sh", "run", rel, "--no-debug", "-j", "8", "-n", "1", "--gpu", "acc"]
    print(f"\n>>> bc_x={bc} ({BC_LABEL[bc]}):  {TR_WINDOW:g} t_r  -> {rel}", flush=True)
    p = subprocess.run(
        cmd,
        cwd=REPO,
        env={**os.environ, **NOBIND},
        capture_output=True,
        text=True,
        check=False,
    )
    log = os.path.join(wd, "sweep_bcx.log")
    open(log, "w").write(p.stdout + "\n===STDERR===\n" + p.stderr)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}); tail of output:")
        print("\n".join((p.stdout + p.stderr).splitlines()[-20:]))
        return False
    for line in reversed(p.stdout.splitlines()):
        if "Total-time" in line or "Time/step" in line:
            print(f"  {line.strip()}")
            break
    print(f"  OK -> {wd}")
    return True


def main():
    args = sys.argv[1:]
    if "--tr" in args:
        i = args.index("--tr")
        global TR_WINDOW
        TR_WINDOW = float(args[i + 1])
        del args[i : i + 2]
    values = [int(a) for a in args] or BCS
    bad = [v for v in values if v not in BC_LABEL]
    if bad:
        sys.exit(f"unknown bc_x code(s) {bad}; choose from {list(BC_LABEL)}")
    print(f"side-wall (bc_x) sweep: {values}  (1 rank GPU each, {TR_WINDOW:g} t_r, w064 Ma=0.1)")
    results = {bc: run_variant(bc) for bc in values}
    print("\n=== side-wall sweep run summary ===")
    for bc, ok in results.items():
        print(f"  bc_x={bc:>4} ({BC_LABEL[bc]:<20}): {'OK' if ok else 'FAILED'}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
