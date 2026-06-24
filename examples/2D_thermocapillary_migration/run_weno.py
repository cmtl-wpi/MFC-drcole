#!/usr/bin/env python3
"""WENO-order study on the TC1 Ma=0.1 w064 conduction droop.

Memory records that the TC1 rise-velocity droop is numerical interface diffusion: WENO smears
the passively-advected color band and the diffuse-interface CSF cannot re-sharpen it, so the
Marangoni force weakens over time. Ma, grid, and smooth_coeff were already swept; WENO order is
the one untested numerical-diffusion knob. Hypothesis: higher order -> sharper band -> less droop.

Four w064 variants, 10 t_r, identical otherwise (dt depends only on grid+Ma, so all share the
same timestepping -- a controlled comparison). mp_weno is only valid at order 5, so the o5_nompw
control separates "WENO order" from "the mp_weno monotonicity limiter".

    name      weno_order  mp_weno   role
    o3        3           F         low order (more diffusive)
    o5        5           T         baseline (matches the committed case_Ma_0p1.py)
    o5_nompw  5           F         mp_weno control (isolates the order effect at fixed order)
    o7        7           F         high order (less diffusive)

Runs on the GPU (OpenACC, --gpu acc) at -np 1 -- this box has a single Tesla V100, so the four
variants run sequentially. pre_process/post_process stay on CPU; only simulation is offloaded.
The first run may build pre_process for the analytic IC; the rest reuse the cached build.

Usage:  python3 run_weno.py            # run all four
        python3 run_weno.py o7 o3      # run a subset
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASE = "case_Ma_0p1.py"
RANKS = 1  # single Tesla V100 -> one GPU rank; the four variants run sequentially
TR_WINDOW = 10.0  # capillary-thermal times to run (matches the existing droop diagnostic window)
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}

# (run-dir leaf name, weno_order, mp_weno)
VARIANTS = [
    ("o5", 5, "T"),
    ("o3", 3, "F"),
    ("o5_nompw", 5, "F"),
    ("o7", 7, "F"),
]


def make_case(dst, weno_order, mp_weno):
    """Copy the committed case into dst and patch weno_order, mp_weno, and the t_r window."""
    text = open(os.path.join(HERE, CASE)).read()
    text, n_win = re.subn(r"round\(2\.0 \* t_r / mydt\)", f"round({TR_WINDOW} * t_r / mydt)", text, count=1)
    text, n_ord = re.subn(r'("weno_order":\s*)\d+', rf"\g<1>{weno_order}", text, count=1)
    text, n_mpw = re.subn(r'("mp_weno":\s*)"[TF]"', rf'\g<1>"{mp_weno}"', text, count=1)
    if not (n_win == n_ord == n_mpw == 1):
        sys.exit(f"  PATCH FAILED for {dst}: window={n_win} order={n_ord} mp_weno={n_mpw} (expected 1 each)")
    open(dst, "w").write(text)


def run_variant(name, weno_order, mp_weno):
    wd = os.path.join(HERE, "runs", "tc1", "ma0p1", "w064", "weno", name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    dst = os.path.join(wd, CASE)
    make_case(dst, weno_order, mp_weno)

    rel = os.path.relpath(dst, REPO)
    print(f"\n>>> {name}: weno_order={weno_order} mp_weno={mp_weno} ranks={RANKS} (GPU)  -> {rel}", flush=True)
    p = subprocess.run(
        ["./mfc.sh", "run", rel, "--gpu", "acc", "-n", str(RANKS)],
        cwd=REPO,
        env={**os.environ, **NOBIND},
        capture_output=True,
        text=True,
        check=False,
    )
    log = os.path.join(wd, "run_weno.log")
    open(log, "w").write(p.stdout + "\n===STDERR===\n" + p.stderr)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}); tail of output:")
        print("\n".join((p.stdout + p.stderr).splitlines()[-20:]))
        return False
    print(f"  OK -> {wd}")
    return True


def main():
    want = set(sys.argv[1:])
    variants = [v for v in VARIANTS if not want or v[0] in want]
    print(f"WENO study: running {[v[0] for v in variants]} ({RANKS} ranks each, {TR_WINDOW:g} t_r)")
    results = {name: run_variant(name, wo, mpw) for name, wo, mpw in variants}
    print("\n=== WENO study run summary ===")
    for name, ok in results.items():
        print(f"  {name:>10}: {'OK' if ok else 'FAILED'}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
