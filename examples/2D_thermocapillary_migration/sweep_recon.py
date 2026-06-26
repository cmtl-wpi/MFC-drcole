#!/usr/bin/env python3
"""Reconstruction-scheme sweep on the TC2 Ma=20 case: MUSCL vs MUSCL+THINC (int_comp).

Isolates THINC interface compression (`int_comp`, the only active anti-diffusion knob in MFC)
on the one case where it can engage. THINC compresses ONLY the volume fraction (eqn_idx%adv)
and partial densities (eqn_idx%cont) of MIXED interface cells (ic_eps <= alpha <= 1-ic_eps,
monotone profile); it NEVER touches the color function that drives the surface-tension CSF.
TC1 (case_Ma_0p1, identical fluids) has alpha1~=1 everywhere -> THINC is a structural no-op.
TC2 (case_Ma_20, distinct fluids, props 0.5x bulk) has a genuine 0->1 volume-fraction split
(~1300 THINC-eligible mixed cells), so THINC actually fires. This sweep asks: does int_comp=T
measurably change TC2 vs int_comp=F -- in the volume-fraction band, the color band, and U*(t*)?

To isolate int_comp, BOTH variants run MUSCL (recon_type=2, muscl_order=2, muscl_lim=4); the
ONLY difference is int_comp F vs T. The committed case ships WENO7, so each variant drops the
WENO-only keys (case_validator forbids them with MUSCL: weno_order must be unset/0; weno_eps,
mapped_weno, weno_avg, weno_Re_flux, etc. are prohibited) and adds the MUSCL keys. Patching is
done on the fully-built `data` dict (we replace the final print) so no key is missed.

    muscl        recon_type=2, muscl_order=2, muscl_lim=4, int_comp=F
    muscl_thinc  recon_type=2, muscl_order=2, muscl_lim=4, int_comp=T, ic_eps=1e-4, ic_beta=1.6
                 (1e-4 / 1.6 are the code defaults dflt_ic_eps / dflt_ic_beta; an older run used
                 ic_beta=2.0 -- 1.6 here.)

recon_type/int_comp/ic_* are runtime namelist params (NOT part of the compiled analytic IC, whose
rho/pres/cf strings are recon-independent), so both variants share ONE build: the first compiles,
the second reuses it. `./mfc.sh run` handles the build; do NOT pass --no-build. After EACH run we
GREP the generated simulation.inp and CONFIRM recon_type=2 and the intended int_comp -- an earlier
recon run silently emitted an inp with NO recon_type (ran WENO despite asking for MUSCL); this guard
makes that impossible to miss. Runs on the single Tesla V100 (--gpu acc --no-debug, -n 1), sequential.

case_Ma_20 has a known NaN risk (6-eq pressure-relaxation 0/0 when the bulk phase drains inside the
drop; fixed in src via max(rho_K_s,sgm_eps)). MUSCL changes the dynamics, so watch MFC.out / the tail.

Prereq: source the NVHPC toolchain first (this box has no modules; GPU is the only backend here):
    source ../../.nighthawk_gpu_env.sh && python3 sweep_recon.py

Usage:  python3 sweep_recon.py                  # both variants (muscl, muscl_thinc), sequential
        python3 sweep_recon.py muscl             # just the int_comp=F baseline
        python3 sweep_recon.py muscl_thinc       # just the THINC variant
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASE = "case_Ma_20.py"

# Don't let prterun pin to cores (single V100, single-user 8-core/16-thread Ryzen; run.py's
# "16-255" taskset range is invalid here). NOBIND keeps prterun from pinning all ranks to one core.
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}

# WENO-only keys the committed case sets that MUSCL (recon_type=2) forbids -- popped from `data`.
WENO_KEYS = ("weno_order", "weno_eps", "mapped_weno", "null_weights", "mp_weno", "weno_avg", "weno_Re_flux")

# variant -> the int_comp / ic_* overrides layered on the common MUSCL keys.
VARIANTS = {
    "muscl": {"int_comp": '"F"'},
    "muscl_thinc": {"int_comp": '"T"', "ic_eps": "1e-4", "ic_beta": "1.6"},
}


def make_case(dst, variant):
    """Copy the committed case into dst and replace its final print with a patch block that
    drops the WENO keys, sets the common MUSCL keys, and applies the variant's int_comp/ic_*.
    Operating on the fully-built `data` dict (not line regex) guarantees no WENO key is left
    behind to trip the validator and no MUSCL key is missed."""
    text = open(os.path.join(HERE, CASE)).read()
    overrides = VARIANTS[variant]
    lines = ["", "# --- sweep_recon.py: WENO7 -> MUSCL (recon_type=2); isolate int_comp ---"]
    lines.append(f"for _k in {WENO_KEYS!r}:")
    lines.append("    data.pop(_k, None)")
    lines.append('data["recon_type"] = 2')
    lines.append('data["muscl_order"] = 2')
    lines.append('data["muscl_lim"] = 4')
    for k, v in overrides.items():
        lines.append(f'data["{k}"] = {v}')
    lines.append("print(json.dumps(data))")
    patch = "\n".join(lines) + "\n"
    text, n = re.subn(r"(?m)^print\(json\.dumps\(data\)\)\s*$", patch, text, count=1)
    if n != 1:
        sys.exit(f"  PATCH FAILED for {dst}: replaced {n} print() lines (expected 1)")
    open(dst, "w").write(text)


def verify_inp(wd, variant):
    """CRITICAL guard: confirm the generated simulation.inp actually contains recon_type=2 and the
    intended int_comp. If recon_type is missing the run silently fell back to WENO -- abort."""
    inp = os.path.join(wd, "simulation.inp")
    if not os.path.isfile(inp):
        print(f"  VERIFY FAILED: no simulation.inp in {wd}")
        return False
    kv = {}
    for line in open(inp):
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip().lower()] = v.strip().rstrip(",")
    recon = kv.get("recon_type")
    int_comp = kv.get("int_comp")
    want_int = "T" if variant == "muscl_thinc" else "F"
    ok = recon == "2" and int_comp == want_int
    detail = f"recon_type={recon!r} int_comp={int_comp!r}"
    if "ic_beta" in kv or "ic_eps" in kv:
        detail += f" ic_eps={kv.get('ic_eps')!r} ic_beta={kv.get('ic_beta')!r}"
    if "weno_order" in kv:
        detail += f"  weno_order={kv['weno_order']!r}(!)"
    print(f"  inp check [{variant}]: {detail}  -> {'OK' if ok else 'MISMATCH'}")
    if not ok:
        print("  STOP: simulation.inp does not match the intended MUSCL/int_comp config; do not trust results.")
    return ok


def run_variant(variant):
    wd = os.path.join(HERE, "runs", "recon_tc2", variant)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    dst = os.path.join(wd, CASE)
    make_case(dst, variant)

    rel = os.path.relpath(dst, REPO)
    # -j 8: the first variant's analytic IC build runs while the GPU is idle; the second reuses it
    # (recon/int_comp are runtime-only). Explicit --gpu acc so the run does not inherit a stale CPU
    # lock.yaml default. NO --no-build (the analytic IC still needs the first build).
    cmd = ["./mfc.sh", "run", rel, "--no-debug", "-j", "8", "-n", "1", "--gpu", "acc"]
    print(f"\n>>> {variant}:  recon_type=2 int_comp={'T' if variant == 'muscl_thinc' else 'F'}  -> {rel}", flush=True)
    p = subprocess.run(cmd, cwd=REPO, env={**os.environ, **NOBIND}, capture_output=True, text=True, check=False)
    log = os.path.join(wd, "sweep_recon.log")
    open(log, "w").write(p.stdout + "\n===STDERR===\n" + p.stderr)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}); tail of output:")
        print("\n".join((p.stdout + p.stderr).splitlines()[-25:]))
        verify_inp(wd, variant)
        return False
    # Echo the final progress line so the step rate / wall time is visible, and flag any NaN.
    tail = p.stdout.splitlines()
    if any("NaN" in ln or "nan" in ln for ln in tail[-40:]):
        print("  WARNING: 'NaN' appears in the run tail -- inspect MFC.out before trusting results.")
    for line in reversed(tail):
        if "Total-time" in line or "Time/step" in line or "sec" in line.lower():
            print(f"  {line.strip()}")
            break
    ok = verify_inp(wd, variant)
    print(f"  {'OK' if ok else 'INP MISMATCH'} -> {wd}")
    return ok


def main():
    args = sys.argv[1:]
    values = args or list(VARIANTS)
    bad = [v for v in values if v not in VARIANTS]
    if bad:
        sys.exit(f"unknown variant(s) {bad}; choose from {list(VARIANTS)}")
    print(f"recon sweep (TC2 Ma=20): {values}  (1 rank GPU each, MUSCL, isolate int_comp)")
    results = {v: run_variant(v) for v in values}
    print("\n=== recon sweep run summary ===")
    for v, ok in results.items():
        print(f"  {v:<12}: {'OK' if ok else 'FAILED / INP MISMATCH'}")
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
