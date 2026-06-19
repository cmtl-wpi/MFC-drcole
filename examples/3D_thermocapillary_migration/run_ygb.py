#!/usr/bin/env python3
"""Run the 3D thermocapillary u_YGB validation sweep, measure each run, collect the results.

Each run is one leaf dir runs/ygb/<geom>/<W>/<grid>/<Ma>/ (MFC writes its output next to the
case copy), selected by the env vars YGB_GEOM/YGB_W/YGB_NX/YGB_MA/YGB_TR that case_ygb.py reads.

Unlike the 2D run.py, this driver NEVER deletes a populated leaf: if a leaf already has
restart_data/ it is SKIPPED (pass --force to rerun it). A partial sweep is therefore safe to
re-invoke -- you won't clobber checkpoints from a run that's still going or already done.

Selectors (around the converged corner cube / W=8 / Nx=80 / Ma=0.5):
  smoke        one tiny cube run -- machinery check (serial; minutes for a short t_r)
  anchor       Samareh-box geometry -- reproduce the confined v_t/v_YGB ~ 0.95 with this physics
  confinement  cube, W in {6,8,10,12} -- the headline series; extrapolate 1/W -> 0
  grid         cube, Nx in {64,96,128} at W=10 -- Richardson refinement in dx
  ma           cube, Ma in {1.0,0.5,0.25,0.1} at W=10 -- finite-Ma extrapolation to Ma -> 0
  all          anchor + confinement + grid + ma

Usage:
    python3 run_ygb.py [smoke|anchor|confinement|grid|ma|all] [--force]   (default: smoke)

The heavy selectors are a multi-hour-per-run, multi-day campaign -- run under nohup / in the
background. Build happens once up front (the grid/W/Ma don't change the binary).
"""

import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs", "ygb")
RESULTS = os.path.join(HERE, "results")
CASE = "case_ygb.py"

PIN = ["taskset", "-c", "16-255"]  # keep off cores 0-15 (a neighbour's job may live there)
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}  # don't let prterun pin to cores

# Each run is (geom, W_in_D, Nx, Ma, n_tr). The leaf path is derived from (geom, W, Nx, Ma).
SWEEPS = {
    "smoke": [("cube", 5, 40, 1.0, 0.3)],
    "anchor": [("samareh", 5, 64, 1.0, 2.0)],
    "confinement": [("cube", w, 80, 0.5, 3.0) for w in (6, 8, 10, 12)],
    "grid": [("cube", 10, nx, 0.5, 3.0) for nx in (64, 96, 128)],
    "ma": [("cube", 10, 80, ma, 3.0) for ma in (1.0, 0.5, 0.25, 0.1)],
}
SWEEPS["all"] = SWEEPS["anchor"] + SWEEPS["confinement"] + SWEEPS["grid"] + SWEEPS["ma"]


def fmt(value):
    """Encode a number for a run-dir name: decimal point -> 'p' (e.g. 7.5 -> '7p5', 0.5 -> '0p5')."""
    return ("%g" % value).replace(".", "p")


def leaf_dir(geom, W, Nx, Ma):
    """runs/ygb/<geom>/w<W>/nx<Nx>/ma<Ma>/ -- the encoding IS the source of truth for the sweep coords."""
    return os.path.join(RUNS, geom, f"w{fmt(W)}", f"nx{Nx:03d}", f"ma{fmt(Ma)}")


def ranks_for(Nx):
    """MPI ranks for an Nx-per-width grid: cube the per-dim split that keeps blocks >= 32 cells.

    weno5 needs a comfortable rank-block (>= ~25 cells) in every split dimension; >=32 leaves margin.
    Nx=40 -> 1, 64 -> 8, 80 -> 8, 96 -> 27, 128 -> 64. MFC factors -n into a valid decomposition or
    aborts, so passing a cube whose factorization exists guarantees a legal split.
    """
    per_dim = max(1, Nx // 32)
    return per_dim**3


def run_one(geom, W, Nx, Ma, n_tr, force):
    """Run one variant into its leaf dir. Skip a populated leaf unless force. Returns True on success."""
    wd = leaf_dir(geom, W, Nx, Ma)
    if os.path.isdir(os.path.join(wd, "restart_data")) and not force:
        print(f">>> SKIP (already populated): {os.path.relpath(wd, HERE)}  -- pass --force to rerun", flush=True)
        return True
    if force and os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd, exist_ok=True)
    shutil.copy(os.path.join(HERE, CASE), os.path.join(wd, CASE))

    ranks = ranks_for(Nx)
    env = {**os.environ, **NOBIND, "YGB_GEOM": geom, "YGB_W": str(W), "YGB_NX": str(Nx), "YGB_MA": str(Ma), "YGB_TR": str(n_tr)}
    rel = os.path.relpath(os.path.join(wd, CASE), REPO)
    print(f"\n>>> {os.path.relpath(wd, HERE)}  geom={geom} W={W} Nx={Nx} Ma={Ma} t_r={n_tr} ranks={ranks}", flush=True)
    p = subprocess.run(PIN + ["./mfc.sh", "run", rel, "-n", str(ranks)], cwd=REPO, env=env, check=False)
    if p.returncode != 0:
        print(f"  RUN FAILED (exit {p.returncode}) -- see {os.path.join(wd, 'MFC.out')}")
        return False
    return True


def measure(wd):
    """Run measure.py <wd> and return its RESULT_JSON dict (or None on failure)."""
    m = subprocess.run([sys.executable, os.path.join(HERE, "measure.py"), wd], capture_output=True, text=True, check=False)
    if m.returncode != 0:
        print(f"  MEASURE FAILED: {m.stderr[-400:]}")
        return None
    for line in m.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


def main():
    args = sys.argv[1:]
    force = "--force" in args
    sel = next((a for a in args if not a.startswith("-")), "smoke")
    if sel not in SWEEPS:
        sys.exit(f"unknown selector {sel!r}; choose from {', '.join(SWEEPS)}")

    os.makedirs(RESULTS, exist_ok=True)
    spath = os.path.join(RESULTS, "ygb_summary.json")
    summary = json.load(open(spath)) if os.path.isfile(spath) else {}

    for geom, W, Nx, Ma, n_tr in SWEEPS[sel]:
        wd = leaf_dir(geom, W, Nx, Ma)
        if not run_one(geom, W, Nx, Ma, n_tr, force):
            continue
        res = measure(wd)
        if res is not None:
            key = os.path.relpath(wd, RUNS)  # e.g. cube/w8/nx080/ma0p5
            summary[key] = res
            json.dump(summary, open(spath, "w"), indent=2)  # checkpoint after each run
            print(f"  {key}: v_t/v_YGB plateau={res['ratio_plateau']:+.3f}  W={res['W']:g}  cells/D={res['cells_per_D']:.1f}  rises={res['rises']}")

    print(f"\n=== {sel}: {len(summary)} runs in {spath} ===")


if __name__ == "__main__":
    main()
