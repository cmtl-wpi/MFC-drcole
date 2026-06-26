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
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs", "ygb")
RESULTS = os.path.join(HERE, "results")
CASE = "case_ygb.py"

PIN = ["taskset", "-c", "16-255"]  # keep off cores 0-15 (a neighbour's job may live there)
NOBIND = {"OMPI_MCA_hwloc_base_binding_policy": "none"}  # don't let prterun pin to cores

# Each run is (geom, W_in_D, Nx, Ma, n_tr). The leaf path is derived from (geom, W, Nx, Ma).
# The confinement sweep holds the DROP RESOLUTION fixed (cells_per_D = Nx/W = 8) by scaling Nx with
# the box width, so it isolates confinement instead of confounding it with grid coarsening. The grid
# sweep then refines cells_per_D at the fixed W=10 corner; the ma sweep varies Ma there.
SWEEPS = {
    "smoke": [("cube", 5, 40, 1.0, 0.3)],
    "anchor": [("samareh", 5, 64, 1.0, 2.0)],
    "confinement": [("cube", w, 8 * w, 0.5, 3.0) for w in (8, 10, 12, 14)],  # Nx = 8*W -> 8 cells/D
    "grid": [("cube", 10, nx, 0.5, 3.0) for nx in (64, 96, 128)],  # + the W10/Nx80 corner = 4 grids
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
    """MPI ranks for an Nx-per-width grid: cube the per-dim split that keeps blocks >= ~26 cells.

    weno5 needs >= ~25 cells per rank-block in every split dimension (verified: Nx=80 -> 27 ranks,
    ~26.7 cells/block, decomposes fine). Nx=40 -> 1, 64 -> 8, 80 -> 27, 96 -> 27, 128 -> 64. MFC
    factors -n into a valid decomposition or aborts, so a cube whose factorization exists is legal.
    """
    per_dim = max(1, Nx // 26)
    return per_dim**3


def launch(geom, W, Nx, Ma, n_tr, force):
    """Start one variant as a background process. Return (popen, wd, ranks, key), or None if skipped.

    Skips a populated leaf unless force. Uses --case-optimization --no-build: the simulation binary
    is grid- and Ma-independent (m/n/p and k_therm are runtime, not in CASE_OPT_PARAMS) so ONE sim
    build serves the whole sweep, but the analytic IC bakes a dx-dependent interface width
    (w_if = 0.75*dx) into the patch strings, so each distinct grid needs its own pre_process build.
    Pre-build every grid (case-optimized) before sweeping; --no-build then finds each by slug.
    Per-leaf stdout -> <leaf>/run.log so concurrent runs don't mix.
    """
    wd = leaf_dir(geom, W, Nx, Ma)
    key = os.path.relpath(wd, RUNS)
    if os.path.isdir(os.path.join(wd, "restart_data")) and not force:
        print(f">>> SKIP (already populated): {key}  -- pass --force to rerun", flush=True)
        return None
    if force and os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd, exist_ok=True)
    shutil.copy(os.path.join(HERE, CASE), os.path.join(wd, CASE))

    ranks = ranks_for(Nx)
    env = {**os.environ, **NOBIND, "YGB_GEOM": geom, "YGB_W": str(W), "YGB_NX": str(Nx), "YGB_MA": str(Ma), "YGB_TR": str(n_tr)}
    rel = os.path.relpath(os.path.join(wd, CASE), REPO)
    log = open(os.path.join(wd, "run.log"), "w")
    print(f">>> LAUNCH {key}  Nx={Nx} Ma={Ma} t_r={n_tr} ranks={ranks}", flush=True)
    p = subprocess.Popen(
        PIN + ["./mfc.sh", "run", rel, "-n", str(ranks), "-t", "pre_process", "simulation", "--mpi", "--case-optimization", "--no-build"], cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT
    )
    return p, wd, ranks, key


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
    # Run several leaves at once up to a total-rank (core) budget so the sweep saturates the box
    # instead of idling cores. Each run's ranks come from ranks_for(Nx). Default leaves headroom on
    # a 256-core host (cores 0-15 are pinned off via PIN).
    maxcores = int(os.environ.get("YGB_MAXCORES", "192"))

    os.makedirs(RESULTS, exist_ok=True)
    spath = os.path.join(RESULTS, "ygb_summary.json")
    summary = json.load(open(spath)) if os.path.isfile(spath) else {}

    queue = list(SWEEPS[sel])
    inflight = {}  # popen -> (wd, ranks, key)
    used = 0
    while queue or inflight:
        # Fill the core budget. Always allow one run if the box is idle (a single run < maxcores).
        while queue and (not inflight or used + ranks_for(queue[0][2]) <= maxcores):
            handle = launch(*queue.pop(0), force=force)
            if handle is None:
                continue  # skipped (already populated)
            p, wd, ranks, key = handle
            inflight[p] = (wd, ranks, key)
            used += ranks
        if not inflight:
            continue
        # Wait for at least one in-flight run to finish, then measure it and free its cores.
        done = []
        while not done:
            time.sleep(5)
            done = [p for p in inflight if p.poll() is not None]
        for p in done:
            wd, ranks, key = inflight.pop(p)
            used -= ranks
            if p.returncode != 0:
                print(f"  RUN FAILED ({key}, exit {p.returncode}) -- see {os.path.join(wd, 'run.log')}", flush=True)
                continue
            res = measure(wd)
            if res is not None:
                summary[key] = res
                json.dump(summary, open(spath, "w"), indent=2)  # checkpoint after each run
                print(f"  DONE {key}: v_t/v_YGB plateau={res['ratio_plateau']:+.3f}  W={res['W']:g}  cells/D={res['cells_per_D']:.1f}  rises={res['rises']}  ({res['t_end_tr']:.1f} t_r)", flush=True)

    print(f"\n=== {sel}: {len(summary)} runs in {spath} ===", flush=True)


if __name__ == "__main__":
    main()
