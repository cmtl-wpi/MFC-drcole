#!/usr/bin/env bash
# Run one AMR+ST experiment case into a canonical run dir, with the Intel toolchain
# and the pinned worktree binary, then stamp a provenance manifest.
#
#   run_case.sh <case_file> <run_id> <ranks> [-- <case args...>]
#
# <case_file> is relative to cases/ (e.g. case_laplace.py). Output lands in
# runs/<run_id>/. Case args after `--` are forwarded to the case script by mfc.sh.
set -uo pipefail

EXP=/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st
WORKTREE=/home/daveygravy/repos/MFC-amr-st
PINNED_SHA=ace2285a7e72fdce9dc3fbc3a5629e1b9d1a89b7

CASE_FILE="$1"; RUN_ID="$2"; RANKS="$3"; shift 3
CASE_ARGS=("$@")   # includes a leading `--` if present

RUN_DIR="$EXP/runs/$RUN_ID"
mkdir -p "$RUN_DIR"
cp "$EXP/cases/$CASE_FILE" "$RUN_DIR/case.py"

source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export FC=ifx CC=icx CXX=icpx

cd "$WORKTREE"
t0=$(date +%s)
echo "[run] $RUN_ID  ranks=$RANKS  case=$CASE_FILE  args=${CASE_ARGS[*]}"
echo "[run] SHA=$(git rev-parse HEAD)  (pinned $PINNED_SHA)"
./mfc.sh run "$RUN_DIR/case.py" -n "$RANKS" "${CASE_ARGS[@]}"
rc=$?
t1=$(date +%s)
wall=$((t1 - t0))
echo "[run] exit=$rc wall=${wall}s"

# ---- provenance manifest ------------------------------------------------------
last_step=$(ls "$RUN_DIR"/restart_data/lustre_*.dat 2>/dev/null \
            | sed -E 's/.*lustre_([0-9]+)\.dat/\1/' | sort -n | tail -1)
slug=$(ls -td "$WORKTREE"/build/install/cpu-* 2>/dev/null | head -1 | xargs -n1 basename)
python3 - "$RUN_DIR" "$RUN_ID" "$RANKS" "$rc" "$wall" "${last_step:-none}" "$slug" "$PINNED_SHA" <<'PY'
import json, sys, socket, os
run_dir, run_id, ranks, rc, wall, last_step, slug, sha = sys.argv[1:9]
m = {
    "run_id": run_id,
    "git_commit": sha,
    "git_branch": "amr-st-experiment (worktree @ PR1628 sbryngelson:up/mega)",
    "invocation": " ".join(sys.argv[0:1]) + f"  (mfc.sh run, -n {ranks})",
    "build_slug": slug,
    "host": socket.gethostname(),
    "n_ranks": int(ranks),
    "mpi_impl": "openmpi (ifx/icx toolchain)",
    "final_step_reached": last_step,
    "wall_time_s": int(wall),
    "exit_status": int(rc),
}
with open(os.path.join(run_dir, "run_manifest.json"), "w") as f:
    json.dump(m, f, indent=2)
print("[run] manifest ->", os.path.join(run_dir, "run_manifest.json"))
PY
exit $rc
