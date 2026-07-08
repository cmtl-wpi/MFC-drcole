#!/usr/bin/env bash
# Durable self-healing run for the Phase-1 seam test (WF-RUN-SHARED). Launched inside a
# detached screen so it lives outside the Claude session process tree (the reaper kills
# detached session-child processes on long idle). MFC's deterministic restart makes a
# resumed trajectory identical to an uninterrupted one, so crash/reap recovery is lossless.
#
#   run_durable.sh <run_id> <ranks> <n_saves> <cpuset>
#
# The run dir's case.py must already exist (write it with gen_case.py first). Only 4 simple
# positional args cross the screen boundary -- no case-arg list to be re-tokenized. On a
# crash/reap the loop re-runs from the latest checkpoint (rewriting the case's _NSTART),
# aborting after 3 no-progress rounds.
set -uo pipefail
EXP=/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st
WORKTREE=/home/daveygravy/repos/MFC-amr-st
PINNED_SHA=ace2285a7e72fdce9dc3fbc3a5629e1b9d1a89b7

RUN_ID="$1"; RANKS="$2"; NSAVES="$3"; CPUSET="$4"
RUN_DIR="$EXP/runs/$RUN_ID"
log="$RUN_DIR/durable.log"

source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export FC=ifx CC=icx CXX=icpx

[ -f "$RUN_DIR/case.py" ] || { echo "[durable] MISSING $RUN_DIR/case.py (run gen_case.py)" >>"$log"; exit 1; }
echo "[durable $(date '+%F %T')] start $RUN_ID ranks=$RANKS cpuset=$CPUSET" >>"$log"

latest() { ls "$RUN_DIR"/restart_data/lustre_[0-9]*.dat 2>/dev/null \
           | sed -E 's/.*lustre_([0-9]+)\.dat/\1/' | sort -n | tail -1; }

cd "$WORKTREE"
stall=0; t0=$(date +%s)
while :; do
    m=$(latest); m=${m:-0}
    if [ "$m" -ge "$NSAVES" ]; then echo "[durable $(date '+%T')] DONE at save $m" >>"$log"; break; fi
    echo "[durable $(date '+%T')] run from save $m / $NSAVES" >>"$log"
    sed -i -E "s/^_NSTART = .*/_NSTART = $m/" "$RUN_DIR/case.py"        # bake restart index
    taskset -c "$CPUSET" ./mfc.sh run "$RUN_DIR/case.py" -n "$RANKS" \
        -t pre_process simulation >>"$log" 2>&1 || true
    m2=$(latest); m2=${m2:-0}
    if [ "$m2" -le "$m" ]; then
        stall=$((stall + 1))
        echo "[durable $(date '+%T')] no progress (stall $stall) at save $m2" >>"$log"
        [ "$stall" -ge 3 ] && { echo "[durable] ABORT: 3 stalls" >>"$log"; break; }
        sleep 5
    else
        stall=0
    fi
done

wall=$(( $(date +%s) - t0 )); fin=$(latest)
slug=$(ls -td "$WORKTREE"/build/install/cpu-*/bin/simulation | head -1 | xargs -n1 dirname | xargs -n1 dirname | xargs -n1 basename)
python3 - "$RUN_DIR" "$RUN_ID" "$RANKS" "$wall" "${fin:-none}" "$slug" "$PINNED_SHA" <<'PY'
import json, sys, socket, os
rd, rid, ranks, wall, fin, slug, sha = sys.argv[1:8]
json.dump({"run_id": rid, "git_commit": sha,
          "git_branch": "amr-st-experiment (PR1628 @ ace2285a)",
          "invocation": f"mfc.sh run -n {ranks} -t pre_process simulation (standalone case.py)",
          "build_slug": slug, "host": socket.gethostname(), "n_ranks": int(ranks),
          "mpi_impl": "intel-mpi (ifx/icx)", "final_save_reached": fin,
          "wall_time_s": int(wall)},
         open(os.path.join(rd, "run_manifest.json"), "w"), indent=2)
PY
echo "[durable $(date '+%F %T')] exit; manifest written" >>"$log"
