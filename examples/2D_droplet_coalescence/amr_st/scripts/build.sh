#!/usr/bin/env bash
# Build the AMR+ST experiment binaries from the pinned worktree (PR 1628 @ ace2285a)
# with the Intel oneAPI toolchain (ifx ~2x faster than gfortran on this EPYC).
#
# Intel-env gotchas baked in:
#   - setvars.sh returns rc 3 on re-init; sourcing under `set -e` aborts silently -> `|| true`.
#   - MFC's CMake honors FC/CC/CXX; point them at the Intel compilers explicitly.
set -uo pipefail

WORKTREE=/home/daveygravy/repos/MFC-amr-st
JOBS="${1:-32}"

source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export FC=ifx CC=icx CXX=icpx

cd "$WORKTREE"
echo "[build] worktree HEAD: $(git rev-parse HEAD)"
echo "[build] FC=$(command -v ifx)  jobs=$JOBS"
./mfc.sh build -t pre_process simulation post_process -j "$JOBS"
rc=$?
echo "[build] exit=$rc"
exit $rc
