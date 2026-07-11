#!/usr/bin/env bash
#
# Build and run the Basilisk translation of MFC 3D_droplet_coalescence case p.
#
#   ./run.sh build [LEVEL]           serial build            -> ./collision
#   ./run.sh mpi   [LEVEL]           MPI build               -> ./collision.mpi
#   ./run.sh run   NPROCS [LEVEL]    build MPI + run          (default LEVEL=10)
#
# LEVEL is the finest AMR level (dx_min = 6D / 2^LEVEL): 10 ~ D/170 (production),
# 7-8 for a quick smoke test.
set -euo pipefail

export BASILISK="${BASILISK:-$HOME/basilisk/src}"
export PATH="$BASILISK:$PATH"

cmd="${1:-run}"
LEVEL_DEFAULT=10

build_serial () {
  local lvl="$1"
  qcc -O2 -DMAXLEVEL="$lvl" collision.c -o collision -lm
  echo "built ./collision  (MAXLEVEL=$lvl)"
}

build_mpi () {
  local lvl="$1"
  # -D_GNU_SOURCE: Intel mpicc's strict mode otherwise hides madvise() from
  # <sys/mman.h>, which Basilisk's octree memory index needs.
  CC99='mpicc -std=c99 -D_GNU_SOURCE' \
    qcc -D_MPI=1 -O2 -DMAXLEVEL="$lvl" collision.c -o collision.mpi -lm
  echo "built ./collision.mpi  (MAXLEVEL=$lvl)"
}

case "$cmd" in
  build) build_serial "${2:-$LEVEL_DEFAULT}" ;;
  mpi)   build_mpi    "${2:-$LEVEL_DEFAULT}" ;;
  run)
    nproc="${2:?usage: ./run.sh run NPROCS [LEVEL]}"
    lvl="${3:-$LEVEL_DEFAULT}"
    build_mpi "$lvl"
    echo "running on $nproc ranks..."
    mpirun -np "$nproc" ./collision.mpi
    ;;
  *) echo "usage: ./run.sh {build|mpi|run} ..."; exit 1 ;;
esac
