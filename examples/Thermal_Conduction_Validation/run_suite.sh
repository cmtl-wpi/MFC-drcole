#!/usr/bin/env bash
# Run the energy-coupled conduction cases and archive each run's restart_data +
# simulation.inp into runs/<label>/, which is what validate.py reads. The cheap
# 1D convergence sweeps run single-rank; the costly 2D/3D benchmarks run
# multi-rank. All runs stage in cases/, so this is strictly sequential.
#
#   bash examples/Thermal_Conduction_Validation/run_suite.sh          # everything
#   bash examples/Thermal_Conduction_Validation/run_suite.sh conv     # only the conv sweeps
set -uo pipefail
EX="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"  # this example's directory
ROOT="$(cd "$EX/../.." && pwd)"                      # repo root (holds mfc.sh)
cd "$ROOT"
WHAT=${1:-all}

run_and_archive() { # <case.py> <label> <nranks>
    local case=$1 label=$2 nr=${3:-1}
    echo ">>> $label ($case, n=$nr)"
    ./mfc.sh run "$EX/cases/$case" -t pre_process simulation -n "$nr" || {
        echo "FAILED $label"
        return 1
    }
    rm -rf "$EX/runs/$label"
    mkdir -p "$EX/runs/$label"
    cp -r "$EX/cases/restart_data" "$EX/runs/$label/"
    cp "$EX/cases/simulation.inp" "$EX/runs/$label/"
    echo "<<< archived $label"
}

# Spatial convergence: vary grid, dt capped at the acoustic CFL inside the case.
for N in 32 64 128 256 512; do
    CONV_N=$N run_and_archive case_conv.py "convx_$N" 1
done

# Temporal convergence: fixed N=32, vary dt to a common final time t*.
TSTAR=$(python3 -c "import math; print(0.3/(0.05*(2*math.pi)**2))")
for NS in 256 512 1024 2048 4096; do
    DT=$(python3 -c "print($TSTAR/$NS)")
    CONV_N=32 CONV_DT=$DT CONV_NSTEPS=$NS run_and_archive case_conv.py "convt_$NS" 1
done

if [ "$WHAT" = all ]; then
    run_and_archive case_2d_mode.py 2d_mode 8
    run_and_archive case_3d_hotspot.py 3d_hotspot 8
fi

echo "ALL DONE"
