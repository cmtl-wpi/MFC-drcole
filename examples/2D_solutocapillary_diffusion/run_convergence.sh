#!/usr/bin/env bash
# 2D circle surface-diffusion convergence sweep: record the m=1 rate vs R/dx across resolutions.
# The x/R IC is resolution-independent, so only the first NX builds. Release build, multi-rank.
#   examples/2D_solutocapillary_diffusion/run_convergence.sh
cd "$(git rev-parse --show-toplevel)" || exit 1
EX=examples/2D_solutocapillary_diffusion
OUT="$EX/convergence.dat"
echo "# R/dx  full-field  band-only  exact(D_s/R^2)  drift%" > "$OUT"
for NX in 48 96 192 256; do
    case "$NX" in 96) NP=4 ;; 192) NP=16 ;; 256) NP=16 ;; *) NP=1 ;; esac
    echo ">>> NX=$NX NP=$NP $(date +%H:%M:%S)"
    MFC_NX="$NX" ./mfc.sh run "$EX/case.py" --no-debug -n "$NP" -t pre_process simulation || { echo "run $NX failed"; exit 1; }
    python3 "$EX/measure.py" >> "$OUT" || { echo "measure $NX failed"; exit 1; }
    echo "=== after NX=$NX ==="; cat "$OUT"
done
python3 "$EX/plot_convergence.py"
