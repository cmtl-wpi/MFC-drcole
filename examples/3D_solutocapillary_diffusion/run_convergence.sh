#!/usr/bin/env bash
# Sphere surface-diffusion convergence sweep: run the l=1 case at several resolutions and record the
# measured decay rate vs R/dx. The analytic IC (z/R) is resolution-independent, so only the first NX
# builds; the rest reuse it. Release build (--no-debug); multi-rank for the finer grids.
#   examples/3D_solutocapillary_diffusion/run_convergence.sh
cd "$(git rev-parse --show-toplevel)" || exit 1
EX=examples/3D_solutocapillary_diffusion
OUT="$EX/convergence.dat"
echo "# R/dx  measured  exact(2 D_s/R^2)  drift%" > "$OUT"
for NX in 32 48 64 96; do
    case "$NX" in 64) NP=8 ;; 96) NP=27 ;; *) NP=1 ;; esac
    echo ">>> NX=$NX  NP=$NP"
    MFC_NX="$NX" ./mfc.sh run "$EX/case.py" --no-debug -n "$NP" -t pre_process simulation || { echo "run NX=$NX failed"; exit 1; }
    python3 "$EX/measure.py" >> "$OUT" || { echo "measure NX=$NX failed"; exit 1; }
done
python3 "$EX/plot_convergence.py"
echo "=== convergence.dat ==="; cat "$OUT"
