#!/usr/bin/env bash
# Xu 2006 coverage sweep: run X = 0, 0.1, 0.3 (constant surf_val = namelist, so one build serves all)
# and print Taylor D + inclination theta + surfactant-mass conservation vs time. Higher coverage ->
# lower interfacial sigma (Langmuir) -> more deformation. Run from repo root. weno5 needs NX=128, NP=8.
cd "$(git rev-parse --show-toplevel)" || exit 1
EX=examples/2D_Xu2006_surfactant_shear
for X in 0 0.1 0.3; do
    echo ">>> coverage X=$X"
    MFC_NX=128 MFC_SURF="$X" ./mfc.sh run "$EX/case.py" --no-debug -n 8 -t pre_process simulation || exit 1
    python3 "$EX/measure.py" "$EX" | tail -3
done
