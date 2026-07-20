#!/usr/bin/env bash
# Coverage sweep for the Marangoni-coupling deformation test: run surf_val=0,1,2 (namelist constant, so
# one build serves all three) and print Taylor D(t) for each. Higher coverage -> lower interfacial sigma
# -> more deformation. Run from repo root.
cd "$(git rev-parse --show-toplevel)" || exit 1
EX=examples/2D_solutocapillary_deformation
for S in 0 1 2; do
    echo ">>> coverage surf_val=$S"
    MFC_NX=128 MFC_SURF="$S" ./mfc.sh run "$EX/case.py" --no-debug -n 4 -t pre_process simulation || exit 1
    python3 "$EX/measure.py" "$EX" | tail -3
done
