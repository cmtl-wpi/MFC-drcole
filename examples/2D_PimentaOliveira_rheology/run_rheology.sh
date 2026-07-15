#!/usr/bin/env bash
# Pimenta & Oliveira M3 rheology: run coverage X = 0, 0.1, 0.3 at Ca=0.1 and extract the interfacial
# stress decomposition ([eta_c], [eta_m], N1). Namelist-only sweep -> one build serves all. weno5 needs
# NX=128, NP=8. Run from repo root. Gate: X=0 -> [eta_m]~0; higher X -> [eta_m] grows; N1>0.
cd "$(git rev-parse --show-toplevel)"
EX=examples/2D_PimentaOliveira_rheology
: >"$EX/results.jsonl"
for X in 0 0.1 0.3; do
    echo ">>> coverage X=$X"
    # tolerate a SIGTERM (143) at toolchain cleanup: the simulation completes and writes all frames first.
    MFC_NX=128 MFC_SURF="$X" ./mfc.sh run "$EX/case.py" --no-debug -n 8 -t pre_process simulation >/dev/null || true
    RES=$(python3 "$EX/measure_m3.py" "$EX")
    python3 -c "import json;d=json.loads('''$RES''');d.update(X=$X);print(json.dumps(d))" >>"$EX/results.jsonl"
    tail -1 "$EX/results.jsonl"
done
echo "done -> $EX/results.jsonl"
