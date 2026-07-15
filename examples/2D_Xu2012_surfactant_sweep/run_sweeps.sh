#!/usr/bin/env bash
# Xu 2012 M2 property sweeps: vary one control group at a time about the baseline (Ca=0.3, Pe=10,
# lambda=1, Re=1, X=0.1) and record quasi-steady D, theta, mass, surfactant non-uniformity. All are
# namelist params -> one build serves every point (gdot=1 shear IC unchanged). Baseline (ca=0.3) is
# shared by all four panels. weno5 needs NX=128, NP=8. Run from repo root. ~8 runs, few min each.
set -e
cd "$(git rev-parse --show-toplevel)"
EX=examples/2D_Xu2012_surfactant_sweep
: >"$EX/results.jsonl"

run() { # group  swept_value  ENV=...
    local G=$1 XV=$2
    shift 2
    echo ">>> $G=$XV  ($*)"
    env MFC_NX=128 MFC_SURF=0.1 "$@" ./mfc.sh run "$EX/case.py" --no-debug -n 8 -t pre_process simulation >/dev/null
    local RES
    RES=$(python3 "$EX/measure_m2.py" "$EX")
    python3 -c "import json;d=json.loads('''$RES''');d.update(group='$G',x=$XV);print(json.dumps(d))" >>"$EX/results.jsonl"
    tail -1 "$EX/results.jsonl"
}

# Ca sweep (gate: D increases with Ca)
run ca 0.2 MFC_CA=0.2
run ca 0.3 MFC_CA=0.3   # BASELINE (shared by all panels)
run ca 0.4 MFC_CA=0.4
# Viscosity-ratio sweep (gate: D decreases as drop gets more viscous)
run lam 0.5 MFC_LAMBDA=0.5
run lam 2.0 MFC_LAMBDA=2
# Reynolds sweep (gate: D increases with Re at low Re)
run re 2.0 MFC_RE=2
# Peclet sweep (measure the surfactant-distribution + D response)
run pe 1.0 MFC_PE=1
run pe 100.0 MFC_PE=100
echo "done -> $EX/results.jsonl"
