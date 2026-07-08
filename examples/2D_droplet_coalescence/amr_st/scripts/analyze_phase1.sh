#!/usr/bin/env bash
# Run the full Phase-1 analysis: per-run metrics (coarse, coarse-in-band control, fine,
# amr) + the seam-current verdict. Safe to re-run as runs accumulate checkpoints.
set -u
A=/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st/analysis
R=/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st/runs

# AMR block cell indices (for the coarse same-band control)
inp="$R/laplace__amr/simulation.inp"
bx0=$(grep -E "amr_block_beg\(1\)" "$inp" | grep -oE "[0-9]+$")
bx1=$(grep -E "amr_block_end\(1\)" "$inp" | grep -oE "[0-9]+$")
by0=$(grep -E "amr_block_beg\(2\)" "$inp" | grep -oE "[0-9]+$")
by1=$(grep -E "amr_block_end\(2\)" "$inp" | grep -oE "[0-9]+$")

python3 "$A/seam_analysis.py" "$R/laplace__coarse" --label coarse      2>/dev/null | grep RESULT_JSON >/dev/null
python3 "$A/seam_analysis.py" "$R/laplace__coarse" --label coarse_band --band-block "$bx0,$bx1,$by0,$by1" 2>/dev/null | grep RESULT_JSON >/dev/null
[ -d "$R/laplace__fine/restart_data" ] && python3 "$A/seam_analysis.py" "$R/laplace__fine" --label fine 2>/dev/null | grep RESULT_JSON >/dev/null
python3 "$A/seam_analysis.py" "$R/laplace__amr" --label amr         2>/dev/null | grep RESULT_JSON >/dev/null
python3 "$A/compare_phase1.py"
