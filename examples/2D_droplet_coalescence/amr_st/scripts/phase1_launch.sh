#!/usr/bin/env bash
# Launch the Phase-1 trio in detached screens (survive the session reaper).
# Reduced config (see README): c_l=12, +-3.2R, R/50, 2 tau, 100 saves.
#   coarse : 320^2 uniform         (32 ranks, socket-0 cores 8-39)
#   fine   : 640^2 uniform         (32 ranks, socket-1 cores 64-95)
#   amr    : 320^2 + static block  (np=1 -- forced by the multi-rank ST halo bug, core 0)
# Cases are pre-generated as standalone case.py (gen_case.py); run_durable takes only 4
# simple positional args so nothing complex crosses the screen boundary. Staggered so the
# three mfc.sh runs don't collide on the shared build/staging dir.
set -u
D=/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st/scripts
R=/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st/runs
C="--c-l 12 --domain-R 3.2 --n-periods 2 --n-saves 100"

python3 "$D/gen_case.py" laplace__coarse "$R/laplace__coarse" -- --variant coarse $C
python3 "$D/gen_case.py" laplace__fine   "$R/laplace__fine"   -- --variant fine   $C
python3 "$D/gen_case.py" laplace__amr    "$R/laplace__amr"    -- --variant amr --block-R 1.5 $C

screen -dmS amrst_coarse bash "$D/run_durable.sh" laplace__coarse 32 100 8-39
sleep 45
screen -dmS amrst_fine   bash "$D/run_durable.sh" laplace__fine   32 100 64-95
sleep 45
screen -dmS amrst_amr    bash "$D/run_durable.sh" laplace__amr    1  100 0
sleep 2
echo "launched:"; screen -ls | grep amrst
