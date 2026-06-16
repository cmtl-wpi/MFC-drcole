#!/usr/bin/env python3
"""Ma -> 0 convergence of the TC1 (Samareh Sec. 4.1.1) thermocapillary rise velocity with bulk
conduction.

Samareh's Test Case 1 is the Ma = 0 (invariant-temperature) limit: infinite thermal diffusivity
holds the imposed linear T field fixed, and the converged 2D cylinder ratio is v_t/v_YGB ~ 0.80
(their Fig 5). We realise that limit physically by turning bulk Fourier conduction ON at a sequence
of DECREASING thermal Marangoni numbers Ma (smaller Ma = larger conductivity = faster relaxation
back to the imposed gradient = closer to invariant T). The test is whether the quasi-steady plateau
v_t/v_YGB relaxes toward 0.80 as Ma -> 0.

For each run dir this calls measure.py, reads its RESULT_JSON plateau, and plots v_t/v_YGB vs Ma
with the Samareh 0.80 reference and a small-Ma linear extrapolation to Ma = 0. Runs that have not
produced snapshots yet are skipped, so it can be run incrementally as the sweep completes.

Usage: python3 plot_ma_convergence.py
"""
import json
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
SAMAREH = 0.80  # Samareh's converged 2D cylinder ratio (Fig 5)

# (Ma, run dir) -- the conduction Ma -> 0 convergence sweep (all w128, slip-wall, tr=2.0)
SWEEP = [
    (0.30, "tc1_cond_w128"),
    (0.10, "tc1_cond_w128_Ma010"),
    (0.05, "tc1_cond_w128_Ma005"),
    (0.03, "tc1_cond_w128_Ma003"),
]


def measure(run_dir):
    """Run measure.py on a run dir; return its RESULT_JSON dict, or None if not ready."""
    proc = subprocess.run(
        [sys.executable, os.path.join(HERE, "measure.py"), run_dir],
        capture_output=True, text=True, check=False,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


points = []
for ma, name in SWEEP:
    run_dir = os.path.join(RUNS, name)
    result = measure(run_dir) if os.path.isdir(run_dir) else None
    if result is None:
        print(f"Ma={ma:<5} {name}: no result yet")
        continue
    points.append((ma, result["ratio_plateau"], result["ratio_final"], result["slope_per_tr"], result["t_end_tr"], result["rises"]))
    print(f"Ma={ma:<5} plateau={result['ratio_plateau']:+.3f}  final={result['ratio_final']:+.3f}  drift={result['slope_per_tr']:+.3f}/t_r  ran {result['t_end_tr']:.2f} t_r  rises={result['rises']}")

if not points:
    sys.exit("no runs measured yet -- let the sweep produce snapshots first")

points.sort()  # ascending Ma
ma = np.array([p[0] for p in points])
plateau = np.array([p[1] for p in points])

fig, ax = plt.subplots(figsize=(7.2, 4.8))
ax.axhline(SAMAREH, ls="--", color="C3", lw=1.4, label=f"Samareh 2D = {SAMAREH:.2f}")
ax.plot(ma, plateau, "o-", color="C0", ms=7, lw=1.6, label=r"MFC conduction (plateau $v_t/v_{\mathrm{YGB}}$)")

# Linear extrapolation to Ma=0 from the two smallest-Ma points (the near-invariant end).
if len(ma) >= 2:
    order = np.argsort(ma)
    lo = order[:2]
    slope, intercept = np.polyfit(ma[lo], plateau[lo], 1)
    ax.plot([0, ma.max()], [intercept, slope * ma.max() + intercept], ":", color="0.5", lw=1.2, label=rf"linear extrap. $\to$ Ma=0: {intercept:.3f}")
    ax.plot(0, intercept, "*", color="k", ms=15, zorder=5)

ax.set_xlabel(r"thermal Marangoni number  $Ma$   (smaller $\to$ closer to invariant $T$)")
ax.set_ylabel(r"quasi-steady  $v_t / v_{\mathrm{YGB}}$")
ax.set_xlim(left=-0.012)
ax.set_title(r"TC1 (Samareh 4.1.1): conduction $Ma\to0$ convergence to the invariant-$T$ limit")
ax.grid(alpha=0.3)
ax.legend(loc="best", fontsize=9)
fig.tight_layout()
out = os.path.join(HERE, "figures", "tc1_ma_convergence.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150)
print("saved ->", out)
