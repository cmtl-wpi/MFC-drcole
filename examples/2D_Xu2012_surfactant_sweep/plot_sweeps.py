#!/usr/bin/env python3
# Plot the M2 property-sweep responses from results.jsonl: D vs Ca, D vs viscosity ratio, D vs Re, and
# (D + surfactant non-uniformity) vs Pe. The shared baseline (Ca=0.3, Pe=10, lambda=1, Re=1) anchors
# every panel. Run from repo root after run_sweeps.sh.
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EX = "examples/2D_Xu2012_surfactant_sweep"
rows = [json.loads(l) for l in open(os.path.join(EX, "results.jsonl"))]
base = next(r for r in rows if r["group"] == "ca" and abs(r["x"] - 0.3) < 1e-9)


def series(group, extra_x, extra_pt):
    pts = [(r["x"], r) for r in rows if r["group"] == group]
    pts.append((extra_x, extra_pt))  # inject the shared baseline point
    pts = sorted({round(x, 4): (x, r) for x, r in pts}.values())
    return [x for x, _ in pts], [r["D"] for _, r in pts], [r["surf_nonunif"] for _, r in pts]


fig, ax = plt.subplots(2, 2, figsize=(11, 8))

cx, cy, _ = series("ca", 0.3, base)
ax[0, 0].plot(cx, cy, "o-", color="C0")
ax[0, 0].set(xlabel="Ca", ylabel="Taylor deformation D", title="Ca ↑  →  more elongation")

lx, ly, _ = series("lam", 1.0, base)
ax[0, 1].plot(lx, ly, "s-", color="C3")
ax[0, 1].set(xlabel="viscosity ratio λ = μ_drop/μ_matrix", ylabel="D", title="λ ↑  →  less deformation")

rx, ry, _ = series("re", 1.0, base)
ax[1, 0].plot(rx, ry, "^-", color="C2")
ax[1, 0].set(xlabel="Re", ylabel="D", title="Re ↑  →  more elongation")

px, py, pn = series("pe", 10.0, base)
ax[1, 1].plot(px, pn, "d-", color="C1")
ax[1, 1].set_xscale("log")
ax[1, 1].set(xlabel="Pe = γ̇R²/D_s", ylabel="Γ non-uniformity (P90/median)", title="Pe ↑  →  surfactant more tip-concentrated")

for a in ax.flat:
    a.grid(alpha=0.3)
fig.suptitle("Xu 2012 (M2): finite-Re surfactant drop in shear — property-ratio responses", fontsize=12)
fig.tight_layout()
out = os.path.join(EX, "figures", "sweeps.png")
fig.savefig(out, dpi=130)
print("wrote", out)
