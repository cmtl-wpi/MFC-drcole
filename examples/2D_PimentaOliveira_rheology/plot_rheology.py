#!/usr/bin/env python3
# Plot the M3 rheology decomposition from results.jsonl: capillary [eta_c], Marangoni [eta_m], and total
# [eta]=[eta_c]+[eta_m] vs surfactant coverage X, plus the first normal-stress difference N1. The gate is
# visible here: [eta_m] is ~0 for the clean drop (X=0) and grows with coverage, while N1 stays positive.
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EX = "examples/2D_PimentaOliveira_rheology"
rows = sorted((json.loads(l) for l in open(os.path.join(EX, "results.jsonl"))), key=lambda r: r["X"])
X = [r["X"] for r in rows]
ec = [r["eta_capillary"] for r in rows]
em = [r["eta_marangoni"] for r in rows]
et = [r["eta_intrinsic"] for r in rows]
N1 = [r["N1_star"] for r in rows]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
w = 0.05
xi = np.arange(len(X))
ax[0].bar(xi - w * 6, ec, width=w * 5, label="capillary [η_c]", color="C0")
ax[0].bar(xi - w, em, width=w * 5, label="Marangoni [η_m]", color="C1")
ax[0].bar(xi + w * 4, et, width=w * 5, label="total [η]=[η_c]+[η_m]", color="C2")
ax[0].set_xticks(xi)
ax[0].set_xticklabels([f"X={x}" for x in X])
ax[0].set(ylabel="intrinsic viscosity", title="Interfacial-stress decomposition")
ax[0].legend(fontsize=9)
ax[0].grid(axis="y", alpha=0.3)

ax[1].plot(X, em, "o-", color="C1", label="Marangoni [η_m]")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set(xlabel="coverage X = Γ₀/Γ∞", ylabel="[η_m]", title="Marangoni stress grows with coverage")
ax[1].grid(alpha=0.3)
axt = ax[1].twinx()
axt.plot(X, N1, "s--", color="C3", label="N₁ (normal-stress diff)")
axt.set_ylabel("N₁*", color="C3")
axt.tick_params(axis="y", labelcolor="C3")

fig.suptitle("Pimenta & Oliveira (M3): surfactant-drop rheology in shear (Ca=0.1)", fontsize=12)
fig.tight_layout()
out = os.path.join(EX, "figures", "rheology.png")
fig.savefig(out, dpi=130)
print("wrote", out)
