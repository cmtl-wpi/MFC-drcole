#!/usr/bin/env python3
# Plot the sphere surface-diffusion convergence: measured l=1 decay rate / exact 2 D_s/R^2 vs R/dx.
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.loadtxt(os.path.join(HERE, "convergence.dat"))
rdx, meas, exact = d[:, 0], d[:, 1], d[:, 2]
ratio = meas / exact

fig, ax = plt.subplots(figsize=(6.2, 4.4))
ax.axhline(1.0, color="#333", lw=2, label=r"exact  $2 D_s/R^2$")
ax.plot(rdx, ratio, "o-", ms=8, mfc="#d1495b", mec="k", mew=0.6, label="MFC (l=1)")
for x, y in zip(rdx, ratio):
    ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(6, 6))
ax.set_xlabel(r"interface resolution  $R/\Delta x$")
ax.set_ylabel(r"measured rate / exact")
ax.set_title(r"Sphere surface diffusion: convergence to $l(l+1)D_s/R^2$ ($l=1$)")
ax.set_ylim(0, 1.15)
ax.grid(True, alpha=0.25)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "convergence.png"), dpi=130)
print("wrote figures/convergence.png")
