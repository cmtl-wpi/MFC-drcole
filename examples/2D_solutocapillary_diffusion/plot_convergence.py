#!/usr/bin/env python3
# Plot the 2D circle surface-diffusion "convergence": the whole-field moment (full-field) and the
# band-only moment BRACKET the exact rate rather than converging to it -- exposing that the moment
# estimator is biased on a staircased curved interface (see README).
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.loadtxt(os.path.join(HERE, "convergence.dat"))
rdx, full, band, exact = d[:, 0], d[:, 1], d[:, 2], d[:, 3]

fig, ax = plt.subplots(figsize=(6.4, 4.6))
ax.axhline(1.0, color="#333", lw=2, label=r"exact  $D_s/R^2$")
ax.fill_between(rdx, full / exact, band / exact, color="#d1495b", alpha=0.12)
ax.plot(rdx, full / exact, "o-", ms=7, color="#2e6f95", mec="k", mew=0.5, label="full-field moment (biased low)")
ax.plot(rdx, band / exact, "s-", ms=7, color="#d1495b", mec="k", mew=0.5, label="band-only moment (biased high)")
ax.set_xlabel(r"interface resolution  $R/\Delta x$")
ax.set_ylabel(r"measured rate / exact")
ax.set_title("2D circle: the moment estimator brackets the exact rate")
ax.grid(True, alpha=0.25)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figures", "convergence.png"), dpi=130)
print("wrote figures/convergence.png")
