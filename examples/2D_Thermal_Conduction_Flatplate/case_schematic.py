#!/usr/bin/env python3
"""Schematic of the 2D_Thermal_Conduction_Flatplate case setup: domain, boundary
conditions, initial state, and key parameters. Values mirror case.py.

    python3 examples/2D_Thermal_Conduction_Flatplate/case_schematic.py
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.special import erf

HERE = os.path.dirname(os.path.abspath(__file__))

# -- setup values (mirror case.py) --
L = 50.0  # domain size [mm]  (0.05 m)
ncell = 199  # cells per direction
T_inf, T_wall = 1125.0, 600.0
p_atm = 101325.0
gam, R_air = 1.4, 287.0
k_air = 0.07
Re = 100000
rho_inf = p_atm / (R_air * T_inf)
cp = gam * R_air / (gam - 1.0)
rho_film = rho_inf * T_inf / (0.5 * (T_inf + T_wall))
alpha_film = k_air / (rho_film * cp)

cmap = plt.cm.inferno
c_free = cmap(0.999)  # free-stream (hot) color
band = 9.0  # schematic thermal-layer thickness [mm] (exaggerated for clarity)

fig = plt.figure(figsize=(13.5, 6.6))
gs = fig.add_gridspec(1, 2, width_ratios=[2.1, 1.12])
ax = fig.add_subplot(gs[0])
ax.set_aspect("equal")
ax.set_xlim(-11, 60)
ax.set_ylim(-11, 61)
ax.axis("off")

# -- domain: uniform hot free stream + exaggerated near-wall thermal layer --
ax.add_patch(Rectangle((0, band), L, L - band, facecolor=c_free, edgecolor="none"))
grad = erf(np.linspace(0, 2.0, 64))[:, None] * np.ones((1, 2))  # cold wall -> hot
ax.imshow(grad, extent=(0, L, 0, band), origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=1, zorder=1)
ax.add_patch(Rectangle((0, 0), L, L, fill=False, edgecolor="k", lw=1.8, zorder=4))
ax.plot([0, L], [band, band], "w--", lw=1.0, zorder=3)

# -- bottom: the flat plate (no-slip isothermal wall) --
ax.add_patch(Rectangle((0, -2.6), L, 2.6, facecolor="0.8", edgecolor="k", hatch="////", lw=1.5, zorder=4))
ax.plot([0, L], [0, 0], "k-", lw=3.0, zorder=5)
ax.text(L / 2, -6.6, "Flat plate — no-slip isothermal wall,  $T_{wall}=600$ K   (bc_y%beg = −16)", ha="center", va="center", fontsize=10.5, fontweight="bold")

# -- edge boundary-condition arrows --
# left edge: subsonic inflow (arrows pointing in, +x)
for f in np.linspace(0.18, 0.82, 4):
    ax.add_patch(FancyArrowPatch((-6.5, L * f), (-1.0, L * f), arrowstyle="-|>", mutation_scale=13, color="C0", lw=1.6, zorder=6))
ax.text(-9.5, L / 2, "Subsonic inflow\n(bc_x%beg = −7)", rotation=90, ha="center", va="center", fontsize=10, color="C0")

# right edge: outflow (arrows pointing out, +x)
for f in np.linspace(0.18, 0.82, 4):
    ax.add_patch(FancyArrowPatch((L + 1.0, L * f), (L + 6.5, L * f), arrowstyle="-|>", mutation_scale=13, color="0.35", lw=1.6, zorder=6))
ax.text(L + 9.5, L / 2, "Outflow\n(bc_x%end = −3)", rotation=90, ha="center", va="center", fontsize=10, color="0.3")

# top edge: outflow (arrows pointing out, +y)
for f in np.linspace(0.15, 0.85, 5):
    ax.add_patch(FancyArrowPatch((L * f, L + 1.0), (L * f, L + 6.5), arrowstyle="-|>", mutation_scale=13, color="0.35", lw=1.6, zorder=6))
ax.text(L / 2, L + 9.0, "Outflow  (bc_y%end = −3)", ha="center", va="center", fontsize=10, color="0.3")

# -- initial-condition annotation --
ax.text(
    L / 2,
    L * 0.62,
    "Initial condition ($t=0$)\nQuiescent:  $u=v=0$\n$T_\\infty=1125$ K,  $p=1$ atm\n(uniform)",
    ha="center",
    va="center",
    fontsize=11,
    bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85, edgecolor="0.5"),
)

# -- thermal-layer callout --
ax.annotate(
    "Thermal boundary layer\n$\\delta_T=\\sqrt{\\alpha t}\\approx 3$ mm at 5 ms\n(uniform in $x$, grows in time)",
    xy=(L * 0.80, band * 0.5),
    xytext=(L * 0.62, band + 9.5),
    fontsize=9.5,
    ha="center",
    arrowprops=dict(arrowstyle="->", color="k", lw=1.2),
    bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff4e6", edgecolor="0.5"),
)
ax.text(2.0, band + 1.2, "(layer thickness exaggerated)", fontsize=8, style="italic", color="0.4")

# -- coordinate axes (at the physical origin; light so they read over the dark band) --
ax.add_patch(FancyArrowPatch((1.2, 1.2), (13, 1.2), arrowstyle="-|>", mutation_scale=14, color="white", lw=2.0, zorder=8))
ax.add_patch(FancyArrowPatch((1.2, 1.2), (1.2, 13), arrowstyle="-|>", mutation_scale=14, color="white", lw=2.0, zorder=8))
ax.text(14.6, 1.2, "$x$", fontsize=12, va="center", color="white", fontweight="bold", zorder=8)
ax.text(1.2, 15.0, "$y$", fontsize=12, ha="center", color="k", fontweight="bold")
# 10 mm scale bar in the clear free-stream region
ax.plot([30, 40], [46, 46], "k-", lw=2.5, zorder=6)
ax.text(35, 47.8, "10 mm", ha="center", fontsize=9)

ax.set_title("2D Thermal-Conduction Flat Plate — case setup", fontsize=13, fontweight="bold", pad=10)

# -- parameter panel --
axp = fig.add_subplot(gs[1])
axp.axis("off")
lines = [
    ("Geometry", ""),
    ("  domain", f"{L:.0f} × {L:.0f} mm"),
    ("  grid", f"{ncell}×{ncell} cells  (Δ = {L / (ncell + 1):.3f} mm)"),
    ("Fluid — air (ideal gas)", ""),
    ("  γ,  R", "1.4,  287 J/(kg·K)"),
    ("  cp,  cv", f"{cp:.0f},  {R_air / (gam - 1):.0f} J/(kg·K)"),
    ("  π∞ (stiffened gas)", "0"),
    ("Free stream", ""),
    ("  T∞,  p", "1125 K,  1 atm"),
    ("  ρ∞", f"{rho_inf:.3f} kg/m³"),
    ("Wall", ""),
    ("  T_wall", "600 K   (ΔT = 525 K)"),
    ("Conduction (Fourier)", ""),
    ("  k", f"{k_air} W/(m·K)"),
    ("  α_film", f"{alpha_film:.2e} m²/s"),
    ("Viscous", ""),
    ("  Re", f"{Re:,}"),
    ("Numerics", ""),
    ("  model_eqns / scheme", "2 / WENO5 + HLLC"),
    ("  dt,  t_end", "1.12e−7 s,  5 ms"),
    ("  steps,  saves", "≈44 600,  10"),
    ("Features", ""),
    ("  thermal_conduction", "T"),
    ("  chemistry", "F"),
]
y = 0.99
for label, val in lines:
    header = val == ""
    axp.text(0.0, y, label, fontsize=10.0 if header else 9.3, fontweight="bold" if header else "normal", color="#7a3b00" if header else "k", transform=axp.transAxes)
    if val:
        axp.text(0.99, y, val, fontsize=9.0, ha="right", transform=axp.transAxes, family="monospace")
    y -= 0.0405
axp.set_title("Case parameters", fontsize=12, fontweight="bold", loc="left")

fig.tight_layout()
out = os.path.join(HERE, "case_setup.png")
fig.savefig(out, dpi=140, bbox_inches="tight")
print(f"wrote {out}")
