#!/usr/bin/env python3
"""Visualize the REAL thermocapillary recirculation from an MFC 2D run, as the data
analogue of figures/thermocapillary_migration.tex (the schematic).

Two panels, both in the drop's co-moving frame (drop velocity subtracted, so the closed
internal circulation is visible rather than a translating blob):
  (a) streamlines colored by flow speed + the c = 0.5 interface  -> the recirculation loops
  (b) raw cell-resolved vorticity omega_z = dv/dx - du/dy        -> the counter-rotating cells

The current case migrates along +y (the 7.5D gradient axis): the cold wall is the floor
(y%beg), the hot wall the ceiling (y%end), and the drop starts 1.5D above the cold floor
(y_drop = -2.25 in the slip-wall box). The viscous time tau = rho*r^2/mu is read from the
data (color-weighted drop density at t=0), so it tracks the case (rho = 0.2 -> tau = 0.5),
not a hard-coded value.

Run-dependent constants are read from <run_dir>/simulation.inp so this can't silently
disagree with the data. Conserved-variable layout (model_eqns=3, num_fluids=2): fields 0,1 =
partial densities, 2 = x-momentum, 3 = y-momentum, then color c (last, or second-to-last when
a thermal_scalar T_s is appended).

Usage:  python3 plot_recirculation.py runs/fig5_2D_w256 [t_over_tau_target=2.6]
Writes: figures/recirculation_2D_w256.png / .pdf
"""

import glob
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

run_dir = sys.argv[1] if len(sys.argv) > 1 else "runs/fig5_2D_w256"
target_ttau = float(sys.argv[2]) if len(sys.argv) > 2 else 2.6
here = os.path.dirname(os.path.abspath(__file__))
run_dir = run_dir if os.path.isabs(run_dir) else os.path.join(here, run_dir)


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


P = read_namelist(os.path.join(run_dir, "simulation.inp"))
nx = int(P["m"]) + 1
ny = int(P["n"]) + 1
dt = float(P["dt"])
mu = 1.0 / float(P["fluid_pp(1)%re(1)"])
dsigma_dT = float(P["sigma_dtdt"])
ts_mode = str(P.get("thermal_scalar", "F")).strip(". ").upper().startswith("T")
r = 0.5
gradT = 2.0 / 15.0
v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu

restart = os.path.join(run_dir, "restart_data")
xb = np.fromfile(os.path.join(restart, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
yb = np.fromfile(os.path.join(restart, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
x = 0.5 * (xb[:-1] + xb[1:])
y = 0.5 * (yb[:-1] + yb[1:])

steps = sorted(
    int(re.search(r"lustre_(\d+)\.dat$", f).group(1))
    for f in glob.glob(os.path.join(restart, "lustre_*.dat"))
    if re.search(r"lustre_(\d+)\.dat$", f)
)
cells = nx * ny
nvars = np.fromfile(os.path.join(restart, f"lustre_{steps[0]}.dat"), np.float64).size // cells
c_idx = nvars - 2 if ts_mode else nvars - 1  # color c (T_s is appended after it in ts mode)


def snapshot(step):
    snap = np.fromfile(os.path.join(restart, f"lustre_{step}.dat"), np.float64)
    fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(ny, nx)  # noqa: E731
    rho = fld(0) + fld(1)
    return rho, fld(2) / rho, fld(3) / rho, np.clip(fld(c_idx), 0.0, 1.0)  # rho, u, v, color


# Viscous time tau = rho*r^2/mu from the t=0 drop density (rho = 0.2 here -> tau = 0.5), so the
# time label tracks the case rather than assuming rho = 1.
rho0, _, _, c0 = snapshot(steps[0])
rho_drop = (c0 * rho0).sum() / c0.sum()
tau = rho_drop * r**2 / mu

# Pick the snapshot closest to the requested t/tau (well inside the quasi-steady window).
step = min(steps, key=lambda s: abs(s * dt / tau - target_ttau))
ttau = step * dt / tau
rho, u, v, c = snapshot(step)

# Co-moving frame: subtract the color-weighted mean drop velocity (both components). The drop
# migrates in +y, so v_drop is the migration speed; u_drop ~ 0 (symmetric about x = 0).
u_drop = (c * u).sum() / c.sum()
v_drop = (c * v).sum() / c.sum()
uc = u - u_drop
vc = v - v_drop

# Vorticity on the raw grid (real spacing, no interpolation).
omega = np.gradient(vc, x, axis=1) - np.gradient(uc, y, axis=0)

# Zoom to a window around the drop. It migrates along +y and starts off-center (y ~ -2.25 in the
# slip-wall box), so center the window on the measured drop centroid in BOTH axes.
xc_drop = (c * x[None, :]).sum() / c.sum()
yc_drop = (c * y[:, None]).sum() / c.sum()
win = 2.5 * r  # show the drop plus a ring of external streaming flow
mx = (x > xc_drop - win) & (x < xc_drop + win)
my = (y > yc_drop - win) & (y < yc_drop + win)
xz, yz = x[mx], y[my]
uz, vz = uc[np.ix_(my, mx)], vc[np.ix_(my, mx)]
cz, oz = c[np.ix_(my, mx)], omega[np.ix_(my, mx)]
speed = np.hypot(uz, vz)

# journal styling: serif text + Computer-Modern math, to match the schematic
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "axes.linewidth": 0.8,
})
coldc = (46 / 255, 86 / 255, 149 / 255)  # palette shared with the schematic
hotc = (171 / 255, 57 / 255, 52 / 255)
inkc = (0.12, 0.12, 0.12)
halo = [pe.withStroke(linewidth=2.4, foreground="white")]  # keep labels legible on any field

fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.2, 3.9), sharey=True, constrained_layout=True)

# (a) co-moving streamlines colored by speed + interface
strm = axa.streamplot(
    xz, yz, uz, vz, color=speed / v_YGB, cmap="viridis",
    density=1.3, linewidth=0.75, arrowsize=0.6,
)
cb = fig.colorbar(strm.lines, ax=axa, fraction=0.046, pad=0.03)
cb.set_label(r"$|\mathbf{u}-\mathbf{U}_{\rm drop}|\,/\,v_{\rm YGB}$")
cb.outline.set_linewidth(0.6)
axa.contour(xz, yz, cz, levels=[0.5], colors="white", linewidths=2.2)
axa.contour(xz, yz, cz, levels=[0.5], colors=[inkc], linewidths=0.9)
axa.set_title("(a) drop-frame streamlines")

# (b) cell-resolved vorticity (nondimensionalized by r / v_YGB) + interface
o_nd = oz * r / v_YGB
vmax = np.percentile(np.abs(o_nd), 99)
pm = axb.pcolormesh(xb[mx.nonzero()[0][0] : mx.nonzero()[0][-1] + 2],
                    yb[my.nonzero()[0][0] : my.nonzero()[0][-1] + 2],
                    o_nd, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat", rasterized=True)
cb2 = fig.colorbar(pm, ax=axb, fraction=0.046, pad=0.03)
cb2.set_label(r"$\omega_z\, r / v_{\rm YGB}$")
cb2.outline.set_linewidth(0.6)
axb.contour(xz, yz, cz, levels=[0.5], colors="white", linewidths=2.0)
axb.contour(xz, yz, cz, levels=[0.5], colors=[inkc], linewidths=0.9)
axb.set_title(r"(b) vorticity field")

for ax in (axa, axb):
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/D$")
    # migration direction = toward hot (+y); white underlay keeps it readable on the field
    for col, lw in (("white", 4.0), (inkc, 1.8)):
        ax.annotate("", xy=(xc_drop, yc_drop + 1.0 * r), xytext=(xc_drop, yc_drop - 0.05 * r),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, shrinkA=0, shrinkB=0))
    # cold (floor, y%beg) / hot (ceiling, y%end) orientation, tying back to the schematic
    ax.text(0.04, 0.04, "cold", transform=ax.transAxes, color=coldc, fontsize=8,
            ha="left", va="bottom", path_effects=halo)
    ax.text(0.04, 0.96, "hot", transform=ax.transAxes, color=hotc, fontsize=8,
            ha="left", va="top", path_effects=halo)
axa.set_ylabel(r"$y/D$")
axa.text(xc_drop + 0.28 * r, yc_drop + 0.15 * r, rf"$U={v_drop / v_YGB:.2f}\,v_{{\rm YGB}}$",
         color=inkc, fontsize=9, ha="left", va="bottom", path_effects=halo)

fig.suptitle(rf"2D thermocapillary drop  $\cdot$  {ny / 7.5:.1f} cells/$D$  $\cdot$  $t={ttau:.1f}\,\tau$",
             fontsize=9.5)

out = os.path.join(here, "figures", "recirculation_2D_w256")
fig.savefig(out + ".png", dpi=300)
fig.savefig(out + ".pdf")
print(f"step={step}  t/tau={ttau:.3f}  nvars={nvars}  rho_drop={rho_drop:.3f}  tau={tau:.3f}")
print(f"U_drop/v_YGB={v_drop / v_YGB:+.3f}  (lateral u_drop/v_YGB={u_drop / v_YGB:+.3f})")
print(f"saved -> {out}.png / .pdf")
