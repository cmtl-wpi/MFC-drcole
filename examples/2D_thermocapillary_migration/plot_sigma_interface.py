#!/usr/bin/env python3
"""Surface tension along the droplet interface -- visualizes the Marangoni driver.

MFC never stores sigma, so we rebuild it from the coded closure
    sigma(T) = sigma0 + (dsigma/dT) * (T - T_ref),
with the temperature T recovered per cell from the stiffened-gas EOS (same recovery as
plot_temperature.py). The interface is the thin band of cells where the color function c is
partway between background (0) and droplet (1); we read sigma in those cells and plot it against
angle around the droplet. sigma comes out lowest on the hot (+x) side and highest on the cold
(-x) side -- that tangential gradient is what pushes the droplet up the temperature gradient.

Usage:  python3 plot_sigma_interface.py [case_dir] [step]   (step defaults to the last snapshot)
"""
import glob
import os
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

case_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

# Read every constant from the namelist MFC actually ran (simulation.inp), so nothing is hardcoded
# and the script works on any run/grid. The file is plain "name = value" lines.
params = {}
for line in open(os.path.join(case_dir, "simulation.inp")):
    if "=" in line:
        name, value = line.split("=", 1)
        params[name.strip().lower()] = value.strip().rstrip(",")


def param(name):
    """Look up a parameter by name in the parsed namelist and return it as a float."""
    return float(params[name.lower()])


nx, ny = int(param("m")) + 1, int(param("n")) + 1  # number of interior cells in x, y
dt = param("dt")
cv = param("fluid_pp(1)%cv")

# MFC stores gamma_mfc = 1/(gamma-1) and pi_inf_mfc = gamma*p_inf/(gamma-1); invert for the
# physical gamma and p_inf, which the temperature formula below needs.
gamma_mfc = param("fluid_pp(1)%gamma")
pi_inf_mfc = param("fluid_pp(1)%pi_inf")
gamma = 1.0 + 1.0 / gamma_mfc
p_inf = pi_inf_mfc * (gamma - 1.0) / gamma

# The sigma(T) closure coefficients.
sigma0, dsigma_dT, T_ref = param("sigma"), param("sigma_dtdt"), param("sigma_t_ref")

# Cell-center coordinates: average adjacent cell boundaries (the last nx+1 / ny+1 are the interior).
restart_dir = os.path.join(case_dir, "restart_data")
xb = np.fromfile(os.path.join(restart_dir, "lustre_x_cb.dat"), np.float64)[-(nx + 1):]
yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1):]
x = 0.5 * (xb[:-1] + xb[1:])
y = 0.5 * (yb[:-1] + yb[1:])

# Find the available snapshots. Files are lustre_<step>.dat; skip lustre_x_cb.dat / lustre_y_cb.dat
# (the coordinate files, which carry no step number). Default to the last snapshot.
steps = []
for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")):
    match = re.search(r"lustre_(\d+)\.dat$", f)
    if match:
        steps.append(int(match.group(1)))
steps.sort()
step = int(sys.argv[2]) if len(sys.argv) > 2 else steps[-1]

# A snapshot is a flat array holding (ny*nx) values per conserved variable. For this case the
# variable order is: 0,1 = partial densities, 7,8 = phasic internal energies, last = color c.
columns = np.fromfile(os.path.join(restart_dir, f"lustre_{step}.dat"), np.float64).reshape(-1, ny * nx)


def conserved_var(index):
    """Conserved variable `index` of this snapshot, reshaped from a flat array to the (y, x) grid."""
    return columns[index].reshape(ny, nx)


# Variable layout (see note above): rows 0,1 = partial densities, 7,8 = phasic internal energies.
rho = conserved_var(0) + conserved_var(1)
rho_e = conserved_var(7) + conserved_var(8)               # internal energy density (no kinetic part)
pressure = (rho_e - pi_inf_mfc) / gamma_mfc                # invert rho*e = gamma_mfc*p + pi_inf_mfc
T = (pressure + p_inf) / ((gamma - 1.0) * rho * cv)        # stiffened-gas temperature
sigma = sigma0 + dsigma_dT * (T - T_ref)                   # the coded closure, evaluated per cell
c = np.clip(conserved_var(-1), 0.0, None)                  # color function (last variable)

# Droplet center = color-weighted average position (the color centroid).
xc = (c * x[None, :]).sum() / c.sum()
yc = (c * y[:, None]).sum() / c.sum()

# Interface cells: c between background and droplet. Read sigma there and find each cell's angle
# about the center (0 deg points to +x = hot, +-180 deg points to -x = cold).
X, Y = np.meshgrid(x, y)
interface = (c > 0.2) & (c < 0.8)
angle = np.degrees(np.arctan2(Y[interface] - yc, X[interface] - xc))
sigma_on_interface = sigma[interface]

fig, (ax_field, ax_line) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.05, 1]})

# Left panel: the sigma field with the interface outlined.
mesh = ax_field.pcolormesh(x, y, sigma, cmap="viridis", shading="auto")
fig.colorbar(mesh, ax=ax_field, label=r"$\sigma(T)$")
ax_field.contour(x, y, c, levels=[0.5], colors="w", linewidths=1.2)
ax_field.set_aspect("equal")
ax_field.set_xlabel("x")
ax_field.set_ylabel("y")
ax_field.set_title(rf"$\sigma(T)$ field + interface (step {step}, t = {step * dt:.2f})")

# Right panel: sigma in the interface cells vs angle -- the Marangoni gradient, made literal.
ax_line.scatter(angle, sigma_on_interface, s=12, color="C0")
ax_line.axvline(0, ls=":", color="C3")
ax_line.text(6, sigma_on_interface.min(), "hot (+x)", color="C3", fontsize=9)
ax_line.set_xlabel("angle around interface (deg)")
ax_line.set_ylabel(r"$\sigma$ in the interface cells")
ax_line.set_title(r"Low $\sigma$ on the hot side drives the Marangoni pull")
ax_line.grid(alpha=0.3)

fig.tight_layout()
out = os.path.join(case_dir, "viz", f"sigma_interface_{step}.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150)
print(f"sigma on interface spans [{sigma_on_interface.min():.4f}, {sigma_on_interface.max():.4f}]")
print(f"saved -> {out}")
