#!/usr/bin/env python3
"""Plot the temperature field of the thermocapillary run (left) and its centerline profile (right).

Temperature is NOT a stored output in MFC -- it is recovered per cell from the stiffened-gas EOS,
    T = (p + p_inf) / ((gamma - 1) * rho * cv),
with the pressure obtained by inverting the conserved internal energy. (At step 0 this reproduces the
imposed IC T(x) = T0 + gradT*x to ~1e-8.) MFC has no bulk heat conduction, so that linear field is a
frozen initial condition the flow slowly advects -- the right panel makes that distortion visible.

Usage:  python3 plot_temperature.py [case_dir] [step]   (step defaults to the last snapshot; e.g. 3000)
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

# Read constants from the namelist MFC ran (simulation.inp), never hardcode them: a stale Nx silently
# mis-slices the data into nonsense (negative "temperatures") with no error. The file regenerates every
# run, so it can't disagree with the data beside it.
params = {}
for line in open(os.path.join(case_dir, "simulation.inp")):
    if "=" in line:
        name, value = line.split("=", 1)
        params[name.strip().lower()] = value.strip().rstrip(",")


def param(name):
    """Look up a parameter by name in the parsed namelist and return it as a float."""
    return float(params[name.lower()])


nx, ny = int(param("m")) + 1, int(param("n")) + 1
dt = param("dt")
cv = param("fluid_pp(1)%cv")

# MFC stores gamma_mfc = 1/(gamma-1) and pi_inf_mfc = gamma*p_inf/(gamma-1); invert for physical values.
gamma_mfc = param("fluid_pp(1)%gamma")
pi_inf_mfc = param("fluid_pp(1)%pi_inf")
gamma = 1.0 + 1.0 / gamma_mfc
p_inf = pi_inf_mfc * (gamma - 1.0) / gamma

# The imposed linear T(x) = T0 + gradT*x lives in case.py's analytic density string, not the namelist,
# so it is set here only to draw the "frozen initial field" reference line.
T0, gradT = 10.0, 1.0

restart_dir = os.path.join(case_dir, "restart_data")
ncell = nx * ny

# Cell-center coordinates from the boundary files (last nx+1 / ny+1 boundaries are the interior cells).
xb = np.fromfile(os.path.join(restart_dir, "lustre_x_cb.dat"), np.float64)[-(nx + 1):]
yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1):]
x = 0.5 * (xb[:-1] + xb[1:])
y = 0.5 * (yb[:-1] + yb[1:])

# Available snapshots: lustre_<step>.dat (skip the lustre_x_cb / lustre_y_cb coordinate files).
steps = []
for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")):
    match = re.search(r"lustre_(\d+)\.dat$", f)
    if match:
        steps.append(int(match.group(1)))
steps.sort()
if not steps:
    sys.exit(f"no lustre_<step>.dat snapshots in {restart_dir!r} -- run the case first")
step = int(sys.argv[2]) if len(sys.argv) > 2 else steps[-1]
if step not in steps:
    sys.exit(f"step {step} not available; choose from {steps[0]}..{steps[-1]}")


def temperature(snapshot_step):
    """Reconstruct (T, color) on the (y, x) grid for a snapshot, via the stiffened-gas EOS."""
    flat = np.fromfile(os.path.join(restart_dir, f"lustre_{snapshot_step}.dat"), np.float64)
    if flat.size % ncell != 0:
        sys.exit(f"snapshot {snapshot_step}: {flat.size} values not divisible by ncell={ncell} -- grid mismatch.")
    # One conserved variable per row, ncell cells per row.
    cons_var = flat.reshape(-1, ncell)
    # Conserved-variable layout for this case: rows 0,1 = partial densities (one per fluid),
    # rows 7,8 = phasic internal energies, last row = color function.
    rho = (cons_var[0] + cons_var[1]).reshape(ny, nx)
    rho_e = (cons_var[7] + cons_var[8]).reshape(ny, nx)
    # Invert MFC's stored internal energy to pressure, then apply the stiffened-gas EOS for T.
    p = (rho_e - pi_inf_mfc) / gamma_mfc
    T = (p + p_inf) / ((gamma - 1.0) * rho * cv)
    color = np.clip(cons_var[-1].reshape(ny, nx), 0.0, None)
    return T, color


T, c = temperature(step)
T_initial, _ = temperature(steps[0])          # step-0 field, for the centerline comparison
t_phys = step * dt

# Guard: if the constants disagreed with the data, T comes out non-physical -- fail loudly, don't plot garbage.
if not np.all(np.isfinite(T)) or T.min() <= 0.0:
    sys.exit(f"reconstructed T is non-physical (min = {T.min():.4f}); constants likely mismatch the data.")

print(f"step {step}  (t = {t_phys:.4f})   T range: [{T.min():.4f}, {T.max():.4f}]")

fig, (ax_field, ax_line) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0]})

# Left: the temperature field, with the droplet interface (c=0.5) and its centroid marked.
vmin = min(T0 + gradT * x[0], T.min())
vmax = max(T0 + gradT * x[-1], T.max())
mesh = ax_field.pcolormesh(x, y, T, cmap="coolwarm", vmin=vmin, vmax=vmax, shading="auto")
fig.colorbar(mesh, ax=ax_field, label="temperature $T$")
ax_field.contour(x, y, c, levels=[0.5], colors="k", linewidths=1.2)
# Droplet centroid = color-weighted average position.
x_centroid = (c * x[None, :]).sum() / c.sum()
y_centroid = (c * y[:, None]).sum() / c.sum()
ax_field.plot(x_centroid, y_centroid, "k+", ms=10, mew=2)
ax_field.set_aspect("equal")
ax_field.set_xlabel("$x$")
ax_field.set_ylabel("$y$")
ax_field.set_title(f"Temperature field (step {step}, $t={t_phys:.3f}$)\ndroplet centroid $x={x_centroid:+.4f}$")

# Right: centerline T(x) now vs the frozen initial linear field -- the gap is the advective distortion.
mid = ny // 2
ax_line.plot(x, T0 + gradT * x, "--", color="0.5", label=r"initial $T(x)=T_0+\nabla T\,x$ (frozen IC)")
ax_line.plot(x, T_initial[mid], ":", color="C0", alpha=0.7, label=f"centerline, step {steps[0]}")
ax_line.plot(x, T[mid], "-", color="C3", label=f"centerline, step {step}")
ax_line.set_xlabel("$x$")
ax_line.set_ylabel(r"temperature $T$ (at $y \approx 0$)")
ax_line.set_title("Centerline temperature profile")
ax_line.legend(loc="upper left", fontsize=8)
ax_line.grid(alpha=0.3)

fig.tight_layout()
viz_dir = os.path.join(case_dir, "viz")
os.makedirs(viz_dir, exist_ok=True)
out_png = os.path.join(viz_dir, f"temperature_{step}.png")
fig.savefig(out_png, dpi=150)
print(f"saved figure -> {out_png}")
