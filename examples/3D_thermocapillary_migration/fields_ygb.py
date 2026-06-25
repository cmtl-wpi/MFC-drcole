#!/usr/bin/env python3
"""Sanity-plot one 3D thermocapillary run before trusting a multi-hour sweep.

Renders an x-y midplane (z = Lz/2) of the EOS-derived temperature T with the drop's color
interface overlaid, plus the centerline T(y) profile against the imposed linear field and the
isothermal-wall temperatures. Confirms the IC is what case_ygb.py intended: a linear, wall-pinned
T and a coherent spherical color patch.

NOTE: named fields_ygb.py (not plot*.py) on purpose -- the repo root .gitignore hides examples/**/p*.

Usage:  python3 fields_ygb.py <run_dir> [step]   (default step: the last snapshot)
"""

import glob
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

run_dir = sys.argv[1] if len(sys.argv) > 1 else "."


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().lower()] = v.strip().rstrip(",")
    return out


P = read_namelist(os.path.join(run_dir, "simulation.inp"))
nx, ny, nz = int(P["m"]) + 1, int(P["n"]) + 1, int(P["p"]) + 1
cells = nx * ny * nz
rd = os.path.join(run_dir, "restart_data")

xb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
x, y = 0.5 * (xb[:-1] + xb[1:]), 0.5 * (yb[:-1] + yb[1:])

steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(rd, "lustre_[0-9]*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
step = int(sys.argv[2]) if len(sys.argv) > 2 else steps[-1]

snap = np.fromfile(os.path.join(rd, f"lustre_{step}.dat"), np.float64)
nvars = snap.size // cells


def field(i):
    return snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)


# Temperature is recovered from the mixture stiffened-gas EOS (no T_s scalar); color is last.
gm = [float(P[f"fluid_pp({i})%gamma"]) for i in (1, 2)]  # = 1/(gamma_i - 1)
pin = [float(P[f"fluid_pp({i})%pi_inf"]) for i in (1, 2)]
cvs = [float(P[f"fluid_pp({i})%cv"]) for i in (1, 2)]
gsmin = [1.0 / g + 1.0 for g in gm]
# model_eqns=3 conserved layout: cont(nf), mom(dim), total-E(1), volume fractions(nf),
# phasic internal energies(nf), ..., color(last). Derive the offsets so we don't mis-index the
# total-E slot (a hardcoded 5,6 reads E as a volume fraction -> garbage T + a fake drop "blob").
nf = int(P.get("num_fluids", 2))
dim = 3 if nz > 1 else 2
i_adv = nf + dim + 1  # first volume fraction (after cont, mom, total energy)
i_ie = i_adv + nf  # first phasic internal energy
a1, a2 = field(i_adv), field(i_adv + 1)  # volume fractions
rho_e = field(i_ie) + field(i_ie + 1)  # mixture internal energy (sum of phasic internal energies)
Gamma = a1 * gm[0] + a2 * gm[1]
pi_mix = a1 * pin[0] + a2 * pin[1]
mCP = field(0) * cvs[0] * gsmin[0] + field(1) * cvs[1] * gsmin[1]
T = ((Gamma + 1.0) * (rho_e - pi_mix) / Gamma + pi_mix) / mCP  # EOS temperature
color = np.clip(field(nvars - 1), 0.0, 1.0)
kz = nz // 2  # z midplane
T_mid, c_mid = T[kz], color[kz]  # (ny, nx)

dt = float(P["dt"])
t_r = 7.5
Tw_in, Tw_out = float(P["bc_y%twall_in"]), float(P["bc_y%twall_out"])

fig, (ax_f, ax_l) = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [1.1, 1.0]})
mesh = ax_f.pcolormesh(x, y, T_mid, cmap="coolwarm", shading="auto")
fig.colorbar(mesh, ax=ax_f, label=r"temperature $T$")
ax_f.contour(x, y, c_mid, levels=[0.5], colors="k", linewidths=1.4)
yc = (c_mid * y[:, None]).sum() / c_mid.sum()
xc = (c_mid * x[None, :]).sum() / c_mid.sum()
ax_f.plot(xc, yc, "k+", ms=12, mew=2)
ax_f.set(aspect="equal", xlabel="x", ylabel="y (rise axis)", title=f"$T$ + color, z-midplane (step {step}, $t/t_r$={step * dt / t_r:.3f})\ndrop centroid y={yc:+.4f}")

mid_x = nx // 2
ax_l.plot(T_mid[:, mid_x], y, "-", color="C3", lw=1.8, label="MFC $T$ centerline")
ax_l.plot(float(P["bc_y%twall_in"]) + (float(P["bc_y%twall_out"]) - float(P["bc_y%twall_in"])) * (y - y[0]) / (y[-1] - y[0]), y, ":", color="0.4", lw=1.3, label="imposed linear field")
ax_l.plot([Tw_in, Tw_out], [y[0], y[-1]], "ko", ms=6, mfc="none", mew=1.4, label="isothermal walls")
ax_l.set(xlabel=r"$T$", ylabel="y (rise axis)", title="centerline $T(y)$ vs imposed gradient")
ax_l.legend(loc="best", fontsize=9)
ax_l.grid(alpha=0.3)
fig.tight_layout()

viz = os.path.join(run_dir, "viz")
os.makedirs(viz, exist_ok=True)
out = os.path.join(viz, f"fields_{step}.png")
fig.savefig(out, dpi=150)
print(f"T range [{T.min():.4f}, {T.max():.4f}]  walls [{Tw_in:.4f}, {Tw_out:.4f}]  drop y={yc:+.4f}")
print(f"saved -> {out}")
