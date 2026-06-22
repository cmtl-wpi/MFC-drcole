#!/usr/bin/env python3
"""Checks MFC's heat conduction against the exact 3D Gaussian-hotspot solution.

A Gaussian temperature blob spreads in a cube as the heat-equation Green's
function, T0 + A (sigma0^2/sigma^2)^{3/2} exp(-r^2/2 sigma^2), sigma^2 =
sigma0^2 + 2 alpha t. It must stay spherically symmetric (no directional bias).
This script does NOT run MFC -- run the case first, then run this to analyze it:

  ./mfc.sh run examples/3D_Thermal_Conduction_hotspot/case.py -n 16
  python3 examples/3D_Thermal_Conduction_hotspot/validate.py

It recovers temperature from the EOS, compares to the exact Gaussian, checks the
x/y/z center lines collapse, and writes figures/heat_3d_hotspot.png + summary.json.

How error is reported: L1 = average error, L2 = typical (root-mean-square) error,
Linf = the single worst cell.
"""

import glob
import json
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")

# Physical constants -- mirror case.py
L, ALPHA = 1.0, 0.05
TWALL = 10.0
HS_A, HS_SIG = 5.0, 0.08  # hot spot amplitude, initial width


def exact_3d_hotspot(X, Y, Z, t):
    s2 = HS_SIG**2 + 2.0 * ALPHA * t
    r2 = (X - L / 2) ** 2 + (Y - L / 2) ** 2 + (Z - L / 2) ** 2
    return TWALL + HS_A * (HS_SIG**2 / s2) ** 1.5 * np.exp(-r2 / (2.0 * s2))


def read_run(rundir):
    """Read one run's settings and saved data from rundir.

    Returns: the time step, the grid size, the cell-center coordinates, the list
    of saved step numbers, and a load(step) function that hands back the raw field
    array for that step.
    """
    inp = {}
    for line in open(os.path.join(rundir, "simulation.inp")):
        if "=" in line:
            k, v = line.split("=", 1)
            inp[k.strip().lower()] = v.strip().rstrip(",")
    dt = float(inp["dt"])
    m, n, p = int(inp["m"]) + 1, int(inp.get("n", "0")) + 1, int(inp.get("p", "0")) + 1

    def coord(axis, ncell):
        beg, end = float(inp.get(f"{axis}_domain%beg", "0")), float(inp.get(f"{axis}_domain%end", "1"))
        return beg + (np.arange(ncell) + 0.5) * (end - beg) / ncell

    xc, yc, zc = coord("x", m), coord("y", n), coord("z", p)
    rd = os.path.join(rundir, "restart_data")
    steps = []
    for f in glob.glob(os.path.join(rd, "lustre_*.dat")):
        match = re.search(r"lustre_(\d+)\.dat$", os.path.basename(f))
        if match:
            steps.append(int(match.group(1)))
    steps.sort()
    ncell = m * n * p

    def load(s):
        a = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        return a.reshape(a.size // ncell, p, n, m)  # x innermost

    return dt, (m, n, p), (xc, yc, zc), steps, load


def ndims_of(fields):
    return sum(1 for d in fields.shape[1:] if d > 1)


def temperature(fields):
    """Recover temperature from MFC's raw conserved fields via the stiffened-gas
    EOS: back out the pressure, then T = (p + p_inf) / ((gam - 1)*rho*cv).
    """
    nd = ndims_of(fields)
    gam, p_inf, cv = 2.0, 100.0, 12.5
    rho, E = fields[0], fields[1 + nd]
    ke = 0.5 * sum((fields[1 + i] / rho) ** 2 for i in range(nd))
    pres = (gam - 1.0) * rho * (E / rho - ke) - gam * p_inf
    return ((pres + p_inf) / ((gam - 1.0) * rho * cv)).squeeze()


def velocity_mag(fields):
    """Flow speed |u| in every cell (the spurious flow from weak conduction-driven
    acoustics; should stay near zero)."""
    nd = ndims_of(fields)
    rho = fields[0]
    vel = [fields[1 + i] / rho for i in range(nd)]
    return np.sqrt(sum(u**2 for u in vel)).squeeze()


def norms(num, ex):
    e = np.abs(num - ex)
    return {"L1": float(e.mean()), "L2": float(np.sqrt((e**2).mean())), "Linf": float(e.max())}


def main():
    rd = os.path.join(HERE, "restart_data")
    if not os.path.isdir(rd):
        raise SystemExit("no restart_data/ here -- run first:\n  ./mfc.sh run examples/3D_Thermal_Conduction_hotspot/case.py -n 16")
    dt, (m, n, p), (xc, yc, zc), steps, load = read_run(HERE)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    R = np.sqrt((X - L / 2) ** 2 + (Y - L / 2) ** 2 + (Z - L / 2) ** 2)
    t = np.array([s * dt for s in steps])
    T_last = np.transpose(temperature(load(steps[-1])), (2, 1, 0))
    ex_last = exact_3d_hotspot(X, Y, Z, t[-1])
    umax = max(float(velocity_mag(load(s)).max()) for s in steps)
    nrm = norms(T_last, ex_last)

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1.0])
    kz = p // 2
    for col, (lab, fld) in enumerate([("MFC", T_last), ("analytic", ex_last)]):
        ax = fig.add_subplot(gs[0, col])
        im = ax.pcolormesh(xc, yc, fld[:, :, kz].T, cmap="inferno", shading="auto", vmin=float(ex_last.min()), vmax=float(ex_last.max()))
        ax.set(title=f"{lab}  z=L/2,  t={t[-1]:.3f}", xlabel="x", ylabel="y", aspect="equal")
        if col == 1:
            fig.colorbar(im, ax=ax, shrink=0.8, label="T")
    ic = m // 2
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(xc - L / 2, T_last[:, ic, ic], "C0-", label="along x")
    ax2.plot(yc - L / 2, T_last[ic, :, ic], "C1--", label="along y")
    ax2.plot(zc - L / 2, T_last[ic, ic, :], "C2:", label="along z")
    ax2.set(title="lines along x, y, z lie on top of each other (round spreading)", xlabel="distance from center", ylabel="T")
    ax2.legend(fontsize=9)
    ax3 = fig.add_subplot(gs[1, :2])
    rr = R.ravel()
    sel = rr < 0.45
    ax3.plot(rr[sel][::37], T_last.ravel()[sel][::37], ".", ms=2, color="0.5", label="MFC (cells)")
    rf = np.linspace(0, 0.45, 200)
    s2 = HS_SIG**2 + 2 * ALPHA * t[-1]
    ax3.plot(rf, TWALL + HS_A * (HS_SIG**2 / s2) ** 1.5 * np.exp(-(rf**2) / (2 * s2)), "k-", lw=2, label="analytic Gaussian")
    ax3.set(xlabel="r (from center)", ylabel="T", title=f"radial collapse at t={t[-1]:.3f}")
    ax3.legend(fontsize=9)
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.axis("off")
    ax4.text(0.0, 0.95, f"grid {m}³\nmax|u| = {umax:.1e}\n\nL1 = {nrm['L1']:.3e}\nL2 = {nrm['L2']:.3e}\nL∞ = {nrm['Linf']:.3e}", va="top", family="monospace", fontsize=10)
    fig.suptitle("3D heat: a hot blob spreading out evenly in all directions", fontweight="bold")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heat_3d_hotspot.png"), dpi=130)
    plt.close(fig)
    print(f"  3D hotspot: L2={nrm['L2']:.3e}  L∞={nrm['Linf']:.3e}  max|u|={umax:.1e}")
    json.dump({"heat_3d_hotspot": {"N": int(m), "max_u": umax, "norms": nrm}}, open(SUMMARY, "w"), indent=2)
    print(f"  wrote {os.path.relpath(os.path.join(FIG, 'heat_3d_hotspot.png'), HERE)}, summary.json")


if __name__ == "__main__":
    main()
