#!/usr/bin/env python3
"""Checks MFC's heat conduction against the exact 1D Dirichlet sine-decay solution.

A sine-shaped temperature bump cools between fixed-temperature walls. The exact
answer is T_wall + A sin(pi x/L) exp(-alpha (pi/L)^2 t). This script does NOT run
MFC -- run the case first, then run this to analyze the result:

  ./mfc.sh run examples/Thermal_Conduction_1D/case.py -n 4
  python3 examples/Thermal_Conduction_1D/validate.py

It reads back the temperature MFC saved (restart_data/ in this directory), compares
to the exact solution, reports the error, and writes figures/heat_1d.png + summary.json.

How error is reported: L1 = average error, L2 = typical (root-mean-square) error,
Linf = the single worst cell.
"""

import glob
import json
import math
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")

# Physical constants -- mirror case.py
L, ALPHA = 1.0, 0.05
TWALL, AMP = 10.0, 3.0  # 1D Dirichlet sine: Twall + AMP*sin(pi x/L)
KAPPA = math.pi / L


def exact_1d(x, t):
    return TWALL + AMP * np.sin(KAPPA * x) * math.exp(-ALPHA * KAPPA**2 * t)


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
    # MFC saves each step as a file named lustre_<step>.dat -- collect the step numbers
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
    EOS. MFC carried density and energy, so back out the pressure and then
    T = (p + p_inf) / ((gam - 1)*rho*cv).
    """
    nd = ndims_of(fields)
    gam, p_inf, cv = 2.0, 100.0, 12.5
    rho, E = fields[0], fields[1 + nd]
    ke = 0.5 * sum((fields[1 + i] / rho) ** 2 for i in range(nd))
    pres = (gam - 1.0) * rho * (E / rho - ke) - gam * p_inf
    return ((pres + p_inf) / ((gam - 1.0) * rho * cv)).squeeze()


def velocity_mag(fields):
    """Flow speed |u| in every cell. MFC stores momentum (rho*u), not velocity, so
    divide each momentum component by density. In a pure-conduction run this should
    stay near zero -- it is the spurious flow from the weak conduction-driven acoustics.
    """
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
        raise SystemExit("no restart_data/ here -- run first:\n  ./mfc.sh run examples/Thermal_Conduction_1D/case.py -n 4")
    dt, (m, _, _), (xc, _, _), steps, load = read_run(HERE)
    t = np.array([s * dt for s in steps])
    nums = [temperature(load(s)) for s in steps]
    umax = max(float(np.atleast_1d(velocity_mag(load(s))).max()) for s in steps)
    nrm = [norms(nu, exact_1d(xc, ti)) for nu, ti in zip(nums, t)]
    L2 = np.array([d["L2"] for d in nrm])
    Li = np.array([d["Linf"] for d in nrm])
    out = {
        "N": int(m),
        "max_u": umax,
        "L1_final": float(nrm[-1]["L1"]),
        "L2_final": float(L2[-1]),
        "Linf_final": float(Li[-1]),
        "Linf_peak": float(Li.max()),
        "amplitude": AMP,
        "n_outputs": len(steps),
    }
    print(f"  1D: peak L∞={Li.max():.3e}  final L2={L2[-1]:.3e}  max|u|={umax:.2e}")

    fig, (axT, axE) = plt.subplots(1, 2, figsize=(12, 5))
    # left: MFC T(x) profiles (markers) on the analytic curves (lines), colored by time
    idx = np.unique(np.linspace(0, len(t) - 1, 4).round().astype(int))
    tcol = plt.cm.viridis(np.linspace(0.12, 0.88, len(idx)))
    xf = np.linspace(0, L, 400)
    for c, j in zip(tcol, idx):
        axT.plot(xf, exact_1d(xf, t[j]), "-", color=c, lw=1.6, zorder=1)
        axT.plot(xc[::8], nums[j][::8], "o", color=c, ms=4, mfc="none", mew=1.1, zorder=2)
    axT.legend(
        [Line2D([], [], color="0.3", lw=1.6), Line2D([], [], color="0.3", marker="o", ls="none", mfc="none")],
        ["exact", "MFC (EOS T)"],
        loc="upper right",
        fontsize=9,
    )
    axT.set(xlabel="x", ylabel="T", title="temperature T(x) at t = " + ", ".join(f"{t[j]:.2f}" for j in idx))
    # right: error norms over time
    axE.semilogy(t, L2, "C3-", lw=1.8, label="L2")
    axE.semilogy(t, Li, "C3--", lw=1.3, label="L∞")
    axE.set(xlabel="t", ylabel="error size", title=f"error over time  (worst L∞ {Li.max():.2e})")
    axE.legend(fontsize=9)
    fig.suptitle("1D heat: a sine bump cooling between fixed-temperature walls (EOS temperature)", fontweight="bold")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heat_1d.png"), dpi=130)
    plt.close(fig)
    json.dump({"heat_1d": out}, open(SUMMARY, "w"), indent=2)
    print(f"  wrote {os.path.relpath(os.path.join(FIG, 'heat_1d.png'), HERE)}, summary.json")


if __name__ == "__main__":
    main()
