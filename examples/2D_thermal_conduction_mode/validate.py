#!/usr/bin/env python3
"""Checks MFC's heat conduction against the exact 2D periodic-mode solution.

A sine wave in a doubly-periodic box, T0 + A sin(kx) sin(ky), k = 2*pi/L, cools in
place as exp(-2 alpha k^2 t). This script does NOT run MFC -- run the case first,
then run this to analyze the result:

  ./mfc.sh run examples/2D_thermal_conduction_mode/case.py -n 8
  python3 examples/2D_thermal_conduction_mode/validate.py

It recovers temperature from the EOS, compares to the exact solution, fits the
in-place cooling rate, and writes figures/heat_2d_mode.png + summary.json.

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

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")

# IC constants -- mirror the analytic patch expression in case.py (everything
# else -- grid, domain, EOS, k_therm -- comes from simulation.inp)
T0, AMP = 10.0, 3.0  # background temperature, mode amplitude


def exact_2d_mode(X, Y, t, kw, alpha):
    rate = alpha * 2.0 * kw**2
    return T0 + AMP * np.sin(kw * X) * np.sin(kw * Y) * math.exp(-rate * t)


def read_run(rundir):
    """Read one run's settings and saved data from rundir.

    Returns: the time step, the grid size, the cell-center coordinates, the list
    of saved step numbers, a load(step) function that hands back the raw field
    array for that step, and the parsed simulation.inp dictionary.
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

    return dt, (m, n, p), (xc, yc, zc), steps, load, inp


def ndims_of(fields):
    return sum(1 for d in fields.shape[1:] if d > 1)


def temperature(fields, g, p_inf, cv):
    """Recover temperature from MFC's raw conserved fields via the stiffened-gas
    EOS: back out the pressure, then T = (p + p_inf) / ((g - 1)*rho*cv).
    """
    nd = ndims_of(fields)
    rho, E = fields[0], fields[1 + nd]
    ke = 0.5 * sum((fields[1 + i] / rho) ** 2 for i in range(nd))
    pres = (g - 1.0) * rho * (E / rho - ke) - g * p_inf
    return ((pres + p_inf) / ((g - 1.0) * rho * cv)).squeeze()


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
        raise SystemExit("no restart_data/ here -- run first:\n  ./mfc.sh run examples/2D_thermal_conduction_mode/case.py -n 8")
    dt, (m, n, _), (xc, yc, _), steps, load, inp = read_run(HERE)
    # EOS + conduction constants from simulation.inp (MFC stores gamma as 1/(g-1)
    # and pi_inf as g*p_inf/(g-1), so physical p_inf = pi_inf_stored*(g-1)/g)
    g = 1.0 + 1.0 / float(inp["fluid_pp(1)%gamma"])
    p_inf = float(inp["fluid_pp(1)%pi_inf"]) * (g - 1.0) / g
    cv = float(inp["fluid_pp(1)%cv"])
    k_therm = float(inp["fluid_pp(1)%k_therm"])
    L = float(inp["x_domain%end"]) - float(inp["x_domain%beg"])
    kw = 2.0 * math.pi / L  # periodic single-mode wavenumber
    # background density from the step-0 fields: p is uniform at t=0, so rho*T is
    # uniform and equals rho0*T0 (the analytic-patch normalization in case.py)
    f0 = load(steps[0])
    rho0 = float((f0[0].squeeze() * temperature(f0, g, p_inf, cv)).mean()) / T0
    alpha = k_therm / (rho0 * g * cv)  # decay-rate diffusivity k/(rho*cp), cp = g*cv
    X, Y = np.meshgrid(xc, yc, indexing="ij")  # [m,n]
    t = np.array([s * dt for s in steps])
    rate_an = alpha * 2.0 * kw**2
    L1s, L2s, Lis, amp, umax = [], [], [], [], 0.0
    field_last = ex_last = None
    # per step: error norms vs analytic, plus the mode amplitude (projection of
    # (T - T0) onto sin*sin) whose exponential decay gives the cooling rate.
    for i, s in enumerate(steps):
        fields = load(s)
        T = np.transpose(temperature(fields, g, p_inf, cv), (1, 0))  # [y,x]->[x,y]
        ex = exact_2d_mode(X, Y, t[i], kw, alpha)
        nr = norms(T, ex)
        L1s.append(nr["L1"])
        L2s.append(nr["L2"])
        Lis.append(nr["Linf"])
        amp.append(((T - T0) * np.sin(kw * X) * np.sin(kw * Y)).sum())
        umax = max(umax, float(velocity_mag(fields).max()))
        if i == len(steps) - 1:
            field_last, ex_last = T, ex
    L1s, L2s, Lis, amp = map(np.array, (L1s, L2s, Lis, amp))
    rate = -np.polyfit(t, np.log(np.abs(amp / amp[0])), 1)[0]
    rate_err = 100 * abs(rate / rate_an - 1)

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1.0])
    vmin, vmax = float(ex_last.min()), float(ex_last.max())
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.pcolormesh(xc, yc, field_last.T, vmin=vmin, vmax=vmax, cmap="inferno", shading="auto")
    ax0.set(title=f"MFC  T(x,y),  t={t[-1]:.3f}", xlabel="x", ylabel="y", aspect="equal")
    ax1 = fig.add_subplot(gs[0, 1])
    im = ax1.pcolormesh(xc, yc, ex_last.T, vmin=vmin, vmax=vmax, cmap="inferno", shading="auto")
    ax1.set(title="analytic", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(im, ax=[ax0, ax1], shrink=0.8, label="T")
    jy = n // 4  # probe the antinode (y=L/4); y=L/2 is a node of sin(ky)
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(xc, ex_last[:, jy], "k-", lw=2, label="analytic")
    ax2.plot(xc[::3], field_last[::3, jy], "C3o", ms=4, mfc="none", label="MFC")
    ax2.set(title="line  y=L/4 (antinode)", xlabel="x", ylabel="T")
    ax2.legend(fontsize=9)
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.semilogy(t, L1s, "C0-o", ms=3, label="L1")
    ax3.semilogy(t, L2s, "C1-s", ms=3, label="L2")
    ax3.semilogy(t, Lis, "C3-^", ms=3, label="L∞")
    ax3.set(xlabel="t", ylabel="error norm", title="L1/L2/L∞ vs t")
    ax3.legend(fontsize=9)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(t, amp / amp[0], "C0-", lw=2, label="MFC")
    ax4.plot(t, np.exp(-rate_an * t), "k--", lw=1.5, label=f"e^(−{rate_an:.2f}t)")
    ax4.set(xlabel="t", ylabel="mode amplitude", title=f"cooling rate: measured {rate:.3f} vs exact {rate_an:.3f} ({rate_err:.1f}%)")
    ax4.legend(fontsize=9)
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis("off")
    ax5.text(0.0, 0.95, f"grid {m}×{n}\nmax|u| = {umax:.1e}\n\npeak L∞ = {Lis.max():.3e}\nfinal L2 = {L2s[-1]:.3e}\n\nrate err = {rate_err:.2f}%", va="top", family="monospace", fontsize=10)
    fig.suptitle("2D heat: a sine wave in a periodic box, cooling in place", fontweight="bold")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heat_2d_mode.png"), dpi=130)
    plt.close(fig)
    print(f"  2D mode: peak L∞={Lis.max():.3e}  final L2={L2s[-1]:.3e}  rate {rate:.3f} vs {rate_an:.3f} ({rate_err:.2f}%)  max|u|={umax:.1e}")
    out = {"heat_2d_mode": {"N": int(m), "max_u": umax, "measured_rate": rate, "analytic_rate": rate_an, "rate_error_pct": rate_err, "Linf_peak": float(Lis.max()), "L2_final": float(L2s[-1])}}
    open(SUMMARY, "w").write(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {os.path.relpath(os.path.join(FIG, 'heat_2d_mode.png'), HERE)}, summary.json")


if __name__ == "__main__":
    main()
