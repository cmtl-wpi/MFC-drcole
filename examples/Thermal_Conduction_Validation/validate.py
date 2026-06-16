#!/usr/bin/env python3
"""Checks MFC's heat conduction against known textbook answers.

This script does NOT run MFC -- run the cases first, then run this to analyze
the results. For each benchmark it reads back the temperature MFC saved (from
runs/<label>/), compares it to the exact pen-and-paper solution of the heat
equation, reports how big the error is, and saves a figure plus a line in
summary.json.

Each benchmark reads an archived run under runs/<label>/, which must already hold
that case's restart_data/ and simulation.inp. To produce one, run the matching
case with MFC and copy its restart_data/ + simulation.inp into runs/<label>/:

  benchmark    reads runs/<label>          from case
  1d           1d                          case_1d.py
  2d           2d_mode                     case_2d_mode.py
  3d-mode      3d_mode                     case_3d_mode.py
  3d-hotspot   3d_hotspot                  case_3d_hotspot.py
  conv-x       convx_{32,64,128,256,512}   case_conv.py (one run per grid)
  conv-t       convt_{50,100,...,800,3200} case_conv.py (one run per step count)

  python3 examples/Thermal_Conduction_Validation/validate.py 1d          # 1D bar, fixed-temperature ends
  python3 examples/Thermal_Conduction_Validation/validate.py 2d          # 2D box, a sine wave cooling off
  python3 examples/Thermal_Conduction_Validation/validate.py 3d-mode     # 3D box, a sine wave cooling off
  python3 examples/Thermal_Conduction_Validation/validate.py 3d-hotspot  # 3D box, a hot blob spreading out
  python3 examples/Thermal_Conduction_Validation/validate.py conv-x      # grid convergence study (cell size)
  python3 examples/Thermal_Conduction_Validation/validate.py conv-t      # time-step convergence study

How error is reported: L1 = average error, L2 = typical (root-mean-square)
error, Linf = the single worst cell.
"""

import glob
import json
import math
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")

# Physical constants -- mirror cases/case_*.py
L, ALPHA = 1.0, 0.05
TWALL, AMP = 10.0, 3.0  # 1D Dirichlet sine: Twall + AMP*sin(pi x/L)
KAPPA = math.pi / L
KW = 2.0 * math.pi / L  # periodic single-mode wavenumber (conv, 2D mode, 3D mode)
HS_A, HS_SIG = 5.0, 0.08  # 3D hot spot amplitude, initial width


# Exact textbook solutions of the heat equation, to compare MFC against


def exact_1d(x, t):
    return TWALL + AMP * np.sin(KAPPA * x) * math.exp(-ALPHA * KAPPA**2 * t)


def exact_conv(x, t):
    return TWALL + AMP * np.sin(KW * x) * math.exp(-ALPHA * KW**2 * t)


def exact_2d_mode(X, Y, t):
    rate = ALPHA * 2.0 * KW**2
    return TWALL + AMP * np.sin(KW * X) * np.sin(KW * Y) * math.exp(-rate * t)


def exact_3d_mode(X, Y, Z, t):
    rate = ALPHA * 3.0 * KW**2
    return TWALL + AMP * np.sin(KW * X) * np.sin(KW * Y) * np.sin(KW * Z) * math.exp(-rate * t)


def exact_3d_hotspot(X, Y, Z, t):
    s2 = HS_SIG**2 + 2.0 * ALPHA * t
    r2 = (X - L / 2) ** 2 + (Y - L / 2) ** 2 + (Z - L / 2) ** 2
    return TWALL + HS_A * (HS_SIG**2 / s2) ** 1.5 * np.exp(-r2 / (2.0 * s2))


# Helpers: locate an already-run case, then read back what MFC saved


def find_run(label, case):
    """Return the archived run directory runs/<label>, erroring if MFC has not run it.

    This script does not run MFC. Populate runs/<label>/ first by running the
    matching case and copying its restart_data/ + simulation.inp there.
    """
    rundir = os.path.join(RUNS, label)
    if not os.path.isdir(os.path.join(rundir, "restart_data")):
        raise SystemExit(f"missing {rundir}/restart_data\n  run cases/{case} with MFC, then copy its restart_data/ + simulation.inp into {rundir}/")
    return rundir


def read_run(rundir):
    """Read one run's settings and saved data.

    Returns: the time step, the grid size, the cell-center coordinates, the list
    of saved step numbers, and a load(step) function that hands back the raw
    field array for that step.
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
    """Flow speed |u| in every cell. MFC stores momentum (rho*u), not velocity,
    so divide each momentum component by density to recover the velocity, then
    take the magnitude. In a pure-conduction run this should stay near zero.
    """
    nd = ndims_of(fields)
    rho = fields[0]
    momentum = [fields[1 + i] for i in range(nd)]  # one component per spatial dim
    vel = [mom / rho for mom in momentum]
    speed_sq = sum(u**2 for u in vel)
    return np.sqrt(speed_sq).squeeze()


def norms(num, ex):
    e = np.abs(num - ex)
    return {"L1": float(e.mean()), "L2": float(np.sqrt((e**2).mean())), "Linf": float(e.max())}


def save_summary(key, entry):
    summary = json.load(open(SUMMARY)) if os.path.exists(SUMMARY) else {}
    summary[key] = entry
    json.dump(summary, open(SUMMARY, "w"), indent=2)


def slope(xs, ys):
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


# 1D test: a sine-shaped temperature bump cooling between fixed-temperature walls


def bench_1d():
    rundir = find_run("1d", "case_1d.py")
    dt, (m, _, _), (xc, _, _), steps, load = read_run(rundir)
    t = np.array([s * dt for s in steps])
    nums = [temperature(load(s)) for s in steps]
    umax = max(float(np.atleast_1d(velocity_mag(load(s))).max()) for s in steps)
    nrm = [norms(nu, exact_1d(xc, ti)) for nu, ti in zip(nums, t)]
    L1 = np.array([d["L1"] for d in nrm])
    L2 = np.array([d["L2"] for d in nrm])
    Li = np.array([d["Linf"] for d in nrm])
    out = {
        "N": int(m),
        "max_u": umax,
        "L1_final": float(L1[-1]),
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
        # plot analytical solution
        axT.plot(xf, exact_1d(xf, t[j]), "-", color=c, lw=1.6, zorder=1)
        # plot MFC results
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
    save_summary("heat_1d", out)


# 2D test: a sine wave in a periodic box that should cool off evenly in place


def bench_2d():
    rundir = find_run("2d_mode", "case_2d_mode.py")
    dt, (m, n, _), (xc, yc, _), steps, load = read_run(rundir)
    X, Y = np.meshgrid(xc, yc, indexing="ij")  # [m,n]
    t = np.array([s * dt for s in steps])
    rate_an = ALPHA * 2.0 * KW**2
    L1s, L2s, Lis, amp, umax = [], [], [], [], 0.0
    field_last = ex_last = None
    # per step: error norms vs analytic, plus the mode amplitude (projection of
    # (T - T0) onto sin*sin) whose exponential decay gives the cooling rate.
    for i, s in enumerate(steps):
        fields = load(s)
        T = np.transpose(temperature(fields), (1, 0))  # [y,x]->[x,y]
        ex = exact_2d_mode(X, Y, t[i])
        nr = norms(T, ex)
        L1s.append(nr["L1"])
        L2s.append(nr["L2"])
        Lis.append(nr["Linf"])
        amp.append(((T - TWALL) * np.sin(KW * X) * np.sin(KW * Y)).sum())
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
    save_summary("heat_2d_mode", {"N": int(m), "max_u": umax, "measured_rate": rate, "analytic_rate": rate_an, "rate_error_pct": rate_err, "Linf_peak": float(Lis.max()), "L2_final": float(L2s[-1])})


# 3D test: a sine wave in a box that should cool off evenly in every direction


def bench_3d_mode():
    rundir = find_run("3d_mode", "case_3d_mode.py")
    dt, (m, n, p), (xc, yc, zc), steps, load = read_run(rundir)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    t = np.array([s * dt for s in steps])
    rate_an = ALPHA * 3.0 * KW**2
    L1s, L2s, Lis, amp, umax = [], [], [], [], 0.0
    field_last = ex_last = None
    # per step: error norms vs analytic, plus the mode amplitude (projection of
    # (T - T0) onto sin*sin*sin) whose exponential decay gives the cooling rate.
    for i, s in enumerate(steps):
        fields = load(s)
        T = np.transpose(temperature(fields), (2, 1, 0))  # [z,y,x]->[x,y,z]
        ex = exact_3d_mode(X, Y, Z, t[i])
        nr = norms(T, ex)
        L1s.append(nr["L1"])
        L2s.append(nr["L2"])
        Lis.append(nr["Linf"])
        amp.append(((T - TWALL) * np.sin(KW * X) * np.sin(KW * Y) * np.sin(KW * Z)).sum())
        umax = max(umax, float(velocity_mag(fields).max()))
        if i == len(steps) - 1:
            field_last, ex_last = T, ex
    L1s, L2s, Lis, amp = map(np.array, (L1s, L2s, Lis, amp))
    rate = -np.polyfit(t, np.log(np.abs(amp / amp[0])), 1)[0]
    rate_err = 100 * abs(rate / rate_an - 1)

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1.0])
    kz = p // 2
    vmin, vmax = float(ex_last[:, :, kz].min()), float(ex_last[:, :, kz].max())
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.pcolormesh(xc, yc, field_last[:, :, kz].T, vmin=vmin, vmax=vmax, cmap="inferno", shading="auto")
    ax0.set(title=f"MFC  T(x,y, z=L/2),  t={t[-1]:.3f}", xlabel="x", ylabel="y", aspect="equal")
    ax1 = fig.add_subplot(gs[0, 1])
    im = ax1.pcolormesh(xc, yc, ex_last[:, :, kz].T, vmin=vmin, vmax=vmax, cmap="inferno", shading="auto")
    ax1.set(title="analytic", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(im, ax=[ax0, ax1], shrink=0.8, label="T")
    jy = n // 2
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(xc, ex_last[:, jy, kz], "k-", lw=2, label="analytic")
    ax2.plot(xc[::3], field_last[::3, jy, kz], "C3o", ms=4, mfc="none", label="MFC")
    ax2.set(title="center line  y=z=L/2", xlabel="x", ylabel="T")
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
    ax4.set(xlabel="t", ylabel="wave height", title=f"cooling rate: measured {rate:.3f} vs exact {rate_an:.3f} ({rate_err:.1f}%)")
    ax4.legend(fontsize=9)
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.axis("off")
    ax5.text(0.0, 0.95, f"grid {m}³\nmax|u| = {umax:.1e}\n\npeak L∞ = {Lis.max():.3e}\nfinal L2 = {L2s[-1]:.3e}\n\nrate err = {rate_err:.2f}%", va="top", family="monospace", fontsize=10)
    fig.suptitle("3D heat: a sine wave in a box, cooling evenly in every direction", fontweight="bold")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heat_3d_mode.png"), dpi=130)
    plt.close(fig)
    print(f"  3D mode: peak L∞={Lis.max():.3e}  final L2={L2s[-1]:.3e}  rate {rate:.3f} vs {rate_an:.3f} ({rate_err:.2f}%)  max|u|={umax:.1e}")
    save_summary("heat_3d_mode", {"N": int(m), "max_u": umax, "measured_rate": rate, "analytic_rate": rate_an, "rate_error_pct": rate_err, "Linf_peak": float(Lis.max()), "L2_final": float(L2s[-1])})


# 3D test: a hot blob spreading out (staying round) in a box


def bench_3d_hotspot():
    rundir = find_run("3d_hotspot", "case_3d_hotspot.py")
    dt, (m, n, p), (xc, yc, zc), steps, load = read_run(rundir)
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
    save_summary("heat_3d_hotspot", {"N": int(m), "max_u": umax, "norms": nrm})


# Convergence studies: check the error shrinks at the rate the math predicts


def converge_spatial():
    grids = [32, 64, 128, 256, 512]
    tstar = 0.3 / (ALPHA * KW**2)
    dxs, L2s, Lis = [], [], []
    for N in grids:
        rundir = find_run(f"convx_{N}", "case_conv.py")
        dtr, _, (xc, _, _), steps, load = read_run(rundir)
        nr = norms(temperature(load(steps[-1])), exact_conv(xc, steps[-1] * dtr))
        dxs.append(L / N)
        L2s.append(nr["L2"])
        Lis.append(nr["Linf"])
        print(f"  N={N:4d}  dx={L / N:.4e}  L2={nr['L2']:.4e}  L∞={nr['Linf']:.4e}")
    dxs, L2s, Lis = map(np.array, (dxs, L2s, Lis))
    s2 = slope(dxs, L2s)
    si = slope(dxs, Lis)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.loglog(dxs, L2s, "C0-o", label=f"L2  (slope {s2:.2f})")
    ax.loglog(dxs, Lis, "C3-^", label=f"L∞  (slope {si:.2f})")
    ax.loglog(dxs, L2s[0] * (dxs / dxs[0]) ** 2, "k--", lw=1.3, label="formal slope 2 (operator)")
    ax.set(xlabel="cell size Δx", ylabel="error vs analytic", title=f"Spatial convergence: plateaus at the physics floor (slope {s2:.2f})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "convergence_spatial.png"), dpi=140)
    plt.close(fig)
    print(f"  GRID convergence slope: L2={s2:.3f}  L∞={si:.3f}  (plateaus at the physics floor; operator order is 2)")
    save_summary("convergence_spatial", {"grids": grids, "dx": dxs.tolist(), "L2": L2s.tolist(), "Linf": Lis.tolist(), "slope_L2": s2, "slope_Linf": si, "tstar": tstar})


def converge_temporal():
    # Use a coarse grid on purpose. We compare each run against a very-small-step
    # reference run on the SAME grid, so the cell-size error AND the energy-path
    # physics floor (variable diffusivity, acoustics) cancel out and only the
    # time-step error is left -- which is why this still isolates RK3 order even
    # though the spatial study plateaus. Step counts start at 256: on N=32 the
    # acoustic CFL caps dt, so coarser sweeps would be unstable.
    N = 32
    tstar = 0.3 / (ALPHA * KW**2)
    step_counts = [256, 512, 1024, 2048]
    ref_steps = 4096
    fields, dt_of = {}, {}
    for ns in step_counts + [ref_steps]:
        rundir = find_run(f"convt_{ns}", "case_conv.py")
        dtr, _, (xc, _, _), steps, load = read_run(rundir)
        fields[ns] = temperature(load(steps[-1]))
        dt_of[ns] = dtr
        print(f"  steps={ns:5d}  dt={dtr:.3e}")
    ref = fields[ref_steps]
    dts, L2s = [], []
    for ns in step_counts:
        e = np.abs(fields[ns] - ref)
        dts.append(dt_of[ns])
        L2s.append(float(np.sqrt((e**2).mean())))
    dts, L2s = np.array(dts), np.array(L2s)
    st = slope(dts, L2s)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.loglog(dts, L2s, "C0-o", label=f"L2 vs reference  (slope {st:.2f})")
    ax.loglog(dts, L2s[0] * (dts / dts[0]) ** 3, "k--", lw=1.3, label="reference slope 3")
    ax.set(xlabel="time step Δt", ylabel="error vs tiny-step reference", title=f"Time-step convergence study (grid fixed at N={N})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "convergence_temporal.png"), dpi=140)
    plt.close(fig)
    print(f"  TIME-STEP convergence slope: L2={st:.3f}  (should be ~3)")
    save_summary("convergence_temporal", {"N": N, "dt": dts.tolist(), "L2": L2s.tolist(), "slope_L2": st, "ref_steps": ref_steps, "tstar": tstar})


COMMANDS = {"1d": bench_1d, "2d": bench_2d, "3d-mode": bench_3d_mode, "3d-hotspot": bench_3d_hotspot, "conv-x": converge_spatial, "conv-t": converge_temporal}

if __name__ == "__main__":
    cmds = sys.argv[1:]
    if not cmds or cmds[0] not in COMMANDS:
        print("usage: validate.py {" + "|".join(COMMANDS) + "}  (run MFC first; this only analyzes runs/)")
        sys.exit(1)
    COMMANDS[cmds[0]]()
