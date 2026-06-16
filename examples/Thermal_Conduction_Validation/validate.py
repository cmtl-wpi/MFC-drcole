#!/usr/bin/env python3
"""Checks MFC's heat conduction against known textbook answers.

For each test it runs MFC, reads the temperature MFC computed, compares it to
the exact pen-and-paper solution of the heat equation, reports how big the
error is, and saves a figure plus a line in summary.json.

Run from the repo root (needs mpirun, so disable the command sandbox):
  python3 examples/Thermal_Conduction_Validation/validate.py 1d          # 1D bar, fixed-temperature ends
  python3 examples/Thermal_Conduction_Validation/validate.py 2d          # 2D plate, one hot edge
  python3 examples/Thermal_Conduction_Validation/validate.py 3d-mode     # 3D box, a sine wave cooling off
  python3 examples/Thermal_Conduction_Validation/validate.py 3d-hotspot  # 3D box, a hot blob spreading out
  python3 examples/Thermal_Conduction_Validation/validate.py conv-x      # grid convergence study (cell size)
  python3 examples/Thermal_Conduction_Validation/validate.py conv-t      # time-step convergence study

Add --no-run to re-use the saved runs/ instead of simulating again.

How error is reported: L1 = average error, L2 = typical (root-mean-square)
error, Linf = the single worst cell.
"""

import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASES_REL = os.path.join("examples", "Thermal_Conduction_Validation", "cases")
CASES = os.path.join(HERE, "cases")
RUNS = os.path.join(HERE, "runs")
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")
NO_RUN = "--no-run" in sys.argv

# Physical constants -- mirror cases/case_*.py
L, ALPHA = 1.0, 0.05
TWALL, AMP = 10.0, 3.0  # 1D Dirichlet sine: Twall + AMP*sin(pi x/L)
KAPPA = math.pi / L
KW = 2.0 * math.pi / L  # periodic single-mode wavenumber (conv, 3D mode)
TCOLD, THOT = 10.0, 110.0  # 2D plate: cold edges / top edge (ΔT=100, MFC requires Twall>0)
HS_A, HS_SIG = 5.0, 0.08  # 3D hot spot amplitude, initial width


# Exact textbook solutions of the heat equation, to compare MFC against


def exact_1d(x, t):
    return TWALL + AMP * np.sin(KAPPA * x) * math.exp(-ALPHA * KAPPA**2 * t)


def exact_conv(x, t):
    return TWALL + AMP * np.sin(KW * x) * math.exp(-ALPHA * KW**2 * t)


def exact_2d_plate(X, Y, n_terms=400):
    T = np.zeros_like(X)
    for n in range(1, 2 * n_terms, 2):  # odd n
        kn = n * np.pi / L
        # sinh(kn*Y)/sinh(kn*L) written so every exponent is <= 0 (sinh(kn*L) itself overflows for large n)
        ratio = np.exp(kn * (Y - L)) * (1.0 - np.exp(-2.0 * kn * Y)) / (1.0 - np.exp(-2.0 * kn * L))
        T += (4.0 * (THOT - TCOLD) / (n * np.pi)) * np.sin(kn * X) * ratio
    return TCOLD + T


def exact_3d_mode(X, Y, Z, t):
    rate = ALPHA * 3.0 * KW**2
    return TWALL + AMP * np.sin(KW * X) * np.sin(KW * Y) * np.sin(KW * Z) * math.exp(-rate * t)


def exact_3d_hotspot(X, Y, Z, t):
    s2 = HS_SIG**2 + 2.0 * ALPHA * t
    r2 = (X - L / 2) ** 2 + (Y - L / 2) ** 2 + (Z - L / 2) ** 2
    return TWALL + HS_A * (HS_SIG**2 / s2) ** 1.5 * np.exp(-r2 / (2.0 * s2))


# Helpers: run a case, then read back what MFC saved


def run_case(case, env, ranks, label):
    rundir = os.path.join(RUNS, label)
    if NO_RUN:
        print(f"  [--no-run] {rundir}")
        return rundir
    os.makedirs(RUNS, exist_ok=True)
    full_env = dict(os.environ, **{k: str(v) for k, v in env.items()})
    cmd = ["./mfc.sh", "run", os.path.join(CASES_REL, case), "-n", str(ranks), "-t", "pre_process", "simulation"]
    print(f"  running {case}  ranks={ranks}" + (f"  env={env}" if env else ""))
    res = subprocess.run(cmd, cwd=REPO, env=full_env, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        print(res.stdout[-3000:], res.stderr[-2000:])
        raise RuntimeError(f"MFC run failed: {case} ({label})")
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    shutil.copytree(os.path.join(CASES, "restart_data"), os.path.join(rundir, "restart_data"))
    shutil.copy(os.path.join(CASES, "simulation.inp"), rundir)
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


def temperature(fields, mode):
    """Pull temperature out of MFC's raw field array.

    'scalar' mode: MFC carried temperature directly, so just read it.
    'energy' mode: MFC carried density and energy, so recover temperature from
    the gas law (pressure and density).
    """
    nd = ndims_of(fields)
    if mode == "scalar":
        return fields[nd + 3].squeeze()
    gam, p_inf, cv = 2.0, 100.0, 12.5
    rho, E = fields[0], fields[1 + nd]
    ke = 0.5 * sum((fields[1 + i] / rho) ** 2 for i in range(nd))
    pres = (gam - 1.0) * rho * (E / rho - ke) - gam * p_inf
    return ((pres + p_inf) / ((gam - 1.0) * rho * cv)).squeeze()


def velocity_mag(fields):
    nd = ndims_of(fields)
    return np.sqrt(sum((fields[1 + i] / fields[0]) ** 2 for i in range(nd))).squeeze()


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
    modes = [("scalar", "case_1d.py", "C0", "o"), ("energy", "case_1d_energy.py", "C3", "s")]
    out, data = {}, {}
    for mode, case, color, marker in modes:
        rundir = run_case(case, {}, 1, f"1d_{mode}")
        dt, (m, _, _), (xc, _, _), steps, load = read_run(rundir)
        t = np.array([s * dt for s in steps])
        nums = [temperature(load(s), mode) for s in steps]
        umax = max(float(np.atleast_1d(velocity_mag(load(s))).max()) for s in steps)
        nrm = [norms(nu, exact_1d(xc, ti)) for nu, ti in zip(nums, t)]
        L1 = np.array([d["L1"] for d in nrm])
        L2 = np.array([d["L2"] for d in nrm])
        Li = np.array([d["Linf"] for d in nrm])
        data[mode] = dict(t=t, xc=xc, nums=nums, L2=L2, Li=Li, color=color, marker=marker)
        out[mode] = {
            "N": int(m),
            "max_u": umax,
            "L1_final": float(L1[-1]),
            "L2_final": float(L2[-1]),
            "Linf_final": float(Li[-1]),
            "Linf_peak": float(Li.max()),
            "amplitude": AMP,
            "n_outputs": len(steps),
        }
        print(f"  1D {mode}: peak L∞={Li.max():.3e}  final L2={L2[-1]:.3e}  max|u|={umax:.2e}")

    fig, (axT, axE) = plt.subplots(1, 2, figsize=(12, 5))
    # left: both modes' T(x) profiles on the shared analytic curves (marker distinguishes the mode)
    t = data["scalar"]["t"]
    idx = np.unique(np.linspace(0, len(t) - 1, 4).round().astype(int))
    tcol = plt.cm.viridis(np.linspace(0.12, 0.88, len(idx)))
    xf = np.linspace(0, L, 400)
    for c, j in zip(tcol, idx):
        axT.plot(xf, exact_1d(xf, t[j]), "-", color=c, lw=1.6, zorder=1)
    for mode, _, _, marker in modes:
        d = data[mode]
        for c, j in zip(tcol, idx):
            axT.plot(d["xc"][::8], d["nums"][j][::8], marker, color=c, ms=4, mfc="none", mew=1.1, zorder=2)
    axT.legend(
        [Line2D([], [], color="0.3", lw=1.6), Line2D([], [], color="0.3", marker="o", ls="none", mfc="none"), Line2D([], [], color="0.3", marker="s", ls="none", mfc="none")],
        ["exact", "MFC (direct T)", "MFC (gas-law T)"],
        loc="upper right",
        fontsize=9,
    )
    axT.set(xlabel="x", ylabel="T", title="temperature T(x) at t = " + ", ".join(f"{t[j]:.2f}" for j in idx))
    # right: error norms for both modes (color = mode, solid L2 / dashed L∞)
    for mode, _, color, _ in modes:
        d = data[mode]
        axE.semilogy(d["t"], d["L2"], "-", color=color, lw=1.8, label=f"{mode} L2")
        axE.semilogy(d["t"], d["Li"], "--", color=color, lw=1.3, label=f"{mode} L∞")
    axE.set(xlabel="t", ylabel="error size", title=f"error over time: direct-T vs gas-law  (worst {data['scalar']['Li'].max():.1e} vs {data['energy']['Li'].max():.1e})")
    axE.legend(fontsize=8, ncol=2)
    fig.suptitle("1D heat: a sine bump cooling between fixed-temperature walls", fontweight="bold")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heat_1d.png"), dpi=130)
    plt.close(fig)
    save_summary("heat_1d", out)


# 2D test: a plate with one hot edge, settling to a steady (unchanging) state


def bench_2d():
    rundir = run_case("case_2d_plate.py", {}, 8, "2d_plate")
    dt, (m, n, _), (xc, yc, _), steps, load = read_run(rundir)
    X, Y = np.meshgrid(xc, yc)  # [n,m]
    Tex = exact_2d_plate(X, Y)
    Tnum = temperature(load(steps[-1]), "scalar")
    drift = float(np.abs(Tnum - temperature(load(steps[-2]), "scalar")).max())
    umax = float(velocity_mag(load(steps[-1])).max())
    nrm = norms(Tnum, Tex)
    interior = (Y > 0.05) & (Y < 0.95) & (X > 0.05) & (X < 0.95)
    nrm_int = norms(Tnum[interior], Tex[interior])

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1.0])
    lv = np.linspace(TCOLD, THOT, 21)
    ax0 = fig.add_subplot(gs[0, 0])
    cf = ax0.contourf(X, Y, Tnum, levels=lv, cmap="inferno")
    ax0.set(title="MFC (steady)", xlabel="x", ylabel="y", aspect="equal")
    ax1 = fig.add_subplot(gs[0, 1])
    ax1.contourf(X, Y, Tex, levels=lv, cmap="inferno")
    ax1.set(title="exact (textbook series)", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(cf, ax=[ax0, ax1], shrink=0.8, label="T")
    ax2 = fig.add_subplot(gs[0, 2])
    d = Tnum - Tex
    dm = float(np.abs(d).max()) or 1e-12
    im = ax2.pcolormesh(X, Y, d, vmin=-dm, vmax=dm, cmap="coolwarm", shading="auto")
    ax2.set(title=f"MFC − analytic (max {dm:.2f})", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(im, ax=ax2, shrink=0.8)
    jmid = m // 2
    ax3 = fig.add_subplot(gs[1, :2])
    ax3.plot(yc, Tex[:, jmid], "k-", lw=2, label="analytic")
    ax3.plot(yc[::3], Tnum[::3, jmid], "C3o", ms=4, mfc="none", label="MFC")
    ax3.set(xlabel="y", ylabel="T", title="centerline probe at x = Lx/2")
    ax3.legend()
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.axis("off")
    txt = (
        f"grid {m}×{n}\nmax|u| = {umax:.1e}\nsteady drift = {drift:.1e}\n\n"
        f"full-field:\n  L1 {nrm['L1']:.3f}\n  L2 {nrm['L2']:.3f}\n  L∞ {nrm['Linf']:.3f}\n\n"
        f"interior 5–95%:\n  L1 {nrm_int['L1']:.4f}\n  L2 {nrm_int['L2']:.4f}\n  L∞ {nrm_int['Linf']:.4f}"
    )
    ax4.text(0.0, 0.95, txt, va="top", family="monospace", fontsize=10)
    fig.suptitle("2D heat: steady plate (three edges held at 10, top edge at 110)", fontweight="bold")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "heat_2d_plate.png"), dpi=130)
    plt.close(fig)
    print(f"  2D plate: full L2={nrm['L2']:.3f}  interior L2={nrm_int['L2']:.4f}  L∞(int)={nrm_int['Linf']:.4f}  drift={drift:.1e}  max|u|={umax:.1e}")
    save_summary("heat_2d_plate", {"N": int(m), "TCOLD": TCOLD, "THOT": THOT, "max_u": umax, "steady_drift": drift, "norms_full": nrm, "norms_interior": nrm_int})


# 3D test: a sine wave in a box that should cool off evenly in every direction


def bench_3d_mode():
    rundir = run_case("case_3d_mode.py", {}, 8, "3d_mode")
    dt, (m, n, p), (xc, yc, zc), steps, load = read_run(rundir)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    t = np.array([s * dt for s in steps])
    rate_an = ALPHA * 3.0 * KW**2
    L1s, L2s, Lis, amp, umax = [], [], [], [], 0.0
    field_last = ex_last = None
    for i, s in enumerate(steps):
        fields = load(s)
        T = np.transpose(fields[ndims_of(fields) + 3], (2, 1, 0))  # [z,y,x]->[x,y,z]
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
    rundir = run_case("case_3d_hotspot.py", {}, 8, "3d_hotspot")
    dt, (m, n, p), (xc, yc, zc), steps, load = read_run(rundir)
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing="ij")
    R = np.sqrt((X - L / 2) ** 2 + (Y - L / 2) ** 2 + (Z - L / 2) ** 2)
    t = np.array([s * dt for s in steps])
    T_last = np.transpose(load(steps[-1])[ndims_of(load(steps[-1])) + 3], (2, 1, 0))
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
    diffnum = 0.2
    dxs, L2s, Lis = [], [], []
    for N in grids:
        dx = L / N
        dt = diffnum * dx**2 / ALPHA
        nsteps = int(round(tstar / dt))
        rundir = run_case("case_conv.py", {"CONV_N": N, "CONV_DT": dt, "CONV_NSTEPS": nsteps}, 1, f"convx_{N}")
        dtr, (m, _, _), (xc, _, _), steps, load = read_run(rundir)
        nr = norms(temperature(load(steps[-1]), "scalar"), exact_conv(xc, steps[-1] * dtr))
        dxs.append(dx)
        L2s.append(nr["L2"])
        Lis.append(nr["Linf"])
        print(f"  N={N:4d}  dx={dx:.4e}  L2={nr['L2']:.4e}  L∞={nr['Linf']:.4e}")
    dxs, L2s, Lis = map(np.array, (dxs, L2s, Lis))
    s2 = slope(dxs, L2s)
    si = slope(dxs, Lis)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.loglog(dxs, L2s, "C0-o", label=f"L2  (slope {s2:.2f})")
    ax.loglog(dxs, Lis, "C3-^", label=f"L∞  (slope {si:.2f})")
    ax.loglog(dxs, L2s[-1] * (dxs / dxs[-1]) ** 2, "k--", lw=1.3, label="reference slope 2")
    ax.set(xlabel="cell size Δx", ylabel="error at a fixed time", title="Grid convergence study: error vs cell size (should halve twice per halving)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "convergence_spatial.png"), dpi=140)
    plt.close(fig)
    print(f"  GRID convergence slope: L2={s2:.3f}  L∞={si:.3f}  (should be ~2)")
    save_summary("convergence_spatial", {"grids": grids, "dx": dxs.tolist(), "L2": L2s.tolist(), "Linf": Lis.tolist(), "slope_L2": s2, "slope_Linf": si, "tstar": tstar})


def converge_temporal():
    # Use a coarse grid on purpose. On a fine grid the stable time step is so tiny
    # that the time-stepping error vanishes into round-off and there's nothing to
    # measure. We compare each run against a very-small-step reference run, so the
    # cell-size error cancels out and only the time-step error is left.
    N = 32
    tstar = 0.3 / (ALPHA * KW**2)
    step_counts = [50, 100, 200, 400, 800]
    ref_steps = 3200
    fields = {}
    for ns in step_counts + [ref_steps]:
        dt = tstar / ns
        rundir = run_case("case_conv.py", {"CONV_N": N, "CONV_DT": dt, "CONV_NSTEPS": ns}, 1, f"convt_{ns}")
        dtr, _, (xc, _, _), steps, load = read_run(rundir)
        fields[ns] = temperature(load(steps[-1]), "scalar")
        print(f"  steps={ns:5d}  dt={dt:.3e}")
    ref = fields[ref_steps]
    dts, L2s = [], []
    for ns in step_counts:
        e = np.abs(fields[ns] - ref)
        dts.append(tstar / ns)
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
    cmds = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not cmds or cmds[0] not in COMMANDS:
        print("usage: validate.py {" + "|".join(COMMANDS) + "} [--no-run]")
        sys.exit(1)
    COMMANDS[cmds[0]]()
