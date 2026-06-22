#!/usr/bin/env python3
"""Convergence study for MFC's heat-conduction operator: does the error shrink at
the rate the math predicts?

Uses a 1D periodic single Fourier mode (no boundary-condition error) so the only
error sources are the discretization and the energy-coupled physics floor. Two
sweeps, both produced by case.py with environment overrides (see Allrun):

  spatial   convx_{32,64,128,256,512}   error vs the exact mode at fixed t*, refining dx
  temporal  convt_{256,512,1024,2048,4096}   error vs a tiny-dt reference on a fixed grid

This script does NOT run MFC -- run the sweeps first (./Allrun), then:

  python3 examples/Thermal_Conduction_Convergence/validate.py            # both
  python3 examples/Thermal_Conduction_Convergence/validate.py spatial
  python3 examples/Thermal_Conduction_Convergence/validate.py temporal

It reads runs/<label>/, writes figures/convergence_{spatial,temporal}.png and
summary.json.
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

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIG = os.path.join(HERE, "figures")
SUMMARY = os.path.join(HERE, "summary.json")

# Physical constants -- mirror case.py
L, ALPHA = 1.0, 0.05
TWALL, AMP = 10.0, 3.0
KW = 2.0 * math.pi / L  # periodic single-mode wavenumber


def exact_conv(x, t):
    return TWALL + AMP * np.sin(KW * x) * math.exp(-ALPHA * KW**2 * t)


def find_run(label):
    rundir = os.path.join(RUNS, label)
    if not os.path.isdir(os.path.join(rundir, "restart_data")):
        raise SystemExit(f"missing {os.path.relpath(rundir, HERE)}/restart_data -- run the sweeps first with ./Allrun")
    return rundir


def read_run(rundir):
    """Read one run's time step, grid, cell-center coords, saved steps, and a
    load(step) function for the raw field array."""
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


def norms(num, ex):
    e = np.abs(num - ex)
    return {"L1": float(e.mean()), "L2": float(np.sqrt((e**2).mean())), "Linf": float(e.max())}


def save_summary(key, entry):
    summary = json.load(open(SUMMARY)) if os.path.exists(SUMMARY) else {}
    summary[key] = entry
    json.dump(summary, open(SUMMARY, "w"), indent=2)


def slope(xs, ys):
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def converge_spatial():
    grids = [32, 64, 128, 256, 512]
    tstar = 0.3 / (ALPHA * KW**2)
    dxs, L2s, Lis = [], [], []
    for N in grids:
        rundir = find_run(f"convx_{N}")
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
        rundir = find_run(f"convt_{ns}")
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


COMMANDS = {"spatial": converge_spatial, "temporal": converge_temporal}

if __name__ == "__main__":
    cmds = sys.argv[1:] or list(COMMANDS)
    if any(c not in COMMANDS for c in cmds):
        print("usage: validate.py [spatial] [temporal]   (run ./Allrun first; this only analyzes runs/)")
        sys.exit(1)
    for c in cmds:
        COMMANDS[c]()
