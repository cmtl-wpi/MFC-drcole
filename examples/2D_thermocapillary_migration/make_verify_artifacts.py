#!/usr/bin/env python3
# Turn a finished 1D conduction-verify run into durable validation artifacts: a PNG showing the
# raw T(x) field alongside the extracted scalar, and an entry in results/verify/summary.json.
# One generator for all three harnesses so the measurement matches what each verify_*.py checks:
#   diffusion  -> verify_1d_diffusion.py      (sinusoidal T imposed via density, decays at alpha*kappa^2)
#   scalar     -> verify_1d_thermal_scalar.py (sinusoidal T_s imposed on the scalar, density uniform)
#   conduction -> verify_1d_conduction.py     (linear T => div(k grad T)=0 => static equilibrium, u~0)
#
#   python3 make_verify_artifacts.py {diffusion|scalar|conduction} <rundir>
import glob
import json
import math
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, "results", "verify")

# Physical constants -- mirror the matching verify_*.py (kept inline so this reads standalone).
GAM, P_INF, CV, P0, T0 = 2.0, 100.0, 12.5, 25.0, 10.0


def read_run(rundir):
    """Return (dt, cells, x_cc, steps, load) for an MFC 1D restart run."""
    rd = os.path.join(rundir, "restart_data")
    inp = {}
    for line in open(os.path.join(rundir, "simulation.inp")):
        if "=" in line:
            k, v = line.split("=", 1)
            inp[k.strip().lower()] = v.strip().rstrip(",")
    dt = float(inp["dt"])
    cells = int(inp["m"]) + 1
    xbeg, xend = float(inp["x_domain%beg"]), float(inp["x_domain%end"])
    x = xbeg + (np.arange(cells) + 0.5) * (xend - xbeg) / cells
    steps = sorted(int(g.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (g := re.search(r"lustre_(\d+)\.dat$", os.path.basename(f))))

    def load(s):
        a = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        return a.reshape(a.size // cells, cells)

    return dt, cells, x, steps, load


def temperature_from_eos(d):
    """Stiffened-gas temperature from a conserved-variable record [rho, mom, E, ...]."""
    rho, mom, E = d[0], d[1], d[2]
    u = mom / rho
    e = E / rho - 0.5 * u**2
    pres = (GAM - 1.0) * rho * e - GAM * P_INF
    T = (pres + P_INF) / ((GAM - 1.0) * rho * CV)
    return T, u


def fit_decay(t, amp):
    """Decay rate from ln|A/A0| vs t over the resolvable window."""
    a0 = amp[0]
    good = np.abs(amp) > 1e-12
    rate = -np.polyfit(t[good], np.log(np.abs(amp[good] / a0)), 1)[0]
    return a0, rate


def save_summary(key, entry):
    os.makedirs(ART, exist_ok=True)
    path = os.path.join(ART, "summary.json")
    summary = json.load(open(path)) if os.path.exists(path) else {}
    summary[key] = entry
    json.dump(summary, open(path, "w"), indent=2)
    print(f"  wrote {path} [{key}]")


def artifact_decay(rundir, kind):
    """diffusion (T via EOS) or scalar (T_s read directly): compare the cell-resolved field at every
    output directly against the CLOSED-FORM conduction solution
        T(x,t) = T0 * (1 + eps*sin(kappa x) * exp(-alpha*kappa^2*t)).
    Left panel overlays MFC vs that solution; right panel is the pointwise deviation over the run.
    The mode-projected decay rate is also fit and kept in the summary as a single-number check."""
    L = 1.0
    eps = 0.03 if kind == "diffusion" else 0.05
    alpha = 0.05
    kappa = 2.0 * math.pi / L
    rate_an = alpha * kappa**2  # analytic decay rate of the fundamental mode
    dt, cells, x, steps, load = read_run(rundir)
    basis = np.sin(kappa * x)
    norm = (basis**2).sum()
    t, amp, umax, fields = [], [], [], []
    for s in steps:
        d = load(s)
        if kind == "scalar":
            field, u = d[4], d[1] / d[0]  # T_s carried as conserved var 5 (rho,mom,E,alpha,T_s)
        else:
            field, u = temperature_from_eos(d)
        fields.append(field)
        amp.append((field * basis).sum() / norm)
        umax.append(np.max(np.abs(u)))
        t.append(s * dt)
    t, amp, fields = np.array(t), np.array(amp), np.array(fields)
    _, rate = fit_decay(t, amp)
    rate_err = 100 * abs(rate / rate_an - 1)
    max_u = float(max(umax))

    # closed-form analytic field at each output time, and the pointwise deviation from it
    analytic = T0 * (1.0 + eps * np.exp(-rate_an * t)[:, None] * basis[None, :])
    abserr = np.abs(fields - analytic)
    linf = abserr.max(axis=1) / (eps * T0) * 100.0  # L-inf error, % of the initial amplitude eps*T0
    rms = np.sqrt((abserr**2).mean(axis=1)) / (eps * T0) * 100.0
    max_linf = float(linf.max())

    fig, (axF, axE) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    label = "T(x)" if kind == "diffusion" else "T_s(x)"
    idx = np.unique(np.linspace(0, len(t) - 1, 4).round().astype(int))  # ~4 discrete times across the run
    colors = [plt.cm.viridis(0.12 + 0.76 * k / max(1, len(idx) - 1)) for k in range(len(idx))]
    for c, j in zip(colors, idx):
        axF.plot(x, fields[j], "-", color=c, lw=2.0, zorder=2)  # MFC: solid line
        axF.plot(x, analytic[j], "--", color=c, lw=1.3, zorder=3)  # analytic closed form: dashed, drawn on top
    # two legends: line style says which solution, color says which time
    style_leg = axF.legend(
        [Line2D([], [], color="0.3", lw=2.0, ls="-"), Line2D([], [], color="0.3", lw=1.3, ls="--")],
        ["MFC (solid)", "analytic, T0(1+ε·sin·e$^{-ακ²t}$) (dashed)"],
        loc="lower left", fontsize=8, framealpha=0.9,
    )
    axF.add_artist(style_leg)
    axF.legend([Line2D([], [], color=c, lw=2.0) for c in colors], [f"t={t[j]:.2f}" for j in idx], loc="upper right", fontsize=8, framealpha=0.9, title="time", ncol=2)
    axF.set(xlabel="x", ylabel=label, title="cell-resolved field: MFC (solid) vs closed form (dashed)")
    axE.plot(t, linf, "C3-", lw=1.6, label="max|T$_{sim}$ − T$_{analytic}$|")
    axE.plot(t, rms, "C0-", lw=1.6, label="RMS")
    axE.set(xlabel="t", ylabel="error  (% of initial amplitude εT0)", title=f"pointwise deviation from analytic (peak {max_linf:.2f}%)")
    axE.legend(fontsize=8)
    name = "diffusion" if kind == "diffusion" else "thermal_scalar"
    fig.suptitle(f"verify_1d_{name}: field vs analytic conduction solution   (mode decay {rate:.3f} vs {rate_an:.3f}, {rate_err:.1f}%)", fontweight="bold")
    fig.tight_layout()
    png = os.path.join(ART, f"{kind}.png")
    fig.savefig(png, dpi=130)
    plt.close(fig)

    print(f"  decay rate {rate:.4f} vs {rate_an:.4f} ({rate_err:.1f}%)   field error: peak {max_linf:.2f}%  final {linf[-1]:.2f}%   max|u|={max_u:.2e}")
    save_summary(kind, {"measured_rate": rate, "analytic_rate": rate_an, "rate_error_pct": rate_err, "field_err_peak_pct": max_linf, "field_err_final_pct": float(linf[-1]), "max_u": max_u, "n_outputs": len(steps), "png": os.path.basename(png)})


def artifact_conduction(rundir, label=""):
    """Linear T(x) => div(k grad T)=0 => exact static equilibrium: u must stay ~0, T stays linear.

    label distinguishes runs of the same case (e.g. "_mpi2" for the 2-rank run that exercises the
    isothermal-BC halo guard -- the code path behind the 2D reversal, invisible to a single rank).
    """
    Lx, gradT, r, mu, dsig = 7.5, 2.0 / 15.0, 0.5, 0.1, -0.1
    U_r = (-dsig) * gradT * r / mu  # interfacial velocity scale, for context on how small u is
    dt, cells, x, steps, load = read_run(rundir)
    Tlin = T0 + gradT * x
    t, umax, T0field, Tend = [], [], None, None
    for s in steps:
        T, u = temperature_from_eos(load(s))
        umax.append(np.max(np.abs(u)))
        t.append(s * dt)
        if T0field is None:
            T0field = T
        Tend = T
    t, umax = np.array(t), np.array(umax)
    max_u = float(umax.max())
    dev = np.abs(Tend - Tlin)
    T_max_dev, T_rms_dev = float(dev.max()), float(np.sqrt((dev**2).mean()))

    fig, (axT, axU) = plt.subplots(1, 2, figsize=(11, 4.2))
    axT.plot(x, Tlin, "k-", lw=2.5, alpha=0.4, label="analytic  T0+gradT·x")
    axT.plot(x, T0field, "C0.", ms=4, label="initial")
    axT.plot(x, Tend, "C3.", ms=4, label="final")
    axT.set(xlabel="x", ylabel="T", title=f"linear T preserved (max dev {T_max_dev:.2e})")
    axT.legend()
    axU.plot(t, umax / U_r, "C2-")
    axU.set(xlabel="t", ylabel="max|u| / U_r", title=f"spurious velocity stays ~0  (max {max_u / U_r:.2e}·U_r)")
    axU.axhline(0, color="k", lw=0.6)
    ranks = " (2 ranks: internal MPI boundary)" if label else ""
    fig.suptitle(f"verify_1d_conduction{label}: linear-T static equilibrium (flux + Dirichlet BC){ranks}", fontweight="bold")
    fig.tight_layout()
    png = os.path.join(ART, f"conduction{label}.png")
    fig.savefig(png, dpi=130)
    plt.close(fig)

    print(f"  max|u|={max_u:.3e}  (={max_u / U_r:.2e}·U_r)   T_max_dev={T_max_dev:.3e}  T_rms_dev={T_rms_dev:.3e}")
    save_summary(f"conduction{label}", {"max_u": max_u, "max_u_over_Ur": max_u / U_r, "U_r": U_r, "T_max_dev": T_max_dev, "T_rms_dev": T_rms_dev, "n_outputs": len(steps), "png": os.path.basename(png)})


if __name__ == "__main__":
    test, rundir = sys.argv[1], sys.argv[2]
    label = sys.argv[3] if len(sys.argv) > 3 else ""
    if test == "conduction":
        artifact_conduction(rundir, label)
    else:
        artifact_decay(rundir, test)
