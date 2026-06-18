#!/usr/bin/env python3
"""Plot a derived field from a 2D thermocapillary run -- three views in one tool.

  temperature    EOS-recovered T field + centerline profile vs the frozen initial linear T
  sigma          sigma(T) field + sigma along the interface vs angle (the Marangoni driver)
  recirculation  drop-frame streamlines (colored by speed) + cell-resolved vorticity

Temperature is not stored by MFC; it is recovered per cell from the stiffened-gas EOS
    T = (p + p_inf) / ((gamma - 1) * rho * cv),    p from the conserved internal energy.
All constants come from simulation.inp (never hardcode a grid: a stale Nx silently mis-slices the
data into nonsense). Conserved layout (model_eqns=3, num_fluids=2): 0,1 = partial densities,
2 = x-momentum, 3 = y-momentum, 7,8 = phasic internal energies, color c last (or second-to-last when
a thermal_scalar T_s is appended).

Usage:  python3 plot_fields.py [case_dir] [temperature|sigma|recirculation] [step]
        recirculation's 3rd arg is a t/tau target instead of a step; default field = temperature.
Writes: <case_dir>/viz/<field>_<step>.png   (recirculation also writes figures/..._recirculation.png/.pdf)
"""

import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np

case_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
field = sys.argv[2].lower() if len(sys.argv) > 2 else "temperature"
if field not in ("temperature", "sigma", "recirculation"):
    sys.exit(f"unknown field {field!r}; choose temperature, sigma, or recirculation")
HERE = os.path.dirname(os.path.abspath(__file__))


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


P = read_namelist(os.path.join(case_dir, "simulation.inp"))
param = lambda k: float(P[k.lower()])  # noqa: E731

nx, ny = int(param("m")) + 1, int(param("n")) + 1
dt, cv = param("dt"), param("fluid_pp(1)%cv")
# MFC stores gamma_mfc = 1/(gamma-1) and pi_inf_mfc = gamma*p_inf/(gamma-1); invert to physical.
gamma_mfc, pi_inf_mfc = param("fluid_pp(1)%gamma"), param("fluid_pp(1)%pi_inf")
gamma = 1.0 + 1.0 / gamma_mfc
p_inf = pi_inf_mfc * (gamma - 1.0) / gamma
ts_mode = str(P.get("thermal_scalar", "F")).strip(". ").upper().startswith("T")

restart_dir = os.path.join(case_dir, "restart_data")
xb = np.fromfile(os.path.join(restart_dir, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
x = 0.5 * (xb[:-1] + xb[1:])
y = 0.5 * (yb[:-1] + yb[1:])
ncell = nx * ny

steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
if not steps:
    sys.exit(f"no lustre_<step>.dat snapshots in {restart_dir!r} -- run the case first")
nvars = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size // ncell
c_idx = nvars - 2 if ts_mode else nvars - 1  # color c (T_s appended after it in ts mode)
viz_dir = os.path.join(case_dir, "viz")
os.makedirs(viz_dir, exist_ok=True)


def columns_of(step):
    """A snapshot as (nvars, ncell): row i is conserved variable i flattened."""
    return np.fromfile(os.path.join(restart_dir, f"lustre_{step}.dat"), np.float64).reshape(-1, ncell)


def eos_temperature(cols):
    """(T, color) on the (ny, nx) grid from a (nvars, ncell) snapshot, via the stiffened-gas EOS."""
    rho = (cols[0] + cols[1]).reshape(ny, nx)
    rho_e = (cols[7] + cols[8]).reshape(ny, nx)  # phasic internal energies (no kinetic part)
    p = (rho_e - pi_inf_mfc) / gamma_mfc
    T = (p + p_inf) / ((gamma - 1.0) * rho * cv)
    color = np.clip(cols[c_idx].reshape(ny, nx), 0.0, None)
    return T, color


def plot_temperature(step):
    """EOS temperature field + centerline profile vs the frozen initial linear T."""
    T0, gradT = 10.0, 1.0  # IC T(x)=T0+gradT*x lives in case_Ma_0.py's density string, set here for the reference line
    T, c = eos_temperature(columns_of(step))
    T_initial, _ = eos_temperature(columns_of(steps[0]))
    if not np.all(np.isfinite(T)) or T.min() <= 0.0:
        sys.exit(f"reconstructed T is non-physical (min = {T.min():.4f}); constants likely mismatch the data.")
    print(f"step {step}  (t = {step * dt:.4f})   T range: [{T.min():.4f}, {T.max():.4f}]")
    fig, (ax_field, ax_line) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0]})
    vmin = min(T0 + gradT * x[0], T.min())
    vmax = max(T0 + gradT * x[-1], T.max())
    mesh = ax_field.pcolormesh(x, y, T, cmap="coolwarm", vmin=vmin, vmax=vmax, shading="auto")
    fig.colorbar(mesh, ax=ax_field, label="temperature $T$")
    ax_field.contour(x, y, c, levels=[0.5], colors="k", linewidths=1.2)
    xc = (c * x[None, :]).sum() / c.sum()
    yc = (c * y[:, None]).sum() / c.sum()
    ax_field.plot(xc, yc, "k+", ms=10, mew=2)
    ax_field.set(aspect="equal", xlabel="$x$", ylabel="$y$", title=f"Temperature field (step {step}, $t={step * dt:.3f}$)\ndroplet centroid $x={xc:+.4f}$")
    mid = ny // 2
    ax_line.plot(x, T0 + gradT * x, "--", color="0.5", label=r"initial $T(x)=T_0+\nabla T\,x$ (frozen IC)")
    ax_line.plot(x, T_initial[mid], ":", color="C0", alpha=0.7, label=f"centerline, step {steps[0]}")
    ax_line.plot(x, T[mid], "-", color="C3", label=f"centerline, step {step}")
    ax_line.set(xlabel="$x$", ylabel=r"temperature $T$ (at $y \approx 0$)", title="Centerline temperature profile")
    ax_line.legend(loc="upper left", fontsize=8)
    ax_line.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(viz_dir, f"temperature_{step}.png")
    fig.savefig(out, dpi=150)
    print(f"saved figure -> {out}")


def plot_sigma(step):
    """sigma(T) field + sigma in the interface cells vs angle (the Marangoni driver)."""
    sigma0, dsigma_dT, T_ref = param("sigma"), param("sigma_dtdt"), param("sigma_t_ref")
    T, c = eos_temperature(columns_of(step))
    sigma = sigma0 + dsigma_dT * (T - T_ref)
    xc = (c * x[None, :]).sum() / c.sum()
    yc = (c * y[:, None]).sum() / c.sum()
    X, Y = np.meshgrid(x, y)
    interface = (c > 0.2) & (c < 0.8)
    angle = np.degrees(np.arctan2(Y[interface] - yc, X[interface] - xc))
    sig_if = sigma[interface]
    fig, (ax_field, ax_line) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.05, 1]})
    mesh = ax_field.pcolormesh(x, y, sigma, cmap="viridis", shading="auto")
    fig.colorbar(mesh, ax=ax_field, label=r"$\sigma(T)$")
    ax_field.contour(x, y, c, levels=[0.5], colors="w", linewidths=1.2)
    ax_field.set(aspect="equal", xlabel="x", ylabel="y", title=rf"$\sigma(T)$ field + interface (step {step}, t = {step * dt:.2f})")
    ax_line.scatter(angle, sig_if, s=12, color="C0")
    ax_line.axvline(0, ls=":", color="C3")
    ax_line.text(6, sig_if.min(), "hot (+x)", color="C3", fontsize=9)
    ax_line.set(xlabel="angle around interface (deg)", ylabel=r"$\sigma$ in the interface cells", title=r"Low $\sigma$ on the hot side drives the Marangoni pull")
    ax_line.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(viz_dir, f"sigma_interface_{step}.png")
    fig.savefig(out, dpi=150)
    print(f"sigma on interface spans [{sig_if.min():.4f}, {sig_if.max():.4f}]")
    print(f"saved -> {out}")


def plot_recirculation(target_ttau):
    """Drop-frame streamlines (colored by speed) + cell-resolved vorticity, near the drop."""
    mu = 1.0 / param("fluid_pp(1)%re(1)")
    dsigma_dT = param("sigma_dtdt")
    r, gradT = 0.5, 2.0 / 15.0
    v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu

    def snap_uvc(step):
        cols = columns_of(step)
        rho = (cols[0] + cols[1]).reshape(ny, nx)
        return rho, (cols[2].reshape(ny, nx)) / rho, (cols[3].reshape(ny, nx)) / rho, np.clip(cols[c_idx].reshape(ny, nx), 0.0, 1.0)

    rho0, _, _, c0 = snap_uvc(steps[0])
    rho_drop = (c0 * rho0).sum() / c0.sum()
    tau = rho_drop * r**2 / mu  # viscous time from the t=0 drop density (tracks the case, not rho=1)
    step = min(steps, key=lambda s: abs(s * dt / tau - target_ttau))
    ttau = step * dt / tau
    _, u, v, c = snap_uvc(step)
    u_drop = (c * u).sum() / c.sum()
    v_drop = (c * v).sum() / c.sum()
    uc, vc = u - u_drop, v - v_drop  # co-moving frame
    omega = np.gradient(vc, x, axis=1) - np.gradient(uc, y, axis=0)
    xc_drop = (c * x[None, :]).sum() / c.sum()
    yc_drop = (c * y[:, None]).sum() / c.sum()
    win = 2.5 * r
    mx = (x > xc_drop - win) & (x < xc_drop + win)
    my = (y > yc_drop - win) & (y < yc_drop + win)
    xz, yz = x[mx], y[my]
    uz, vz = uc[np.ix_(my, mx)], vc[np.ix_(my, mx)]
    cz, oz = c[np.ix_(my, mx)], omega[np.ix_(my, mx)]
    speed = np.hypot(uz, vz)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "font.size": 9,
            "axes.titlesize": 9.5,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "axes.linewidth": 0.8,
        }
    )
    coldc, hotc, inkc = (46 / 255, 86 / 255, 149 / 255), (171 / 255, 57 / 255, 52 / 255), (0.12, 0.12, 0.12)
    halo = [pe.withStroke(linewidth=2.4, foreground="white")]
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.2, 3.9), sharey=True, constrained_layout=True)
    strm = axa.streamplot(xz, yz, uz, vz, color=speed / v_YGB, cmap="viridis", density=1.3, linewidth=0.75, arrowsize=0.6)
    cb = fig.colorbar(strm.lines, ax=axa, fraction=0.046, pad=0.03)
    cb.set_label(r"$|\mathbf{u}-\mathbf{U}_{\rm drop}|\,/\,v_{\rm YGB}$")
    cb.outline.set_linewidth(0.6)
    axa.contour(xz, yz, cz, levels=[0.5], colors="white", linewidths=2.2)
    axa.contour(xz, yz, cz, levels=[0.5], colors=[inkc], linewidths=0.9)
    axa.set_title("(a) drop-frame streamlines")
    o_nd = oz * r / v_YGB
    vmax = np.percentile(np.abs(o_nd), 99)
    pm = axb.pcolormesh(xb[mx.nonzero()[0][0] : mx.nonzero()[0][-1] + 2], yb[my.nonzero()[0][0] : my.nonzero()[0][-1] + 2], o_nd, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat", rasterized=True)
    cb2 = fig.colorbar(pm, ax=axb, fraction=0.046, pad=0.03)
    cb2.set_label(r"$\omega_z\, r / v_{\rm YGB}$")
    cb2.outline.set_linewidth(0.6)
    axb.contour(xz, yz, cz, levels=[0.5], colors="white", linewidths=2.0)
    axb.contour(xz, yz, cz, levels=[0.5], colors=[inkc], linewidths=0.9)
    axb.set_title(r"(b) vorticity field")
    for ax in (axa, axb):
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x/D$")
        for col, lw in (("white", 4.0), (inkc, 1.8)):
            ax.annotate("", xy=(xc_drop, yc_drop + 1.0 * r), xytext=(xc_drop, yc_drop - 0.05 * r), arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, shrinkA=0, shrinkB=0))
        ax.text(0.04, 0.04, "cold", transform=ax.transAxes, color=coldc, fontsize=8, ha="left", va="bottom", path_effects=halo)
        ax.text(0.04, 0.96, "hot", transform=ax.transAxes, color=hotc, fontsize=8, ha="left", va="top", path_effects=halo)
    axa.set_ylabel(r"$y/D$")
    axa.text(xc_drop + 0.28 * r, yc_drop + 0.15 * r, rf"$U={v_drop / v_YGB:.2f}\,v_{{\rm YGB}}$", color=inkc, fontsize=9, ha="left", va="bottom", path_effects=halo)
    fig.suptitle(rf"2D thermocapillary drop  $\cdot$  {ny / 7.5:.1f} cells/$D$  $\cdot$  $t={ttau:.1f}\,\tau$", fontsize=9.5)
    out = os.path.join(HERE, "figures", "case1_zero_marangoni_2D_recirculation")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out + ".png", dpi=300)
    fig.savefig(out + ".pdf")
    print(f"step={step}  t/tau={ttau:.3f}  nvars={nvars}  rho_drop={rho_drop:.3f}  tau={tau:.3f}")
    print(f"U_drop/v_YGB={v_drop / v_YGB:+.3f}  (lateral u_drop/v_YGB={u_drop / v_YGB:+.3f})")
    print(f"saved -> {out}.png / .pdf")


if field == "recirculation":
    plot_recirculation(float(sys.argv[3]) if len(sys.argv) > 3 else 2.6)
else:
    step = int(sys.argv[3]) if len(sys.argv) > 3 else steps[-1]
    if step not in steps:
        sys.exit(f"step {step} not available; choose from {steps[0]}..{steps[-1]}")
    (plot_temperature if field == "temperature" else plot_sigma)(step)
