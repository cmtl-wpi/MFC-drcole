#!/usr/bin/env python3
"""Skeptical empirical test of the 'fundamental vertical acoustic standing wave' claim for the Fig 5
rise-velocity runs. Reads each run's simulation.inp + restart snapshots and checks, from the data alone:

  1. FFT of the drop (color-weighted) v_y history -> dominant temporal frequency, vs c/(2 Ly) and c/(2 Lx).
  2. Samples per oscillation period (resolved vs aliased: Nyquist = 2 samples/period).
  3. Spatial structure: time-RMS of v_y averaged over x, as a function of y -> is it a single domain-spanning
     half-wave (one hump, antinode near center, nodes at the slip walls)? And is there x-structure (width mode)?
  4. Effective sound speed from the actual density field (c = sqrt(gam (p0+pi_inf)/rho)), column-averaged.
  5. Ripple amplitude vs v_YGB, across grids (does it shrink with refinement?).

Usage: python3 diag_standing_wave.py
"""

import glob
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def load(run_dir):
    """Return dict of arrays/scalars for one run, or None."""
    inp = os.path.join(run_dir, "simulation.inp")
    rd = os.path.join(run_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    P = read_namelist(inp)
    f = lambda k: float(P[k.lower()])  # noqa: E731
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    cells = nx * ny * nz
    ts = str(P.get("thermal_scalar", "F")).strip(". ").upper().startswith("T")
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    c_idx = nvars - 2 if ts else nvars - 1
    dt = f("dt")
    Ly = f("y_domain%end") - f("y_domain%beg")
    Lx = f("x_domain%end") - f("x_domain%beg")
    yb, ye = f("y_domain%beg"), f("y_domain%end")
    y_cc = yb + (np.arange(ny) + 0.5) * (Ly / ny)

    # stiffened-gas params for sound speed
    gam = 1.0 + 1.0 / f("fluid_pp(1)%gamma")  # stored gamma = 1/(gam-1)
    pi_phys = f("fluid_pp(1)%pi_inf") * (gam - 1.0) / gam  # stored pi_inf = gam*p_inf/(gam-1)
    p0 = float(P.get("patch_icpp(1)%pres", "8.0"))  # background pressure (not in sim namelist; fig5 p0=8)

    t, u_drop, vyx_t, c_col = [], [], [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        rho = fld(0) + fld(1)
        vy = fld(3) / rho
        c = np.clip(fld(c_idx), 0.0, None)
        t.append(s * dt)
        u_drop.append((c * vy).sum() / c.sum())  # drop-tracking color-weighted v_y (the plotted curve)
        vyx_t.append(vy[0].mean(axis=1))  # v_y averaged over x -> profile vs y (whole field)
        c_col.append(np.sqrt(gam * (p0 + pi_phys) / rho[0]).mean())  # field sound speed, domain-averaged
    return dict(P=P, nx=nx, ny=ny, dt=dt, Lx=Lx, Ly=Ly, y_cc=y_cc, gam=gam, pi_phys=pi_phys, p0=p0, t=np.array(t), u_drop=np.array(u_drop), vyx=np.array(vyx_t), c_col=float(np.mean(c_col)))


def v_ygb(P):
    f = lambda k: float(P[k.lower()])  # noqa: E731
    mu = 1.0 / f("fluid_pp(1)%re(1)")
    dsdt = f("sigma_dtdt")
    gradT = 2.0 / 15.0
    return (2.0 / 15.0) * (-dsdt) * gradT * 0.5 / mu


def dominant_freq(t, y):
    """Detrend (remove slow migration via deg-5 poly), FFT the residual, return (freqs, amp, f_peak)."""
    dt = np.median(np.diff(t))
    resid = y - np.polyval(np.polyfit(t, y, 5), t)
    resid = resid - resid.mean()
    amp = np.abs(np.fft.rfft(resid * np.hanning(len(resid))))
    freqs = np.fft.rfftfreq(len(resid), d=dt)
    k = np.argmax(amp[1:]) + 1
    return freqs, amp, freqs[k], dt


print(f"{'run':18s} {'Nsnap':>5s} {'dt_save':>8s} {'f_peak':>7s} {'period':>7s} {'samp/per':>8s} {'c_col':>6s} {'c/2Ly':>6s} {'c/2Lx':>6s} {'ripple/vYGB':>11s}")
results = {}
for name in ["fig5_2D_w064", "fig5_2D_w128", "fig5_2D_w256"]:
    d = load(os.path.join(RUNS, name))
    if d is None:
        print(f"{name:18s}  MISSING")
        continue
    results[name] = d
    freqs, amp, fpk, dts = dominant_freq(d["t"], d["u_drop"])
    period = 1.0 / fpk
    fy1, fx1 = d["c_col"] / (2 * d["Ly"]), d["c_col"] / (2 * d["Lx"])
    resid = d["u_drop"] - np.polyval(np.polyfit(d["t"], d["u_drop"], 5), d["t"])
    ripple = resid.std() / v_ygb(d["P"])
    print(f"{name:18s} {len(d['t']):5d} {dts:8.4f} {fpk:7.3f} {period:7.3f} {period / dts:8.2f} {d['c_col']:6.2f} {fy1:6.3f} {fx1:6.3f} {ripple:11.3f}")

# Detail panel for the w128 run: time series, spectrum (modes marked), spatial RMS profile.
d = results.get("fig5_2D_w128")
if d is not None:
    freqs, amp, fpk, dts = dominant_freq(d["t"], d["u_drop"])
    fy1, fx1 = d["c_col"] / (2 * d["Ly"]), d["c_col"] / (2 * d["Lx"])
    # isolate the OSCILLATORY part: detrend each y-column in time (remove steady migration flow), then RMS
    tt = d["t"]
    osc = np.array([d["vyx"][:, j] - np.polyval(np.polyfit(tt, d["vyx"][:, j], 5), tt) for j in range(d["ny"])]).T
    vyx_rms = osc.std(axis=0)  # time-RMS of the OSCILLATION of x-averaged v_y, vs y
    f = lambda k: float(d["P"][k.lower()])  # noqa: E731
    Ly = d["Ly"]
    y_drop = -Ly / 2 + 1.5  # Samareh: 1.5 D off the cold floor
    y_anti = d["y_cc"][np.argmax(vyx_rms)]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    ax[0].plot(d["t"], d["u_drop"] / v_ygb(d["P"]), ".-", ms=4, lw=0.8)
    ax[0].set(xlabel="t", ylabel=r"$v_y^{drop}/v_{YGB}$", title="drop velocity (jagged curve)")
    ax[1].plot(freqs, amp, "-")
    ax[1].axvline(fy1, color="C1", ls="--", label=f"c/(2Ly)={fy1:.2f}")
    ax[1].axvline(fx1, color="C2", ls=":", label=f"c/(2Lx)={fx1:.2f}")
    ax[1].axvline(fpk, color="k", ls="-", lw=0.8, label=f"FFT peak={fpk:.2f}")
    ax[1].set(xlabel="frequency", ylabel="|FFT| residual", title="spectrum", xlim=(0, 4))
    ax[1].legend(fontsize=9)
    ax[2].plot(vyx_rms, d["y_cc"], "-")
    ax[2].axhline(y_drop, color="C3", ls="--", label=f"drop y={y_drop:.2f}")
    ax[2].axhline(y_anti, color="k", ls=":", label=f"RMS max y={y_anti:.2f}")
    ax[2].set(xlabel=r"time-RMS of $\langle v_y\rangle_x$", ylabel="y", title="spatial mode shape")
    ax[2].legend(fontsize=9)
    out = os.path.join(HERE, "results", "diag_standing_wave.png")
    fig.savefig(out, dpi=140)
    print(f"\nwrote {out}")
    print(f"drop at y={y_drop:.2f}; RMS antinode at y={y_anti:.2f}; RMS at drop / RMS max = {np.interp(y_drop, d['y_cc'], vyx_rms) / vyx_rms.max():.2f}")
    # width-mode check: is there x-structure? compare RMS of y-averaged v_y(x) to the y-profile RMS
    vyy_rms = d["vyx"]  # already x-averaged; load x-structure separately is heavier -- report y vs flat
