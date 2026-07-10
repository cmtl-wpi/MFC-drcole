#!/usr/bin/env python3
"""
Measure the surfactant surface-diffusion mode-decay rate(s) and compare to the exact D_s * k^2.

Reads grid/dt/D_s from the run's simulation.inp (never hard-coded) and the conserved surfactant
density (the last conserved variable) from the parallel-IO restart files.

    python3 measure.py           # single-mode decay -> figures/decay.png       (run case.py first)
    python3 measure.py sweep     # dispersion across wavenumbers -> figures/dispersion.png (run sweep.py first)

The dispersion run (sweep.py) seeds a superposition of modes; because surface diffusion is linear
each Fourier component decays independently at its own eigenvalue rate, so one run yields the whole
rate(k) = D_s k^2 spectrum -- the canonical Laplace-Beltrami mode-decay benchmark (e.g. Xu, Li,
Lowengrub & Zhao, JCP 2006; the surface-FEM literature, Dziuk & Elliott).
"""

import glob
import math
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def inp(name, path):
    m = re.search(rf"^\s*{re.escape(name)}\s*=\s*([^\s,]+)", open(path).read(), re.M)
    if not m:
        raise KeyError(name)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


sim = os.path.join(HERE, "simulation.inp")
nx = int(inp("m", sim)) + 1
ny = int(inp("n", sim)) + 1
dt = inp("dt", sim)
D_s = inp("surf_diff", sim)
x0 = inp("x_domain%beg", sim)
Lx = inp("x_domain%end", sim) - x0
dx = Lx / nx
x_cc = x0 + (np.arange(nx) + 0.5) * dx

files = [f for f in glob.glob(os.path.join(HERE, "restart_data", "lustre_*.dat")) if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])]
files.sort(key=lambda f: int(os.path.basename(f)[7:-4]))
ss = os.path.getsize(files[0]) // 8 // (nx * ny)  # surfactant is the last conserved variable
ts = np.array([int(os.path.basename(f)[7:-4]) * dt for f in files])
surf = np.array([np.fromfile(f, "<f8").reshape(ss, ny * nx)[ss - 1] for f in files])
M0 = surf.sum(axis=1)


def mode_rate(n):
    """Fit the decay rate of azimuthal mode n (k = n*2*pi/Lx) over its clean exponential window."""
    k = n * 2.0 * math.pi / Lx
    C, S = np.tile(np.cos(k * x_cc), ny), np.tile(np.sin(k * x_cc), ny)
    A = np.hypot(surf @ C, surf @ S)
    r = A / A[0]
    w = (r > 0.12) & (r < 0.95)  # avoid the t=0 point and the noise floor
    rate = -np.polyfit(ts[w], np.log(r[w]), 1)[0]
    return k, rate, D_s * k**2, r


if "sweep" in sys.argv:
    modes = [1, 2, 3]
    print(f"total surfactant drift : {100*(M0[-1]/M0[0]-1):+.3f}%")
    ks, meas, exact = [], [], []
    for n in modes:
        k, rate, ana, _ = mode_rate(n)
        ks.append(k)
        meas.append(rate)
        exact.append(ana)
        print(f"  mode n={n}  k={k:6.3f}  measured={rate:8.4f}  exact D_s k^2={ana:8.4f}  err={100*(rate/ana-1):+.2f}%")
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    kk = np.linspace(0, max(ks) * 1.05, 100)
    ax.plot(kk, D_s * kk**2, "-", color="#333", lw=2, label=r"exact  $D_s k^2$")
    ax.plot(ks, meas, "o", ms=9, mfc="#d1495b", mec="k", mew=0.6, label="MFC")
    for k, m in zip(ks, meas):
        ax.annotate(rf"$n={round(k*Lx/(2*math.pi))}$", (k, m), textcoords="offset points", xytext=(8, -10))
    ax.set_xlabel(r"wavenumber  $k$")
    ax.set_ylabel(r"decay rate")
    ax.set_title(r"Surface-diffusion dispersion: rate $=D_s k^2$ across modes")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "dispersion.png"), dpi=130)
    print("wrote figures/dispersion.png")
else:
    k, rate, ana, r = mode_rate(1)
    print(f"total surfactant drift : {100*(M0[-1]/M0[0]-1):+.3f}%")
    print(f"measured decay rate    : {rate:.4f}")
    print(f"exact  D_s k^2         : {ana:.4f}")
    print(f"relative error         : {100*(rate/ana-1):+.2f}%")
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    tt = np.linspace(0, ts[-1], 200)
    ax.semilogy(tt, np.exp(-ana * tt), "-", color="#333", lw=2, label=rf"exact  $e^{{-D_s k^2 t}}$   ($D_s k^2={ana:.3f}$)")
    ax.semilogy(ts, r, "o", ms=6, mfc="#d1495b", mec="k", mew=0.5, label=rf"MFC  (fit ${rate:.3f}$, err ${100*(rate/ana-1):+.1f}\%$)")
    ax.set_xlabel("time  $t$")
    ax.set_ylabel(r"mode amplitude  $A(t)/A(0)$")
    ax.set_title("Insoluble-surfactant surface diffusion: flat-interface mode decay")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "decay.png"), dpi=130)
    print("wrote figures/decay.png")
