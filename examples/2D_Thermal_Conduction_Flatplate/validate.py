#!/usr/bin/env python3
"""Validate the 2D flat-plate thermal-conduction run against the semi-infinite
error-function solution.

A hot quiescent gas (T_inf) is suddenly put in contact with a cold isothermal
wall (T_wall) at y = 0. For pure 1D conduction in a semi-infinite medium the
self-similar solution is

    (T(y,t) - T_wall) / (T_inf - T_wall) = erf( y / (2*sqrt(alpha*t)) ),
    alpha = k / (rho*cp).

This script reads the conserved-variable restart files directly ([var, y, x]
global layout), derives the EOS temperature, and overlays the erf solution.
Grid / dt / EOS are read from simulation.inp so the script stays correct if
the case is retuned.

    ./mfc.sh run examples/2D_Thermal_Conduction_Flatplate/case.py -n 16
    python3 examples/2D_Thermal_Conduction_Flatplate/validate.py
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
from scipy.special import erf

HERE = os.path.dirname(os.path.abspath(__file__))


def read_inp(path):
    inp = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            inp[k.strip().lower()] = v.strip().rstrip(",")
    return inp


def main():
    inp = read_inp(os.path.join(HERE, "simulation.inp"))
    m, n = int(inp["m"]), int(inp["n"])
    ncx, ncy = m + 1, n + 1
    dt = float(inp["dt"])
    Lx = float(inp["x_domain%end"]) - float(inp["x_domain%beg"])
    Ly = float(inp["y_domain%end"]) - float(inp["y_domain%beg"])

    # Stiffened-gas EOS storage form: gamma_p = 1/(g-1), pi_inf_p = g*p_inf/(g-1)
    gamma_p = float(inp["fluid_pp(1)%gamma"])
    pi_inf_p = float(inp["fluid_pp(1)%pi_inf"])
    cv = float(inp["fluid_pp(1)%cv"])
    k_therm = float(inp["fluid_pp(1)%k_therm"])
    g = 1.0 + 1.0 / gamma_p  # recover ratio of specific heats
    R_gas = (g - 1.0) * cv  # specific gas constant
    cp = g * cv
    T_wall = float(inp["bc_y%twall_in"])

    dy = Ly / ncy
    yc = (np.arange(ncy) + 0.5) * dy  # cell-center distance from the wall at y=0

    rd = os.path.join(HERE, "restart_data")
    steps = sorted(int(g_.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (g_ := re.search(r"lustre_(\d+)\.dat$", os.path.basename(f))))

    def temperature(step):
        a = np.fromfile(os.path.join(rd, f"lustre_{step}.dat"), np.float64)
        rec = a.reshape(a.size // (ncx * ncy), ncy, ncx)  # [var, y, x]
        rho, mx, my, E = rec[0], rec[1], rec[2], rec[3]
        rho_e = E - 0.5 * (mx**2 + my**2) / rho
        p = (rho_e - pi_inf_p) / gamma_p
        T = (p + pi_inf_p * (g - 1.0) / g) / (rho * R_gas)  # p_inf physical = pi_inf_p*(g-1)/g
        umag = np.sqrt(mx**2 + my**2) / rho
        return T, umag, p

    T0field, _, _ = temperature(0)
    T_inf = float(T0field.mean())  # uniform free-stream from the IC itself
    rho_inf = float(np.fromfile(os.path.join(rd, "lustre_0.dat"), np.float64).reshape(-1, ncy, ncx)[0].mean())
    alpha_inf = k_therm / (rho_inf * cp)  # free-stream diffusivity
    T_film = 0.5 * (T_inf + T_wall)  # film temperature
    rho_film = rho_inf * (T_inf / T_film)  # constant pressure -> rho ~ 1/T
    alpha_film = k_therm / (rho_film * cp)  # film-temperature diffusivity (the matching one)

    T_last, u_last, _ = temperature(steps[-1])
    t_last = steps[-1] * dt
    xc = (np.arange(ncx) + 0.5) * (Lx / ncx)

    # field comparison: MFC T(x,y) | analytic erf field | difference
    # The analytic solution is 1D in y; tile it across x onto the same grid.
    theta_col = erf(yc / (2.0 * np.sqrt(alpha_film * t_last)))
    T_an_field = (T_wall + (T_inf - T_wall) * theta_col)[:, None] * np.ones((1, ncx))
    diff = T_last - T_an_field

    fig1, ax = plt.subplots(1, 3, figsize=(15, 4.4), layout="constrained")
    im0 = ax[0].pcolormesh(xc * 1e3, yc * 1e3, T_last, vmin=T_wall, vmax=T_inf, cmap="inferno", shading="auto")
    ax[0].set(title=f"MFC  T(x,y)   t={t_last * 1e3:.2f} ms", xlabel="x [mm]", ylabel="y [mm]")
    im1 = ax[1].pcolormesh(xc * 1e3, yc * 1e3, T_an_field, vmin=T_wall, vmax=T_inf, cmap="inferno", shading="auto")
    ax[1].set(title="analytic  T_wall+(T_inf−T_wall)·erf(y/2√(α_film·t))", xlabel="x [mm]", ylabel="y [mm]")
    fig1.colorbar(im1, ax=[ax[0], ax[1]], shrink=0.85, label="T [K]")
    dmax = float(np.abs(diff).max())
    im2 = ax[2].pcolormesh(xc * 1e3, yc * 1e3, diff, vmin=-dmax, vmax=dmax, cmap="coolwarm", shading="auto")
    ax[2].set(title=f"MFC − analytic  (max |Δ| = {dmax:.1f} K)", xlabel="x [mm]", ylabel="y [mm]")
    fig1.colorbar(im2, ax=ax[2], shrink=0.85, label="ΔT [K]")
    for a in ax:
        a.set_ylim(0, 6)  # zoom into the near-wall thermal layer
    fig1.savefig(os.path.join(HERE, "validation_field.png"), dpi=130)

    # profiles vs erf, averaged over the x-interior (away from inflow/outflow)
    i0, i1 = ncx // 4, 3 * ncx // 4  # central half in x (away from inflow/outflow)

    fig2, ax = plt.subplots(1, 2, figsize=(13, 5))
    theta_an = lambda y, t, al: erf(y / (2.0 * np.sqrt(al * t)))
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(steps) - 1))
    print(f"  T_inf={T_inf:.1f} K  T_wall={T_wall:.1f} K  alpha_inf={alpha_inf:.3e}  alpha_film={alpha_film:.3e}")
    print(f"  {'t [ms]':>8} {'delta_mfc[mm]':>13} {'rms err vs erf(alpha_film) [K]':>32}")
    rms_all = []
    for c, s in zip(colors, steps[1:]):
        T, _, _ = temperature(s)
        prof = T[:, i0:i1].mean(axis=1)
        theta = (prof - T_wall) / (T_inf - T_wall)
        t = s * dt
        ax[0].plot(theta, yc * 1e3, color=c, lw=1.3)
        # erf reference using the film-temperature diffusivity
        th_ref = theta_an(yc, t, alpha_film)
        rms = math.sqrt(np.mean((theta - th_ref) ** 2)) * (T_inf - T_wall)
        rms_all.append(rms)
        # thermal-layer thickness: where theta crosses 0.99
        d99 = np.interp(0.99, theta, yc)
        print(f"  {t * 1e3:8.3f} {d99 * 1e3:13.3f} {rms:32.2f}")
    # one erf reference curve at the final time for the legend
    ax[0].plot(theta_an(yc, t_last, alpha_film), yc * 1e3, "k--", lw=2, label=f"erf, α_film={alpha_film:.2e}")
    ax[0].set(xlabel="θ = (T−T_wall)/(T_inf−T_wall)", ylabel="y [mm]", title="T profiles (color=time) vs erf", ylim=(0, 6))
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    # self-similar collapse: T vs eta = y/(2 sqrt(alpha_film t))
    for c, s in zip(colors, steps[1:]):
        T, _, _ = temperature(s)
        prof = T[:, i0:i1].mean(axis=1)
        theta = (prof - T_wall) / (T_inf - T_wall)
        t = s * dt
        eta = yc / (2.0 * np.sqrt(alpha_film * t))
        ax[1].plot(eta, theta, color=c, lw=1.0, alpha=0.8)
    eta_ref = np.linspace(0, 3, 200)
    ax[1].plot(eta_ref, erf(eta_ref), "k--", lw=2, label="erf(η)")
    ax[1].set(xlabel="η = y / (2√(α_film t))", ylabel="θ", title="self-similar collapse", xlim=(0, 3))
    ax[1].legend()
    ax[1].grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(HERE, "validation_profiles.png"), dpi=130)

    Tspread = T_last.max(axis=1) - T_last.min(axis=1)  # 1D-ness check (printed, not plotted)
    print(f"\n  max|u| at final time = {u_last.max():.3f} m/s")
    print(f"  peak x-inhomogeneity  = {Tspread.max():.2f} K")
    print(f"  mean rms error vs erf(alpha_film) = {np.mean(rms_all):.2f} K ({100 * np.mean(rms_all) / (T_inf - T_wall):.1f}% of ΔT)")
    out = {
        "flatplate": {
            "N": ncx,
            "T_inf": T_inf,
            "T_wall": T_wall,
            "alpha_film": alpha_film,
            "mean_rms_vs_erf_K": float(np.mean(rms_all)),
            "mean_rms_pct_of_dT": float(100 * np.mean(rms_all) / (T_inf - T_wall)),
            "max_u_final": float(u_last.max()),
        }
    }
    open(os.path.join(HERE, "summary.json"), "w").write(json.dumps(out, indent=2) + "\n")
    print("  wrote validation_field.png, validation_profiles.png, summary.json")


if __name__ == "__main__":
    main()
