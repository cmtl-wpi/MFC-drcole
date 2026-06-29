#!/usr/bin/env python3
"""All thermocapillary plotting in one tool. Pick a subcommand:

  samareh   The two headline validation overlays:
              figures/case1_fig5.png  -- TC1 Fig 5 (Ma=0), v/v_YGB vs t/t_r vs Samareh VOF
              figures/case2_fig7.png  -- TC2 Fig 7, U* vs t* vs the Nas & Tryggvason transient
            Raw per-snapshot markers, no smoothing (the compressible box's acoustic ring is left visible).

  ma        figures/tc1_ma_convergence.png -- TC1 terminal v_t/v_YGB vs a conduction Marangoni-number
            sweep (smaller Ma -> closer to the invariant-T limit), with a linear extrapolation to Ma=0.
            Calls measure.py per run dir; skips runs that have no snapshots yet.

  fields    A derived field from one run (writes <case_dir>/viz/<field>_<step>.png):
              temperature    EOS-recovered T field + centerline profile vs the frozen initial linear T
              sigma          sigma(T) field + sigma along the interface vs angle (the Marangoni driver)
              recirculation  drop-frame streamlines (colored by speed) + cell-resolved vorticity
            Temperature is not stored by MFC; it is recovered per cell from the stiffened-gas EOS
                T = (p + p_inf) / ((gamma - 1) * rho * cv),    p from the conserved internal energy.

All run-dependent constants are read from each run's simulation.inp, so the plots can't silently
disagree with the data. Conserved layout (model_eqns=3, num_fluids=2): 0,1 = partial densities,
2 = x-momentum, 3 = y-momentum, 7,8 = phasic internal energies, color c last.

Usage (no subcommand runs `samareh`, the headline overlays -- the usual rebuild):
    python3 plot.py [samareh]
    python3 plot.py ma
    python3 plot.py fields [case_dir] [temperature|sigma|recirculation] [step]
        (fields' step defaults to the last snapshot; recirculation's 3rd arg is a t/tau target.)
    python3 plot.py clean [--force]
        (remove orphaned figures/ images no current script produces; dry-run unless --force)
"""

import glob
import json
import os
import re
import subprocess
import sys
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIGS = os.path.join(HERE, "figures")
R = 0.5  # drop radius (D = 1)
GRADT = 2.0 / 15.0  # imposed |dT/dy|, common to TC1/TC2


def read_namelist(path):
    """Parse a Fortran namelist file's plain "name = value" lines into a dict (lowercase keys)."""
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


# samareh: the two headline validation overlays


def color_weighted_vy(run_dir):
    """Per-snapshot color-weighted lab-frame y-velocity history of a slip-wall run.
    Returns (t, u_lab, params) or None. The drop migrates in +y, so u_lab IS the rise velocity."""
    inp = os.path.join(run_dir, "simulation.inp")
    rd = os.path.join(run_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    P = read_namelist(inp)
    f = lambda k: float(P[k.lower()])  # noqa: E731
    nx, ny, nz = int(f("m")) + 1, int(f("n")) + 1, int(f("p")) + 1
    cells = nx * ny * nz
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    c_idx = nvars - 1  # color function is the last conserved variable
    t, u_lab = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)  # noqa: E731
        vy = fld(3) / (fld(0) + fld(1))
        c = np.clip(fld(c_idx), 0.0, None)
        t.append(s * f("dt"))
        u_lab.append((c * vy).sum() / c.sum())
    return np.array(t), np.array(u_lab), P


def v_ygb_ratio(out):
    """(t/t_r, v/v_YGB) from a color_weighted_vy result, using each run's own constants."""
    t, u_lab, P = out
    mu = 1.0 / float(P["fluid_pp(1)%re(1)"])
    dsdt = float(P["sigma_dtdt"])
    t_r = mu / abs(dsdt * GRADT)
    v_YGB = (2.0 / 15.0) * (-dsdt) * GRADT * R / mu
    return t / t_r, u_lab / v_YGB


# Nas & Tryggvason U*(t*) transient, digitized BY EYE from Samareh Fig 7 (the red open triangles; the
# two Samareh grids nearly coincide with it). Accuracy ~ +/-0.005 in U*. Anchors match the paper text:
# broad peak ~0.131 at t*~4-5, terminal ~0.10 at t*=20 (the fine grid is within 1.7% of N&T).
NAS_TRYGGVASON = np.array(
    [
        (0.0, 0.0),
        (1.0, 0.055),
        (2.0, 0.100),
        (3.0, 0.122),
        (4.0, 0.130),
        (5.0, 0.131),
        (6.0, 0.128),
        (7.0, 0.124),
        (8.0, 0.120),
        (10.0, 0.114),
        (12.0, 0.110),
        (14.0, 0.106),
        (16.0, 0.103),
        (18.0, 0.101),
        (20.0, 0.0995),
    ]
)

# Samareh Fig 5(d) VOF curve (sharp-interface analogue of MFC), digitized by eye from the published
# raster (~ +/-0.02 in v/v_YGB); his invariant-T plateau holds flat ~0.82-0.84 out to t/t_r = 10.
SAMAREH_VOF = np.array(
    [
        (0.0, 0.0),
        (0.12, 0.42),
        (0.28, 0.70),
        (0.45, 0.80),
        (0.7, 0.815),
        (1.0, 0.82),
        (2.0, 0.825),
        (3.0, 0.83),
        (4.0, 0.83),
        (5.0, 0.835),
        (6.0, 0.83),
        (7.0, 0.835),
        (8.0, 0.838),
        (9.0, 0.835),
        (10.0, 0.84),
    ]
)

# Seaborn aesthetic (whitegrid + notebook context). Applied via rc_context in each plot; the
# explicit per-curve colors set below (frozen reds, conduction blues) are preserved on top of it.
PLATE_STYLE = {
    **sns.axes_style("ticks"),
    **sns.plotting_context("paper", font_scale=1.3),
    "figure.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
}


def read_case_Ma(run_dir):
    """Read the `Ma = <float>` line from the case copy in a run dir. Returns None when there is
    no such line (the frozen-T case has no Marangoni number), so the label can say Ma=0."""
    for f in sorted(glob.glob(os.path.join(run_dir, "case_*.py"))):
        for line in open(f):
            if line.strip().startswith("Ma ="):
                try:
                    return float(line.split("=", 1)[1].split("#")[0])
                except ValueError:
                    return None
    return None


def samareh_fig5():
    """TC1, Fig 5: grid-convergence sweeps of both the frozen-T (Ma=0) and bulk-conduction (Ma=0.1)
    cases, each at several cells/D, plotted against Samareh's VOF curve.

    Frozen-T pins the imposed linear T field through density (no energy transport); conduction evolves
    the coupled energy equation. Each case is its own color family (frozen reds, conduction blues),
    shaded light->dark = coarse->fine, so the rise/overshoot and the late-time droop can be read as a
    function of resolution. Each curve labels its own cells/D and run length, read from its case copy.
    """
    # (run dir, color); missing runs are skipped. box width = 5D, so cells/D = (m+1)/5, read per run.
    # Grid convergence for BOTH cases: frozen-T (Ma=0) in reds, conduction (Ma=0.1) in blues; within
    # each family light->dark = coarse->fine. Built over the grid ladder so new sweep points (the
    # sweep_grid_2d campaign: w096/w192/w384) appear automatically without editing this list.
    grid_nx = [64, 96, 128, 192, 256, 384]
    shades = np.linspace(0.35, 0.95, len(grid_nx))
    runs = [(f"tc1/frozen/w{nx:03d}", plt.cm.Reds(s)) for nx, s in zip(grid_nx, shades)]
    runs += [(f"tc1/ma0p1/w{nx:03d}/sc050", plt.cm.Blues(s)) for nx, s in zip(grid_nx, shades)]
    series = []
    for name, color in runs:
        run = os.path.join(RUNS, name)
        if not os.path.isdir(os.path.join(run, "restart_data")):
            print(f"  fig5: {name} not found, skipping")
            continue
        out = color_weighted_vy(run)
        if out is None or len(out[0]) < 5:
            print(f"  fig5: {name} not ready ({0 if out is None else len(out[0])} snapshots), skipping")
            continue
        x, y = v_ygb_ratio(out)
        cells_per_D = (int(out[2]["m"]) + 1) / 5.0
        Ma = read_case_Ma(run)
        kind = "frozen-T" if Ma is None else "conduction"
        ma_txt = r"$Ma=0$" if Ma is None else rf"$Ma={Ma:g}$"
        label = rf"{kind} ({ma_txt}), {cells_per_D:.1f}/$D$"
        series.append((x, y, color, label))
    if not series:
        print("  fig5: no runs found")
        return
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(SAMAREH_VOF[:, 0], SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0, mfc="none", mew=1.3, lw=1.0, label=r"Samareh Fig 5(d), VOF ($Ma=0$, digitized)")
        ax.axhline(1.0, color="0.3", lw=1.1, ls="--", zorder=1, label=r"$u_{\mathrm{YGB}}$ (analytic terminal, $\approx 8.89{\times}10^{-3}$)")
        for x, y, color, label in series:
            ax.plot(x, y, "-", color=color, lw=1.7, alpha=0.95, solid_capstyle="round", label=label)
        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$t / t_r$")
        ax.set_ylabel(r"rise velocity   $u / u_{\mathrm{YGB}}$")
        ax.set_title(r"2D thermocapillary rise: grid convergence", fontsize=13, loc="left")
        ax.legend(loc="lower left", fontsize=9, frameon=False, ncol=2, columnspacing=1.2, handlelength=1.6)
        sns.despine(ax=ax)
        dst = os.path.join(FIGS, "case1_fig5.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(series)} runs)")


def samareh_fig7():
    """TC2, Fig 7 (Re=5, Ma=20, Ca=0.01666): MFC migration vs the digitized Nas & Tryggvason transient."""
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(NAS_TRYGGVASON[:, 0], NAS_TRYGGVASON[:, 1], "^--", color="0.0", ms=6.5, mfc="none", mew=1.3, lw=1.0, zorder=5, label="Nas & Tryggvason (digitized)")
        plotted = False
        for name, nx, color in [("tc2/w064", 64, "#4C72B0"), ("tc2/w128", 128, "#DD8452")]:
            out = color_weighted_vy(os.path.join(RUNS, name))
            if out is None or len(out[0]) < 10:  # skip absent or depleted runs (need a real curve)
                if out is not None:
                    print(f"  fig7: skipping {name} -- only {len(out[0])} snapshots on disk (data depleted)")
                continue
            t, u_lab, params = out

            # Build the Marangoni reference scales (velocity U_r, time t_r) from this run's constants.
            mu_b = 1.0 / float(params["fluid_pp(1)%re(1)"])
            sigma_T = float(params["sigma_dtdt"])
            Ly = float(params["y_domain%end"]) - float(params["y_domain%beg"])
            if "bc_y%twall_out" in params:
                gradT = abs(float(params["bc_y%twall_out"]) - float(params["bc_y%twall_in"])) / Ly
            else:
                gradT = 1.0 / Ly
            marangoni_stress = abs(sigma_T * gradT)
            U_r = marangoni_stress * R / mu_b
            t_r = mu_b / marangoni_stress
            ts, us = t / t_r, u_lab / U_r

            ax.plot(ts, us, "o--", color=color, ms=3.5, mew=0, lw=0.9, alpha=0.6, zorder=3, label=f"MFC {nx} cells/width ({nx // 2}/$D$)")
            plotted = True
        if not plotted:
            print("  fig7: no runs found")
            plt.close(fig)
            return
        ax.axhline(0.0, color="0.75", lw=0.8, zorder=1)  # rest baseline (raw scatter dips below it)
        ax.set_xlim(0.0, 20.0)
        ax.set_ylim(-0.025, 0.15)  # slightly past Samareh's 0-0.14 so the raw acoustic scatter is visible
        ax.set_xlabel(r"$t^* = t/t_r$")
        ax.set_ylabel(r"$U^* = U/U_r$")
        ax.set_title("Fig 7 — 2D migration at finite $Ma$ (Re=5, Ma=20, Ca=0.0167)", fontsize=12)
        ax.legend(loc="upper right", fontsize=10, frameon=False)
        dst = os.path.join(FIGS, "case2_fig7.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}")


def cmd_samareh(argv):
    os.makedirs(FIGS, exist_ok=True)
    samareh_fig5()
    samareh_fig7()
    # Also produce the temperature centerline figure for any available run.
    # Runs nest under runs/ by axis (case/Ma/grid/smoothing), so walk for leaf dirs.
    leaf_runs = sorted(d for d, _, files in os.walk(RUNS) if "simulation.inp" in files)
    for run_dir in leaf_runs:
        name = os.path.relpath(run_dir, RUNS)
        rd = os.path.join(run_dir, "restart_data")
        if not os.path.isdir(rd):
            continue
        steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
        if len(steps) < 2:
            continue
        print(f"\n--- temperature field: {name} ---")
        c = load_field_case(run_dir)
        field_temperature(c, steps[-1])
        break  # only the first available run


# ma: TC1 conduction Ma -> 0 convergence sweep

SAMAREH_RATIO = 0.80  # Samareh's converged 2D cylinder ratio (Fig 5)

# (Ma, run dir) -- the conduction Ma -> 0 convergence sweep (all w128, slip-wall, tr=2.0).
# Nested-tree paths; these w128 Ma points are not yet on disk (aspirational sweep).
MA_SWEEP = [
    (0.30, "tc1/ma0p30/w128/sc050"),
    (0.10, "tc1/ma0p10/w128/sc050"),
    (0.05, "tc1/ma0p05/w128/sc050"),
    (0.03, "tc1/ma0p03/w128/sc050"),
]


def measure_run(run_dir):
    """Run measure.py on a run dir; return its RESULT_JSON dict, or None if not ready."""
    proc = subprocess.run(
        [sys.executable, os.path.join(HERE, "measure.py"), run_dir],
        capture_output=True,
        text=True,
        check=False,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    return None


def cmd_ma(argv):
    points = []
    for ma, name in MA_SWEEP:
        run_dir = os.path.join(RUNS, name)
        result = measure_run(run_dir) if os.path.isdir(run_dir) else None
        if result is None:
            print(f"Ma={ma:<5} {name}: no result yet")
            continue
        points.append((ma, result["ratio_plateau"]))
        print(
            f"Ma={ma:<5} plateau={result['ratio_plateau']:+.3f}  final={result['ratio_final']:+.3f}  drift={result['slope_per_tr']:+.3f}/t_r  ran {result['t_end_tr']:.2f} t_r  rises={result['rises']}"
        )
    if not points:
        sys.exit("no runs measured yet -- let the sweep produce snapshots first")

    points.sort()  # ascending Ma
    ma = np.array([p[0] for p in points])
    plateau = np.array([p[1] for p in points])

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.axhline(SAMAREH_RATIO, ls="--", color="C3", lw=1.4, label=f"Samareh 2D = {SAMAREH_RATIO:.2f}")
    ax.plot(ma, plateau, "o-", color="C0", ms=7, lw=1.6, label=r"MFC conduction (plateau $v_t/v_{\mathrm{YGB}}$)")

    # Linear extrapolation to Ma=0 from the two smallest-Ma points (the near-invariant end).
    if len(ma) >= 2:
        lo = np.argsort(ma)[:2]
        slope, intercept = np.polyfit(ma[lo], plateau[lo], 1)
        ax.plot([0, ma.max()], [intercept, slope * ma.max() + intercept], ":", color="0.5", lw=1.2, label=rf"linear extrap. $\to$ Ma=0: {intercept:.3f}")
        ax.plot(0, intercept, "*", color="k", ms=15, zorder=5)

    ax.set_xlabel(r"thermal Marangoni number  $Ma$   (smaller $\to$ closer to invariant $T$)")
    ax.set_ylabel(r"terminal  $v_t / v_{\mathrm{YGB}}$")
    ax.set_xlim(left=-0.012)
    ax.set_title(r"TC1 (Samareh 4.1.1): conduction $Ma\to0$ convergence to the invariant-$T$ limit")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    os.makedirs(FIGS, exist_ok=True)
    out = os.path.join(FIGS, "tc1_ma_convergence.png")
    fig.savefig(out, dpi=150)
    print("saved ->", out)


# fields: a derived field from one run


def load_field_case(case_dir):
    """Read the grid, EOS constants, and snapshot list for one run into a namespace."""
    P = read_namelist(os.path.join(case_dir, "simulation.inp"))
    p = lambda k: float(P[k.lower()])  # noqa: E731
    nx, ny = int(p("m")) + 1, int(p("n")) + 1
    # MFC stores gamma_mfc = 1/(gamma-1) and pi_inf_mfc = gamma*p_inf/(gamma-1); invert to physical.
    gamma_mfc, pi_inf_mfc = p("fluid_pp(1)%gamma"), p("fluid_pp(1)%pi_inf")
    gamma = 1.0 + 1.0 / gamma_mfc
    p_inf = pi_inf_mfc * (gamma - 1.0) / gamma
    restart_dir = os.path.join(case_dir, "restart_data")
    xb = np.fromfile(os.path.join(restart_dir, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
    yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
    steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(restart_dir, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
    if not steps:
        sys.exit(f"no lustre_<step>.dat snapshots in {restart_dir!r} -- run the case first")
    ncell = nx * ny
    nvars = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size // ncell
    viz_dir = os.path.join(case_dir, "viz")
    os.makedirs(viz_dir, exist_ok=True)
    return SimpleNamespace(
        P=P,
        p=p,
        nx=nx,
        ny=ny,
        ncell=ncell,
        dt=p("dt"),
        cv=p("fluid_pp(1)%cv"),
        gamma_mfc=gamma_mfc,
        pi_inf_mfc=pi_inf_mfc,
        gamma=gamma,
        p_inf=p_inf,
        restart_dir=restart_dir,
        xb=xb,
        yb=yb,
        x=0.5 * (xb[:-1] + xb[1:]),
        y=0.5 * (yb[:-1] + yb[1:]),
        steps=steps,
        nvars=nvars,
        c_idx=nvars - 1,
        viz_dir=viz_dir,
    )


def columns_of(c, step):
    """A snapshot as (nvars, ncell): row i is conserved variable i flattened."""
    return np.fromfile(os.path.join(c.restart_dir, f"lustre_{step}.dat"), np.float64).reshape(-1, c.ncell)


def eos_temperature(c, cols):
    """(T, color) on the (ny, nx) grid from a (nvars, ncell) snapshot, via the stiffened-gas EOS."""
    rho = (cols[0] + cols[1]).reshape(c.ny, c.nx)
    rho_e = (cols[7] + cols[8]).reshape(c.ny, c.nx)  # phasic internal energies (no kinetic part)
    p = (rho_e - c.pi_inf_mfc) / c.gamma_mfc
    T = (p + c.p_inf) / ((c.gamma - 1.0) * rho * c.cv)
    color = np.clip(cols[c.c_idx].reshape(c.ny, c.nx), 0.0, None)
    return T, color


def field_temperature(c, step):
    """EOS temperature field + centerline profile along the droplet rise direction (y)."""
    T, color = eos_temperature(c, columns_of(c, step))
    T_initial, _ = eos_temperature(c, columns_of(c, c.steps[0]))
    if not np.all(np.isfinite(T)) or T.min() <= 0.0:
        sys.exit(f"reconstructed T is non-physical (min = {T.min():.4f}); constants likely mismatch the data.")
    print(f"step {step}  (t = {step * c.dt:.4f})   T range: [{T.min():.4f}, {T.max():.4f}]")
    x, y = c.x, c.y
    fig, (ax_field, ax_line) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"width_ratios": [1.05, 1.0]})
    vmin, vmax = min(T_initial.min(), T.min()), max(T_initial.max(), T.max())
    mesh = ax_field.pcolormesh(x, y, T, cmap="coolwarm", vmin=vmin, vmax=vmax, shading="auto")
    fig.colorbar(mesh, ax=ax_field, label="temperature $T$")
    ax_field.contour(x, y, color, levels=[0.5], colors="k", linewidths=1.2)
    xc = (color * x[None, :]).sum() / color.sum()
    yc = (color * y[:, None]).sum() / color.sum()
    ax_field.plot(xc, yc, "k+", ms=10, mew=2)
    ax_field.set(aspect="equal", xlabel="$x$", ylabel="$y$", title=f"Temperature field (step {step}, $t={step * c.dt:.3f}$)\ndroplet centroid $x={xc:+.4f}$, $y={yc:+.4f}$")
    # Centerline along the rise direction (y): T vs y at x ≈ 0 (domain center)
    mid = c.nx // 2
    ax_line.plot(T_initial[:, mid], y, ":", color="C0", alpha=0.7, label=f"centerline $x\\approx 0$, step {c.steps[0]}")
    ax_line.plot(T[:, mid], y, "-", color="C3", label=f"centerline $x\\approx 0$, step {step}")
    ax_line.set(xlabel=r"temperature $T$", ylabel="$y$ (rise direction)", title="Centerline temperature profile\n(along rise direction at $x \\approx 0$)")
    ax_line.legend(loc="best", fontsize=8)
    ax_line.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(FIGS, exist_ok=True)
    out = os.path.join(FIGS, f"temperature_{step}.png")
    fig.savefig(out, dpi=150)
    print(f"saved figure -> {out}")


def field_sigma(c, step):
    """sigma(T) field + sigma in the interface cells vs angle (the Marangoni driver)."""
    sigma0, dsigma_dT, T_ref = c.p("sigma"), c.p("sigma_dtdt"), c.p("sigma_t_ref")
    T, color = eos_temperature(c, columns_of(c, step))
    sigma = sigma0 + dsigma_dT * (T - T_ref)
    x, y = c.x, c.y
    xc = (color * x[None, :]).sum() / color.sum()
    yc = (color * y[:, None]).sum() / color.sum()
    X, Y = np.meshgrid(x, y)
    interface = (color > 0.2) & (color < 0.8)
    angle = np.degrees(np.arctan2(Y[interface] - yc, X[interface] - xc))
    sig_if = sigma[interface]
    fig, (ax_field, ax_line) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.05, 1]})
    mesh = ax_field.pcolormesh(x, y, sigma, cmap="viridis", shading="auto")
    fig.colorbar(mesh, ax=ax_field, label=r"$\sigma(T)$")
    ax_field.contour(x, y, color, levels=[0.5], colors="w", linewidths=1.2)
    ax_field.set(aspect="equal", xlabel="x", ylabel="y", title=rf"$\sigma(T)$ field + interface (step {step}, t = {step * c.dt:.2f})")
    ax_line.scatter(angle, sig_if, s=12, color="C0")
    ax_line.axvline(0, ls=":", color="C3")
    ax_line.text(6, sig_if.min(), "hot (+x)", color="C3", fontsize=9)
    ax_line.set(xlabel="angle around interface (deg)", ylabel=r"$\sigma$ in the interface cells", title=r"Low $\sigma$ on the hot side drives the Marangoni pull")
    ax_line.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(c.viz_dir, f"sigma_interface_{step}.png")
    fig.savefig(out, dpi=150)
    print(f"sigma on interface spans [{sig_if.min():.4f}, {sig_if.max():.4f}]")
    print(f"saved -> {out}")


def field_recirculation(c, target_ttau):
    """Drop-frame streamlines (colored by speed) + cell-resolved vorticity, near the drop."""
    x, y, xb, yb = c.x, c.y, c.xb, c.yb
    mu = 1.0 / c.p("fluid_pp(1)%re(1)")
    dsigma_dT = c.p("sigma_dtdt")
    r, gradT = 0.5, 2.0 / 15.0
    v_YGB = (2.0 / 15.0) * (-dsigma_dT) * gradT * r / mu

    def snap_uvc(step):
        cols = columns_of(c, step)
        rho = (cols[0] + cols[1]).reshape(c.ny, c.nx)
        return rho, (cols[2].reshape(c.ny, c.nx)) / rho, (cols[3].reshape(c.ny, c.nx)) / rho, np.clip(cols[c.c_idx].reshape(c.ny, c.nx), 0.0, 1.0)

    rho0, _, _, c0 = snap_uvc(c.steps[0])
    rho_drop = (c0 * rho0).sum() / c0.sum()
    tau = rho_drop * r**2 / mu  # viscous time from the t=0 drop density (tracks the case, not rho=1)
    step = min(c.steps, key=lambda s: abs(s * c.dt / tau - target_ttau))
    ttau = step * c.dt / tau
    _, u, v, col = snap_uvc(step)
    u_drop = (col * u).sum() / col.sum()
    v_drop = (col * v).sum() / col.sum()
    uc, vc = u - u_drop, v - v_drop  # co-moving frame
    omega = np.gradient(vc, x, axis=1) - np.gradient(uc, y, axis=0)
    xc_drop = (col * x[None, :]).sum() / col.sum()
    yc_drop = (col * y[:, None]).sum() / col.sum()
    win = 2.5 * r
    mx = (x > xc_drop - win) & (x < xc_drop + win)
    my = (y > yc_drop - win) & (y < yc_drop + win)
    xz, yz = x[mx], y[my]
    uz, vz = uc[np.ix_(my, mx)], vc[np.ix_(my, mx)]
    cz, oz = col[np.ix_(my, mx)], omega[np.ix_(my, mx)]
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
        for col_, lw in (("white", 4.0), (inkc, 1.8)):
            ax.annotate("", xy=(xc_drop, yc_drop + 1.0 * r), xytext=(xc_drop, yc_drop - 0.05 * r), arrowprops=dict(arrowstyle="-|>", color=col_, lw=lw, shrinkA=0, shrinkB=0))
        ax.text(0.04, 0.04, "cold", transform=ax.transAxes, color=coldc, fontsize=8, ha="left", va="bottom", path_effects=halo)
        ax.text(0.04, 0.96, "hot", transform=ax.transAxes, color=hotc, fontsize=8, ha="left", va="top", path_effects=halo)
    axa.set_ylabel(r"$y/D$")
    axa.text(xc_drop + 0.28 * r, yc_drop + 0.15 * r, rf"$U={v_drop / v_YGB:.2f}\,v_{{\rm YGB}}$", color=inkc, fontsize=9, ha="left", va="bottom", path_effects=halo)
    fig.suptitle(rf"2D thermocapillary drop  $\cdot$  {c.ny / 7.5:.1f} cells/$D$  $\cdot$  $t={ttau:.1f}\,\tau$", fontsize=9.5)
    os.makedirs(FIGS, exist_ok=True)
    out = os.path.join(FIGS, "case1_zero_marangoni_2D_recirculation")
    fig.savefig(out + ".png", dpi=300)
    fig.savefig(out + ".pdf")
    print(f"step={step}  t/tau={ttau:.3f}  nvars={c.nvars}  rho_drop={rho_drop:.3f}  tau={tau:.3f}")
    print(f"U_drop/v_YGB={v_drop / v_YGB:+.3f}  (lateral u_drop/v_YGB={u_drop / v_YGB:+.3f})")
    print(f"saved -> {out}.png / .pdf")


def cmd_fields(argv):
    case_dir = argv[0] if len(argv) > 0 else HERE
    field = argv[1].lower() if len(argv) > 1 else "temperature"
    if field not in ("temperature", "sigma", "recirculation"):
        sys.exit(f"unknown field {field!r}; choose temperature, sigma, or recirculation")
    c = load_field_case(case_dir)
    if field == "recirculation":
        field_recirculation(c, float(argv[2]) if len(argv) > 2 else 2.6)
        return
    step = int(argv[2]) if len(argv) > 2 else c.steps[-1]
    if step not in c.steps:
        sys.exit(f"step {step} not available; choose from {c.steps[0]}..{c.steps[-1]}")
    (field_temperature if field == "temperature" else field_sigma)(c, step)


# clean: remove orphaned figures (figures no current script produces)

# Figures the scripts here legitimately produce -- KEEP IN SYNC with the savefig calls above and with
# compare_tc3_visc.py. `clean` deletes any other .png/.pdf in figures/ as a stale orphan; source files
# (.tex, etc.) and PRECIOUS run data (restart_data/case.py/*.inp) are never touched.
KEEP_FIGURES = {
    "case1_fig5.png",  # samareh
    "case2_fig7.png",  # samareh
    "tc1_ma_convergence.png",  # ma
    "case1_zero_marangoni_2D_recirculation.png",  # fields recirculation
    "case1_zero_marangoni_2D_recirculation.pdf",
    "case3_large_marangoni_mu_of_T_validation.png",  # compare_tc3_visc.py
    "mechanism_schematic.png",
    "mechanism_schematic.pdf",  # TikZ schematic (source: mechanism_schematic.tex)
}
KEEP_FIGURE_PATTERNS = [re.compile(r"^temperature_\d+\.png$")]  # fields temperature, one per step


def cmd_clean(argv):
    """Remove orphaned figures from figures/ -- ones no current script produces. Dry-run unless --force.
    Only .png/.pdf are candidates; source files and run data are out of scope."""
    force = "--force" in argv
    if not os.path.isdir(FIGS):
        print("figures/: nothing to clean (directory absent)")
        return
    orphans = [
        f for f in sorted(os.listdir(FIGS)) if os.path.isfile(os.path.join(FIGS, f)) and f.endswith((".png", ".pdf")) and f not in KEEP_FIGURES and not any(p.match(f) for p in KEEP_FIGURE_PATTERNS)
    ]
    if not orphans:
        print("figures/: no orphaned figures")
        return
    for f in orphans:
        if force:
            os.remove(os.path.join(FIGS, f))
        print(f"  {'removed' if force else 'would remove'}  figures/{f}")
    print(f"\n{len(orphans)} orphaned figure(s) " + ("removed" if force else "-- re-run with `plot.py clean --force` to delete"))


# recon: late-droop sensitivity to the reconstruction scheme (fixed grid, Ma=0.1, 12.8 cells/D).
# The droop is interface-band smearing during advection, so a less-diffusive / interface-compressing
# scheme reduces it -- shown here at fixed dx (no grid refinement), isolating the scheme's diffusion.
RECON_RUNS = [
    ("recon/muscl", "MUSCL (Van Leer)", "#dd8452"),  # 2nd-order, most diffusive
    ("recon/weno5", "WENO5 (baseline)", "#4c72b0"),
    ("recon/wenoz", "WENO-Z", "#8172b3"),
    ("recon/weno7", "WENO7", "#55a868"),  # higher order, less diffusive
    ("recon/muscl_thinc", "MUSCL + THINC (int_comp)", "#c44e52"),  # active interface compression
]


def samareh_recon():
    """Late-time droop vs reconstruction scheme: the Ma=0.1 conduction case at a fixed grid
    (12.8 cells/D), reconstructed with WENO5/7, WENO-Z, MUSCL, and MUSCL+THINC interface
    compression. Same dx for all, so any difference in the droop is the scheme's interface
    diffusion -- the direct counterpart to the grid sweep, isolating the advection scheme."""
    series = []
    for name, label, color in RECON_RUNS:
        run = os.path.join(RUNS, name)
        if not os.path.isdir(os.path.join(run, "restart_data")):
            print(f"  recon: {name} not found, skipping")
            continue
        out = color_weighted_vy(run)
        if out is None or len(out[0]) < 5:
            print(f"  recon: {name} not ready, skipping")
            continue
        x, y = v_ygb_ratio(out)
        series.append((x, y, color, label))
    if not series:
        print("  recon: no runs found")
        return
    with plt.rc_context(PLATE_STYLE):
        fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
        ax.plot(SAMAREH_VOF[:, 0], SAMAREH_VOF[:, 1], "s--", color="0.0", ms=5.0, mfc="none", mew=1.3, lw=1.0, label=r"Samareh Fig 5(d), VOF ($Ma=0$)")
        ax.axhline(1.0, color="0.4", lw=1.0, ls="--", zorder=1, label=r"$u_{\mathrm{YGB}}$")
        for x, y, color, label in series:
            ax.plot(x, y, "-", color=color, lw=1.7, alpha=0.95, solid_capstyle="round", label=label)
        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$t / t_r$")
        ax.set_ylabel(r"rise velocity   $u / u_{\mathrm{YGB}}$")
        ax.set_title(r"late-time droop vs reconstruction ($Ma=0.1$, 12.8 cells/$D$, fixed grid)", fontsize=12, loc="left")
        ax.legend(loc="lower left", fontsize=9, frameon=False, ncol=2, columnspacing=1.2, handlelength=1.6)
        sns.despine(ax=ax)
        dst = os.path.join(FIGS, "case1_recon.png")
        fig.savefig(dst, dpi=200)
        plt.close(fig)
        print(f"  wrote {dst}  ({len(series)} runs)")


def cmd_recon(argv):
    os.makedirs(FIGS, exist_ok=True)
    samareh_recon()


COMMANDS = {"samareh": cmd_samareh, "ma": cmd_ma, "fields": cmd_fields, "recon": cmd_recon, "clean": cmd_clean}


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "samareh"  # bare `plot.py` rebuilds the headline overlays
    if cmd not in COMMANDS:
        sys.exit(f"usage: python3 plot.py [{'|'.join(COMMANDS)}] [args]   (default: samareh; see the module docstring)")
    COMMANDS[cmd](sys.argv[2:])


if __name__ == "__main__":
    main()
