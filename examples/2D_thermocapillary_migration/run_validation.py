#!/usr/bin/env python3
"""Drive the Samareh-2014 thermocapillary validation sweep and build the comparison figures.

Run groups, all reproducing Samareh's Test Case 1 (zero-Marangoni YGB benchmark):
  grid     : 2D grid convergence, OPEN box, at Samareh's resolutions (64/128/256 cells per box
             width), sigma_T = -0.1 -> converges toward the unbounded-2D analytic 15/16*v_YGB
  ma       : 2D Marangoni-strength sweep at fixed grid (sigma_T = -0.1, -0.05, -0.025) -> U is
             linear in sigma_T, so the ratio stays flat (deficit is grid/window, not strength)
  wall     : Samareh's ACTUAL Test Case 1 geometry (slip walls, drop 1.5D off the cold wall),
             6*tau -> lands on THEIR 2D value ~0.8 at their grid
  3d       : 3D sphere at the coarse grid (12.8 cells/D), quasi-steady window only
  3d_fine   : a finer 3D point (19.2 cells/D); NOT converged -- see `longtime`
  longtime : 6*tau diagnostics (2D + both 3D grids) -> 2D reaches a true plateau, 3D drifts
             unboundedly past v_YGB (frozen-T distortion, no bulk conduction)

Each variant runs a COPY of case.py from its own runs/<name>/ directory (MFC writes output next
to the case file), then measure.py reports the migration velocity (drift-corrected in the open
box, lab-frame in the slip-wall box) and its trend over the final viscous time. Results are
aggregated to results/summary.json and plotted against the three reference values defined below.

Usage (invokes mpirun, so run from a normal shell):
    python3 run_validation.py [grid|ma|wall|3d|3d_fine|longtime|all|remeasure]   (default: all)

`all` = grid + ma + wall + 3d + 3d_fine (the cheap groups; run `longtime` separately).
`remeasure` re-runs measure.py on every existing runs/<name>/ directory and rebuilds
results/summary.json and all figures WITHOUT running any simulation -- use it after editing the
measurement or plotting code.

The analytic density IC is identical across every variant (it depends only on T0 and gradT, which
are fixed), so a single compiled build serves the whole sweep; we still use a full `mfc.sh run`
(not --no-build) each time so an incremental no-op build guards correctness.
"""

import glob
import json
import os
import re
import shutil
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs")
RESULTS = os.path.join(HERE, "results")

# Reference values for v_t/v_YGB (provenance: README "Reference values").
YGB = 1.0  # zero-Ma Stokes SPHERE in an infinite domain (Samareh Eq. 29) -- the normalization
U2D_UNBOUNDED = 15.0 / 16.0  # unbounded 2D CYLINDER analytic, mu*=k*=1: U = |sigma_T|*gradT*r/(8*mu)
SAMAREH_2D = 0.80  # Samareh Fig. 5: 2D in their slip-wall box (3 of 4 methods; the 4th sits ~0.9)
SAMAREH_3D = 0.95  # Samareh Fig. 6: converged 3D VOF in their slip-wall box

# Rank counts respect MFC's decomposition rule: each local block must hold >= num_stcls_min *
# weno_order = 5*5 = 25 cells per split dimension, so small grids cannot over-decompose.
VARIANTS = {
    "grid": [
        dict(name="2D_w064", dim=2, nx=64, dsdt=-0.1, ranks=4),
        dict(name="2D_w128", dim=2, nx=128, dsdt=-0.1, ranks=16),
        dict(name="2D_w256", dim=2, nx=256, dsdt=-0.1, ranks=32),
    ],
    "ma": [
        # 2D_w128 (sigma_T = -0.1) is shared with the grid sweep; add the two weaker points.
        dict(name="2D_w128_dsdt050", dim=2, nx=128, dsdt=-0.05, ranks=16),
        dict(name="2D_w128_dsdt025", dim=2, nx=128, dsdt=-0.025, ranks=16),
    ],
    "wall": [
        dict(name="2D_w064_wall", dim=2, nx=64, dsdt=-0.1, ranks=4, tau=6.0, wall=True),
        dict(name="2D_w128_wall", dim=2, nx=128, dsdt=-0.1, ranks=16, tau=6.0, wall=True),
    ],
    "3d": [
        dict(name="3D_w064", dim=3, nx=64, dsdt=-0.1, ranks=12),
    ],
    "3d_fine": [
        dict(name="3D_w096", dim=3, nx=96, dsdt=-0.1, ranks=45),
    ],
    "longtime": [
        dict(name="2D_w128_t6", dim=2, nx=128, dsdt=-0.1, ranks=16, tau=6.0),
        dict(name="3D_w064_t6", dim=3, nx=64, dsdt=-0.1, ranks=12, tau=6.0),
        dict(name="3D_w096_t6", dim=3, nx=96, dsdt=-0.1, ranks=45, tau=6.0),
    ],
    # Bulk-conduction A/B at 6 tau: does a finite thermal Marangoni number (ma=0.3, the
    # thermal_conduction feature) hold the migration velocity to a plateau where the frozen-T
    # limit (ma=0, no bulk conduction) drifts? The 3D pair is the decisive test -- frozen-T 3D
    # runs away past v_YGB (CONDUCTION_SCOPE.md Part I); conduction should give it a steady state.
    # The 2D pair confirms conduction does not break the already-good 2D plateau.
    "conduction": [
        dict(name="3D_w064_t6_frozen", dim=3, nx=64, dsdt=-0.1, ranks=12, tau=6.0, ma=0.0),
        dict(name="3D_w064_t6_cond", dim=3, nx=64, dsdt=-0.1, ranks=12, tau=6.0, ma=0.3),
        dict(name="2D_w128_t6_frozen", dim=2, nx=128, dsdt=-0.1, ranks=16, tau=6.0, ma=0.0),
        dict(name="2D_w128_t6_cond", dim=2, nx=128, dsdt=-0.1, ranks=16, tau=6.0, ma=0.3),
    ],
}
GROUPS_ALL = ["grid", "ma", "wall", "3d", "3d_fine"]


def measure(workdir):
    """Run measure.py on a case directory and return its RESULT_JSON dict (or None)."""
    meas = subprocess.run([sys.executable, os.path.join(HERE, "measure.py"), workdir], capture_output=True, text=True, check=False)
    print("\n".join(meas.stdout.splitlines()[-4:]))
    if meas.returncode != 0:
        print(f"  MEASURE FAILED: {meas.stderr[-500:]}")
        return None
    for line in meas.stdout.splitlines():
        if line.startswith("RESULT_JSON "):
            return json.loads(line[len("RESULT_JSON ") :])
    print("  no RESULT_JSON in measure output")
    return None


def run_variant(v):
    """Run one MFC case and measure it. Returns the measurement dict (or None on failure)."""
    workdir = os.path.join(RUNS, v["name"])
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    shutil.copy(os.path.join(HERE, "case.py"), os.path.join(workdir, "case.py"))

    env = {
        **os.environ,
        "SAMAREH_DIM": str(v["dim"]),
        "SAMAREH_NX": str(v["nx"]),
        "SAMAREH_DSDT": str(v["dsdt"]),
        "SAMAREH_TAU": str(v.get("tau", 3.0)),
        "SAMAREH_WALL": "1" if v.get("wall") else "0",
        "SAMAREH_MA": str(v.get("ma", 0.3)),  # thermal Marangoni number; 0 = no bulk conduction
    }
    rel_case = os.path.relpath(os.path.join(workdir, "case.py"), REPO)
    print(f"\n>>> {v['name']}: dim={v['dim']} width={v['nx']} dsdt={v['dsdt']} ranks={v['ranks']} tau={v.get('tau', 3.0)} wall={v.get('wall', False)} ma={v.get('ma', 0.3)}", flush=True)
    proc = subprocess.run(["./mfc.sh", "run", rel_case, "-n", str(v["ranks"])], cwd=REPO, env=env, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        print(f"  RUN FAILED (exit {proc.returncode}). Last stderr/stdout:")
        print("\n".join((proc.stdout + proc.stderr).splitlines()[-25:]))
        return None
    return measure(workdir)


def load_summary():
    path = os.path.join(RESULTS, "summary.json")
    return json.load(open(path)) if os.path.isfile(path) else {}


def save_summary(summary):
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(summary, open(os.path.join(RESULTS, "summary.json"), "w"), indent=2)


def ratio_curve(case_dir):
    """Full U(t)/v_YGB history of a run, read straight from the restart data: returns
    (t/tau, ratio) or None if the run is absent. Mirrors measure.py's protocol (lab-frame in the
    slip-wall box, drift-corrected in the open box) -- keep the two in sync."""
    inp = os.path.join(case_dir, "simulation.inp")
    rd = os.path.join(case_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    p = {}
    with open(inp) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                p[k.strip().lower()] = v.strip().rstrip(",")
    f = lambda k: float(p[k.lower()])  # noqa: E731
    nx = int(f("m")) + 1
    dt, mu, dsdt = f("dt"), 1.0 / f("fluid_pp(1)%re(1)"), f("sigma_dtdt")
    wall = int(f("bc_x%beg")) == -2
    tau = 0.25 / mu  # rho*r^2/mu with rho(drop) = 1, r = 0.5
    v_YGB = (2.0 / 15.0) * (-dsdt) * (2.0 / 15.0) * 0.5 / mu  # gradT = 2/15, r = 0.5
    cells = (int(f("m")) + 1) * (int(f("n")) + 1) * (int(f("p")) + 1)
    far = 0.75 * (f("x_domain%end") - f("x_domain%beg")) / 2.0
    xb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
    x = 0.5 * (xb[:-1] + xb[1:])
    is_far = np.abs(x) > far
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    t_tau, ratio = [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        fld = lambda i: snap[i * cells : (i + 1) * cells].reshape(-1, nx)  # noqa: E731
        vx = fld(2) / (fld(0) + fld(1))
        c = np.clip(fld(nvars - 1), 0.0, None)
        u = (c * vx).sum() / c.sum() - (0.0 if wall else vx[:, is_far].mean())
        t_tau.append(s * dt / tau)
        ratio.append(u / v_YGB)
    return np.array(t_tau), np.array(ratio)


def ref_lines(ax, refs):
    """Draw reference levels in neutral styles (data colors are reserved for MFC series)."""
    styles = {
        "ygb": (YGB, ":", "0.45", "YGB sphere = 1 (Samareh Eq. 29)"),
        "u2d": (U2D_UNBOUNDED, "-.", "0.1", "2D unbounded analytic = 15/16"),
        "s2d": (SAMAREH_2D, "--", "0.55", "Samareh 2D $\\approx$ 0.80 (their slip-wall box)"),
        "s3d": (SAMAREH_3D, (0, (5, 2, 1, 2)), "0.55", "Samareh 3D $\\approx$ 0.95 (their slip-wall box)"),
    }
    for key in refs:
        y, ls, color, label = styles[key]
        ax.axhline(y, ls=ls, color=color, lw=1.2, label=label)


def make_figures(summary):
    """Build the four figures from summary.json + the saved restart data."""
    os.makedirs(RESULTS, exist_ok=True)
    stale = os.path.join(RESULTS, "samareh_comparison.png")  # replaced by samareh_geometry.png
    if os.path.isfile(stale):
        os.remove(stale)

    def entries(**conds):
        out = []
        for name, r in summary.items():
            ok = all(abs(r[k] - v) < 1e-9 if isinstance(v, float) else r[k] == v for k, v in conds.items())
            if ok:
                out.append({**r, "name": name})
        return sorted(out, key=lambda r: r["nx_width"])

    def window3(r):  # default-length (3*tau) runs, tolerant of the dt rounding
        return abs(r["t_end_tau"] - 3.0) < 0.2

    # 1) 2D grid convergence (open box) against the unbounded-cylinder analytic.
    g2 = [r for r in entries(dim=2, wall=False, dsigma_dT=-0.1) if window3(r)]
    star = next((r for r in entries(dim=2, wall=False, dsigma_dT=-0.1) if r["t_end_tau"] > 5), None)
    if g2:
        fig, ax = plt.subplots(figsize=(7.0, 4.6))
        ax.plot([r["cells_per_D"] for r in g2], [r["ratio_window"] for r in g2], "o-", color="C0", label=r"MFC 2D open box (3$\tau$-window mean)")
        if star:
            ax.plot(star["cells_per_D"], star["ratio_window"], "*", color="C0", ms=15, label=rf"6$\tau$ run, true plateau ({star['ratio_window']:.2f})")
        ref_lines(ax, ["ygb", "u2d", "s2d"])
        ax.set_xlabel(r"cells per droplet diameter $D/\Delta x$")
        ax.set_ylabel(r"$v_t / v_{\mathrm{YGB}}$")
        ax.set_title(r"2D grid convergence, open box ($\sigma_T = -0.1$)")
        ax.set_ylim(0.6, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS, "grid_convergence.png"), dpi=150)
        plt.close(fig)

    # 2) The geometry experiment: open box vs Samareh's slip-wall box, as U(t) curves.
    curves = [
        ("2D_w128_t6", "open box, centered drop (25.6 cells/$D$, drift-corrected)", "C0"),
        ("2D_w128_wall", "Samareh box: slip walls, drop 1.5$D$ off wall (25.6 cells/$D$)", "C3"),
        ("2D_w064_wall", "Samareh box at their headline grid (12.8 cells/$D$)", "C1"),
    ]
    data = [(lab, col, ratio_curve(os.path.join(RUNS, name))) for name, lab, col in curves]
    data = [d for d in data if d[2] is not None]
    if data:
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        for lab, col, (t, rr) in data:
            ax.plot(t, rr, "-", color=col, lw=1.6, label=lab)
        ref_lines(ax, ["ygb", "u2d", "s2d"])
        ax.set_xlabel(r"$t/\tau$  ($\tau = \rho r^2/\mu$, viscous time)")
        ax.set_ylabel(r"$U / v_{\mathrm{YGB}}$")
        ax.set_title("Boundary conditions select the target: open box vs Samareh's box")
        ax.set_xlim(0.0, 6.05)
        ax.set_ylim(0.0, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS, "samareh_geometry.png"), dpi=150)
        plt.close(fig)

    # 3) Marangoni-strength sweep at fixed grid (width 128, 2D, open): U scales LINEARLY with
    # |sigma_T| (YGB is linear in the slope), so the ratio v_t/v_YGB stays flat.
    ma = [r for r in entries(dim=2, wall=False, nx_width=128) if window3(r)]
    ma.sort(key=lambda r: abs(r["dsigma_dT"]))
    if len(ma) >= 2:
        strengths = np.array([abs(r["dsigma_dT"]) for r in ma])
        U = np.array([r["U_window"] for r in ma])
        vy = np.array([r["v_YGB"] for r in ma])
        smax = strengths.max() * 1.08
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.5, 4.3))
        slope = float(np.mean(U / strengths))  # U/|sigma_T| ~ constant if linear
        axL.plot([0, smax], [0, (vy / strengths).mean() * smax], "--", color="0.5", lw=1.2, label=r"YGB: $v_{\mathrm{YGB}} \propto |\sigma_T|$")
        axL.plot([0, smax], [0, slope * smax], ":", color="C2", lw=1.2, label="fit through origin")
        axL.plot(strengths, U, "s", color="C2", ms=8, label=r"MFC 2D (3$\tau$-window mean)")
        axL.set_xlabel(r"Marangoni strength $|\sigma_T|$  ($|\nabla T|$ fixed)")
        axL.set_ylabel(r"migration velocity $U$")
        axL.set_xlim(0.0, smax)
        axL.set_ylim(bottom=0.0)
        axL.set_title("Linearity of the $\\sigma(T)$ coupling")
        axL.grid(alpha=0.3)
        axL.legend(loc="upper left", fontsize=9)
        axR.plot(strengths, U / vy, "s-", color="C2", label=r"MFC 2D (3$\tau$-window mean)")
        ref_lines(axR, ["ygb", "u2d"])
        axR.annotate(r"Samareh's $\sigma_T = -0.1$", xy=(0.1, (U / vy)[-1]), xytext=(0.055, 0.70), fontsize=9, color="0.3", arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0))
        axR.set_xlabel(r"Marangoni strength $|\sigma_T|$")
        axR.set_ylabel(r"$v_t / v_{\mathrm{YGB}}$")
        axR.set_xlim(0.0, smax)
        axR.set_ylim(0.6, 1.05)
        axR.set_title("Ratio is strength-independent (flat)")
        axR.grid(alpha=0.3)
        axR.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS, "marangoni_sweep.png"), dpi=150)
        plt.close(fig)

    # 4) Long-time behaviour: 2D plateaus, 3D drifts unboundedly (frozen-T runaway).
    curves = [
        ("2D_w128_t6", "2D, 25.6 cells/$D$", "C0", "o"),
        ("3D_w064_t6", "3D, 12.8 cells/$D$", "C1", "s"),
        ("3D_w096_t6", "3D, 19.2 cells/$D$", "C2", "^"),
    ]
    data = [(lab, col, mk, ratio_curve(os.path.join(RUNS, name))) for name, lab, col, mk in curves]
    data = [d for d in data if d[3] is not None]
    if data:
        fig, ax = plt.subplots(figsize=(7.6, 4.8))
        ax.axvspan(0.0, 3.0, color="0.93", zorder=0)
        ax.text(1.5, 0.06, "3$\\tau$ measurement window\n(quasi-steady ratios)", ha="center", fontsize=8, color="0.35")
        for lab, col, mk, (t, rr) in data:
            ax.plot(t, rr, mk + "-", color=col, ms=4, lw=1.4, label=lab)
        ref_lines(ax, ["ygb", "s3d"])
        if any(lab.startswith("3D, 19.2") for lab, *_ in data):
            t, rr = next(d[3] for d in data if d[0].startswith("3D, 19.2"))
            ax.annotate(
                "frozen-$T$ runaway (no bulk conduction):\nnot a physical terminal velocity",
                xy=(t[-1], rr[-1]),
                xytext=(3.1, 1.22),
                fontsize=8,
                color="C2",
                arrowprops=dict(arrowstyle="->", color="C2", lw=1.0),
            )
        ax.set_xlabel(r"$t/\tau$  ($\tau = \rho r^2/\mu$, viscous time)")
        ax.set_ylabel(r"$U / v_{\mathrm{YGB}}$")
        ax.set_title("Long-time behaviour: 2D plateaus, 3D drifts past YGB")
        ax.set_xlim(left=0.0)
        ax.grid(alpha=0.3)
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS, "longtime_drift.png"), dpi=150)
        plt.close(fig)


def print_table(summary):
    print("\n" + "=" * 100)
    print(f"{'variant':>18} {'dim':>3} {'cells/D':>8} {'mode':>5} {'sigma_T':>8} {'t_end/tau':>9} {'U/vYGB win':>10} {'final':>7} {'slope/tau':>9}")
    for name, r in sorted(summary.items(), key=lambda kv: (kv[1]["dim"], kv[1]["wall"], kv[1]["nx_width"], kv[1]["t_end_tau"], abs(kv[1]["dsigma_dT"]))):
        mode = "wall" if r["wall"] else "open"
        print(
            f"{name:>18} {r['dim']:>3} {r['cells_per_D']:>8.1f} {mode:>5} {r['dsigma_dT']:>8.3f} {r['t_end_tau']:>9.1f} {r['ratio_window']:>+10.3f} {r['ratio_final']:>+7.3f} {r['slope_per_tau']:>+9.3f}"
        )
    print("=" * 100)


if __name__ == "__main__":
    group = sys.argv[1] if len(sys.argv) > 1 else "all"

    if group == "remeasure":
        # Rebuild summary + figures from whatever run directories exist; no simulations.
        summary = {}
        for name in sorted(os.listdir(RUNS)) if os.path.isdir(RUNS) else []:
            workdir = os.path.join(RUNS, name)
            if not os.path.isfile(os.path.join(workdir, "simulation.inp")):
                continue
            print(f"\n>>> remeasure {name}", flush=True)
            res = measure(workdir)
            if res is not None:
                summary[name] = res
        save_summary(summary)
    else:
        todo = [v for g in (GROUPS_ALL if group == "all" else [group]) for v in VARIANTS[g]]
        summary = load_summary()
        for v in todo:
            res = run_variant(v)
            if res is not None:
                summary[v["name"]] = res
                save_summary(summary)  # checkpoint after every successful variant

    make_figures(summary)
    print_table(summary)
    print(f"\nwrote {os.path.join(RESULTS, 'summary.json')} and figures in {RESULTS}/")
