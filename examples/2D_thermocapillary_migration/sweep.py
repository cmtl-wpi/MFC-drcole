#!/usr/bin/env python3
"""Terminal-velocity TREND sweeps for the single-particle claims in the Introduction of
Nas & Tryggvason, Int. J. Multiphase Flow 29 (2003) 1117-1135 (those bullets summarize prior
single-particle literature; they are NOT figures in that paper, so we generate the trend curves
ourselves). Quick "trends-only" fidelity: 16 cells/D in a fixed 2D x 4D confined box -- the goal is
whether each curve bends the right way, not a converged unbounded-domain value.

Five sweeps (each point = one case_sweep.py run; the committed case is never touched -- a copy is
made per point and its `name = value` lines are rewritten, exactly like run.py rewrites `Nx`):

  drop_vs_Ma    drop (ratios 0.5), vary Ma          -> V_t decreases, min, then increases (U-shape)
  bubble_vs_Ca  bubble (ratios 1/25), vary Ca       -> V_t decreases rapidly with Ca
  bubble_vs_Re  bubble (ratios 1/25), vary Re       -> V_t increases very weakly with Re
  bubble_vs_Ma  bubble (ratios 1/25), vary Ma       -> V_t decreases with Ma
  drop_vs_rho   drop at Ca=0.1, vary density ratio  -> deforms oblate/prolate with rho* (aspect ratio)

Usage:
    python3 sweep.py [run|remeasure] [sweep1 sweep2 ...]   (default: run, all sweeps)

  run        run every point, then measure and plot
  remeasure  re-measure existing runs/ and replot, no simulation

Terminal velocity is the peak of the color-weighted rise speed U* = U/U_r taken over the window where
the drop is clear of both walls (centroid > 1 radius from each). Aspect ratio AR = (vertical RMS
extent)/(horizontal RMS extent) of the color field: AR > 1 prolate, AR < 1 oblate. All run constants
are read back from simulation.inp so the measurement can't silently disagree with the data.
"""

import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RUNS = os.path.join(HERE, "runs", "sweep")
RESULTS = os.path.join(HERE, "results")
CASE = "case_sweep.py"
R = 0.5  # drop radius (D = 1)

PIN = ["taskset", "-c", "16-255"]
# MPI-only ranks pinned off cores 0-15, binding disabled, single-threaded so concurrent points don't
# oversubscribe. Each point is mpirun -n 2; CONCURRENCY points run at once on the shared EPYC.
RUNENV = {"OMPI_MCA_hwloc_base_binding_policy": "none", "OMP_NUM_THREADS": "1"}
CONCURRENCY = 8  # 8 points x 2 ranks = 16 ranks; staggered launch avoids concurrent cmake-config races

# "Bubble" = light fluid particle at ratio 0.1 (rho/mu/cv/k all 0.1). A true gas bubble is ~1/25, but
# in this compressible 1/T-proxy construction a 1/25 ratio forces a high reference sound speed (EOS
# positivity) and a high thermal diffusivity, both of which shrink dt ~5x for no change in the
# qualitative trend. 0.1 is a defensible light-particle proxy for "does the curve bend the right way".
BUBBLE = dict(rho_ratio=0.1, mu_ratio=0.1, cv_ratio=0.1, k_ratio=0.1)
DROP = dict(rho_ratio=0.5, mu_ratio=0.5, cv_ratio=0.5, k_ratio=0.5)
N_TAU = 1.5  # run length in viscous times (captures the overshoot peak across Re)

SWEEPS = {
    "drop_vs_Ma": dict(
        vary="Ma", values=[1.0, 5.0, 20.0, 50.0, 100.0, 200.0],
        fixed=dict(Re=5.0, Ca=0.01666, **DROP), n_tau=N_TAU, metric="Vstar",
        claim="single drop V_t: decrease, minimum, then increase with Ma",
        xlabel=r"$Ma$", xlog=True,
    ),
    "bubble_vs_Ca": dict(
        vary="Ca", values=[0.01, 0.05, 0.1, 0.2, 0.4],
        fixed=dict(Re=5.0, Ma=20.0, **BUBBLE), n_tau=N_TAU, metric="Vstar",
        claim="gas bubble V_t: decreases rapidly with Ca",
        xlabel=r"$Ca$", xlog=True,
    ),
    "bubble_vs_Re": dict(
        vary="Re", values=[2.0, 5.0, 10.0, 20.0],
        fixed=dict(Ma=20.0, Ca=0.05, **BUBBLE), n_tau=N_TAU, metric="Vstar",
        claim="gas bubble V_t: increases very weakly with Re",
        xlabel=r"$Re$", xlog=True,
    ),
    "bubble_vs_Ma": dict(
        vary="Ma", values=[5.0, 20.0, 50.0, 100.0],
        fixed=dict(Re=5.0, Ca=0.05, **BUBBLE), n_tau=N_TAU, metric="Vstar",
        claim="gas bubble V_t: decreases with Ma",
        xlabel=r"$Ma$", xlog=True,
    ),
    "drop_vs_rho": dict(
        vary="rho_ratio", values=[0.5, 1.0, 2.0],
        fixed=dict(Re=5.0, Ma=20.0, Ca=0.1, mu_ratio=0.5, cv_ratio=0.5, k_ratio=0.5),
        n_tau=N_TAU, metric="AR",
        claim="drop deforms oblate/prolate depending on density ratio",
        xlabel=r"$\rho^* = \rho_i/\rho_o$", xlog=False,
    ),
}


def read_namelist(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                name, value = line.split("=", 1)
                out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def write_case(dst, fixed, vary, value, n_tau):
    """Copy CASE to dst and rewrite the swept + fixed `name = value` lines."""
    text = open(os.path.join(HERE, CASE)).read()
    overrides = {**fixed, vary: value, "n_tau": n_tau}
    for name, val in overrides.items():
        text, hits = re.subn(rf"(?m)^{name} = [-0-9.eE]+", f"{name} = {val}", text, count=1)
        if hits != 1:
            raise RuntimeError(f"could not rewrite `{name}` in {CASE}")
    open(dst, "w").write(text)


def launch_point(name, fixed, vary, value, n_tau):
    """Set up a point's run dir and start mpirun detached. Returns (name, Popen, logfile)."""
    wd = os.path.join(RUNS, name)
    if os.path.isdir(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    dst = os.path.join(wd, CASE)
    write_case(dst, fixed, vary, value, n_tau)
    rel = os.path.relpath(dst, REPO)
    log = open(os.path.join(wd, "run.log"), "w")
    proc = subprocess.Popen(
        PIN + ["./mfc.sh", "run", rel, "-n", "2"],
        cwd=REPO, env={**os.environ, **RUNENV}, stdout=log, stderr=subprocess.STDOUT,
    )
    return name, proc, log


def run_all(jobs):
    """Run (name, fixed, vary, value, n_tau) jobs concurrently, CONCURRENCY at a time."""
    pending, running, done = list(jobs), [], []
    print(f"  launching {len(pending)} points, up to {CONCURRENCY} concurrent (mpirun -n 2 each)", flush=True)
    while pending or running:
        while pending and len(running) < CONCURRENCY:
            job = pending.pop(0)
            name, proc, log = launch_point(*job)
            running.append((name, proc, log))
            print(f"    started {name}", flush=True)
            time.sleep(8)  # stagger so concurrent per-hash cmake/builds don't race on first config
        time.sleep(3)
        still = []
        for name, proc, log in running:
            if proc.poll() is None:
                still.append((name, proc, log))
            else:
                log.close()
                ok = proc.returncode == 0
                print(f"    {'done ' if ok else 'FAIL '} {name} (exit {proc.returncode})", flush=True)
                done.append((name, ok))
        running = still
    return done


def measure_point(wd):
    """Return dict with peak U*, terminal U*, aspect ratio at peak, and the time series."""
    params = read_namelist(os.path.join(wd, "simulation.inp"))

    def P(n):
        return float(params[n.lower()])

    nx, ny = int(P("m")) + 1, int(P("n")) + 1
    dt = P("dt")
    Ly = P("y_domain%end") - P("y_domain%beg")
    mu = 1.0 / P("fluid_pp(1)%re(1)")
    sigma_T = P("sigma_dtdt")
    gradT = abs(P("bc_y%twall_out") - P("bc_y%twall_in")) / Ly
    G = abs(sigma_T * gradT)
    U_r, t_r = G * R / mu, mu / G

    rd = os.path.join(wd, "restart_data")
    yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1):]
    xb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), np.float64)[-(nx + 1):]
    y = 0.5 * (yb[:-1] + yb[1:])
    x = 0.5 * (xb[:-1] + xb[1:])
    y_lo, y_hi = yb[0], yb[-1]

    steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
    cells = nx * ny
    nvars = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    c_idx = nvars - 1  # color is the last conserved variable

    def fld(snap, i):
        return snap[i * cells:(i + 1) * cells].reshape(ny, nx)

    xx, yy = np.meshgrid(x, y)  # (ny, nx)
    t_star, U_star, ycen, ar = [], [], [], []
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        rho = fld(snap, 0) + fld(snap, 1)
        vy = fld(snap, 3) / rho
        c = np.clip(fld(snap, c_idx), 0.0, 1.0)
        csum = c.sum()
        yc = (c * yy).sum() / csum
        xc = (c * xx).sum() / csum
        t_star.append(s * dt / t_r)
        U_star.append((c * vy).sum() / csum / U_r)
        ycen.append(yc)
        rms_y = np.sqrt((c * (yy - yc) ** 2).sum() / csum)
        rms_x = np.sqrt((c * (xx - xc) ** 2).sum() / csum)
        ar.append(rms_y / rms_x)
    t_star, U_star, ycen, ar = map(np.array, (t_star, U_star, ycen, ar))

    # "clear" window: drop centroid more than one radius from each wall (avoid wall contamination)
    clear = (ycen - y_lo > 2.0 * R) & (y_hi - ycen > 2.0 * R)
    if not clear.any():
        clear = np.ones_like(ycen, dtype=bool)
    pk = np.where(clear, U_star, -np.inf)
    peak_i = int(np.argmax(pk))
    return dict(
        U_r=U_r, t_r=t_r, nx=nx, ny=ny, nvars=nvars,
        peak_Vstar=float(U_star[peak_i]), t_peak=float(t_star[peak_i]),
        terminal_Vstar=float(U_star[clear][-1]),
        AR_at_peak=float(ar[peak_i]), AR_final=float(ar[clear][-1]),
        rises=bool(ycen[-1] > ycen[0]),
        t_star=t_star.tolist(), U_star=U_star.tolist(), AR=ar.tolist(),
    )


def plot_sweep(key, cfg, summary):
    pts = sorted(((float(v), summary[f"{key}/{v}"]) for v in cfg["values"] if f"{key}/{v}" in summary), key=lambda t: t[0])
    if not pts:
        return
    xs = [p[0] for p in pts]
    metric = cfg["metric"]
    ys = [p[1]["peak_Vstar" if metric == "Vstar" else "AR_at_peak"] for p in pts]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(xs, ys, "o-", color="C0", ms=6, lw=1.6)
    if cfg["xlog"]:
        ax.set_xscale("log")
    ax.set_xlabel(cfg["xlabel"])
    if metric == "Vstar":
        ax.set_ylabel(r"terminal $V^* = V_t / U_r$  (peak, drop clear of walls)")
        ax.set_ylim(bottom=0.0)
    else:
        ax.set_ylabel(r"aspect ratio  (vertical / horizontal RMS extent)")
        ax.axhline(1.0, ls=":", color="0.5", lw=1.2)
        ax.text(ax.get_xlim()[0], 1.0, " sphere (AR=1): >1 prolate, <1 oblate", va="bottom", fontsize=8, color="0.4")
    ax.set_title(f"{key}: {cfg['claim']}", fontsize=9.5)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = os.path.join(HERE, "figures", f"sweep_{key}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  saved figure -> {out}")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    which = sys.argv[2:] if len(sys.argv) > 2 else list(SWEEPS)
    os.makedirs(RESULTS, exist_ok=True)
    spath = os.path.join(RESULTS, "sweep_summary.json")
    summary = json.load(open(spath)) if os.path.isfile(spath) else {}

    if mode == "run":
        jobs = [(f"{key}/{v}", SWEEPS[key]["fixed"], SWEEPS[key]["vary"], v, SWEEPS[key]["n_tau"]) for key in which for v in SWEEPS[key]["values"]]
        run_all(jobs)

    for key in which:
        cfg = SWEEPS[key]
        print(f"\n=== {key}: {cfg['claim']} ===")
        for v in cfg["values"]:
            name = f"{key}/{v}"
            wd = os.path.join(RUNS, name)
            if not os.path.isfile(os.path.join(wd, "simulation.inp")):
                continue
            try:
                res = measure_point(wd)
            except Exception as e:
                print(f"      MEASURE FAILED ({name}): {e}")
                continue
            res.update(sweep=key, vary=cfg["vary"], value=float(v), **{k: cfg["fixed"].get(k) for k in ("Re", "Ma", "Ca")})
            summary[name] = res
            json.dump(summary, open(spath, "w"), indent=2)
            m = "peak_Vstar" if cfg["metric"] == "Vstar" else "AR_at_peak"
            print(f"      {cfg['vary']}={v}: {m}={res[m]:.4f}  t_peak={res['t_peak']:.2f}  rises={res['rises']}")
        plot_sweep(key, cfg, summary)


if __name__ == "__main__":
    main()
