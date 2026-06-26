#!/usr/bin/env python3
"""Comparison still + tiled animation for the six case1_fig5 runs.

fig5 itself compares them *quantitatively* (rise velocity vs t/t_r). This is the
*visual* companion: the drop interface field (color_function, last conserved var; its
0.5 level is the interface "isosurface") for all six runs, laid out as rows = physics
(frozen-T Ma=0 / conduction Ma=0.1), cols = grid (12.8 / 25.6 / 51.2 cells/D). Same
physical box, same 0..1 color scale, so resolution and frozen-vs-conduction differences
read directly. Constants come from each run's simulation.inp, so the panels can't
silently disagree with the data.

    python3 compare_fig5_fields.py            # still montage at matched t/t_r
                                              #   -> figures/case1_fig5_fields.png
    python3 compare_fig5_fields.py anim       # tiled animation, all six in lockstep on t/t_r
                                              #   -> animations/case1_fig5_tiled.mp4

mfc viz is single-run, so a six-run tiled view is a custom plot (same pattern as plot.py
and compare_tc3_visc.py). The mp4 is encoded with the same imageio/libx264/yuv420p
settings mfc viz uses, so it plays anywhere the per-run clips do.

Conserved layout (model_eqns=3, num_fluids=2): color c is the last variable.
"""

import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ANIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "animations")

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
FIGS = os.path.join(HERE, "figures")
R = 0.5  # drop radius (D = 1)
GRADT = 2.0 / 15.0  # imposed |dT/dy|, common to TC1/TC2

# (run dir, row, col, case label) — the exact six runs samareh_fig5() plots.
PANELS = [
    ("tc1/frozen/w064", 0, 0, "frozen-T  ($Ma=0$)"),
    ("tc1/frozen/w128", 0, 1, "frozen-T  ($Ma=0$)"),
    ("tc1/frozen/w256", 0, 2, "frozen-T  ($Ma=0$)"),
    ("tc1/ma0p1/w064/sc050", 1, 0, "conduction  ($Ma=0.1$)"),
    ("tc1/ma0p1/w128/sc050", 1, 1, "conduction  ($Ma=0.1$)"),
    ("tc1/ma0p1/w256/sc050", 1, 2, "conduction  ($Ma=0.1$)"),
]


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def load_run(run_dir):
    """Grid, snapshot steps, t_r and dt for one run. Returns a dict or None if not ready."""
    inp = os.path.join(run_dir, "simulation.inp")
    rd = os.path.join(run_dir, "restart_data")
    if not (os.path.isfile(inp) and os.path.isdir(rd)):
        return None
    P = read_namelist(inp)
    f = lambda k: float(P[k.lower()])  # noqa: E731
    nx, ny = int(f("m")) + 1, int(f("n")) + 1
    steps = sorted(int(m.group(1)) for ff in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", ff)))
    if not steps:
        return None
    xb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), np.float64)[-(nx + 1) :]
    yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
    mu = 1.0 / f("fluid_pp(1)%re(1)")
    t_r = mu / abs(f("sigma_dtdt") * GRADT)
    return dict(rd=rd, nx=nx, ny=ny, ncell=nx * ny, dt=f("dt"), t_r=t_r, xb=xb, yb=yb, steps=steps, cells_per_D=(nx) / 5.0)


def color_at(run, step):
    """color_function on the (ny, nx) grid for one snapshot (last conserved variable)."""
    snap = np.fromfile(os.path.join(run["rd"], f"lustre_{step}.dat"), np.float64)
    nvars = snap.size // run["ncell"]
    c = snap[(nvars - 1) * run["ncell"] :].reshape(run["ny"], run["nx"])
    return np.clip(c, 0.0, 1.0)


def main():
    runs = {}
    for name, _, _, _ in PANELS:
        r = load_run(os.path.join(RUNS, name))
        if r is None:
            sys.exit(f"run not ready (no snapshots): {name}")
        runs[name] = r

    # Match on t/t_r: use the largest time every run actually reaches, so no panel extrapolates.
    target_ttr = min((r["steps"][-1] * r["dt"]) / r["t_r"] for r in runs.values())

    # First pass: pick each run's matched snapshot, load its color field, and find the drop
    # centroid. The drops are tiny (D=1) in the 5D box, so we zoom to a shared window that
    # still spans every drop — this keeps the (small) rise-height differences visible while
    # making the interface shape/sharpness legible.
    snaps = {}
    yc_all = []
    for name, _, _, _ in PANELS:
        r = runs[name]
        step = min(r["steps"], key=lambda s: abs(s * r["dt"] / r["t_r"] - target_ttr))
        c = color_at(r, step)
        xcc = 0.5 * (r["xb"][:-1] + r["xb"][1:])
        ycc = 0.5 * (r["yb"][:-1] + r["yb"][1:])
        cy = (c * ycc[:, None]).sum() / c.sum()
        snaps[name] = dict(step=step, ttr=step * r["dt"] / r["t_r"], c=c, xcc=xcc, ycc=ycc, cy=cy)
        yc_all.append(cy)

    # Shared crop: x centred on the (centred) drop, y spanning all centroids + a 1.25D margin.
    half_x = 1.6
    y_lo = min(yc_all) - 1.25
    y_hi = max(yc_all) + 1.25

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.8), constrained_layout=True, sharex=True, sharey=True)
    mesh = None
    for name, row, col, case_lbl in PANELS:
        r, s = runs[name], snaps[name]
        ax = axes[row, col]
        mesh = ax.pcolormesh(s["xcc"], s["ycc"], s["c"], cmap="RdBu_r", vmin=0.0, vmax=1.0, shading="auto")
        ax.contour(s["xcc"], s["ycc"], s["c"], levels=[0.5], colors="k", linewidths=1.1)
        ax.plot(0.0, s["cy"], "k+", ms=8, mew=1.4)
        ax.set_aspect("equal")
        ax.set_xlim(-half_x, half_x)
        ax.set_ylim(y_lo, y_hi)
        ax.set_title(rf"{r['cells_per_D']:.1f} cells/$D$   ($t/t_r={s['ttr']:.2f}$)", fontsize=10)
        if col == 0:
            ax.set_ylabel(f"{case_lbl}\n$y$", fontsize=10)
        if row == 1:
            ax.set_xlabel("$x$", fontsize=10)

    cbar = fig.colorbar(mesh, ax=axes, shrink=0.85, pad=0.02)
    cbar.set_label("color function  (0 = ambient, 1 = drop;  0.5 = interface)")
    fig.suptitle(rf"2D thermocapillary rise — interface field, six fig5 runs at $t/t_r \approx {target_ttr:.2f}$", fontsize=13)

    os.makedirs(FIGS, exist_ok=True)
    out = os.path.join(FIGS, "case1_fig5_fields.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"matched t/t_r = {target_ttr:.3f}")
    print(f"wrote {out}")


def animate(n_frames=121, fps=12):
    """Tiled animation: the same 2x3 montage, advanced in lockstep on a shared t/t_r timeline.

    Each run picks its nearest snapshot to the current t/t_r (coarse runs stair-step, which is
    honest — they simply have fewer frames). The lab-frame crop is fixed for the whole clip and
    sized to contain every drop at every time, so the (small) rise is visible as real motion.
    """
    import imageio

    runs = {}
    for name, _, _, _ in PANELS:
        r = load_run(os.path.join(RUNS, name))
        if r is None:
            sys.exit(f"run not ready (no snapshots): {name}")
        r["ttr_steps"] = np.array(r["steps"]) * r["dt"] / r["t_r"]
        r["xcc"] = 0.5 * (r["xb"][:-1] + r["xb"][1:])
        r["ycc"] = 0.5 * (r["yb"][:-1] + r["yb"][1:])
        runs[name] = r

    target_max = min(r["ttr_steps"][-1] for r in runs.values())
    timeline = np.linspace(0.0, target_max, n_frames)

    # Pre-pass: the step each (run, frame) lands on, plus a color/centroid cache (dedups the
    # repeated steps coarse runs land on) and the global centroid range for the fixed crop.
    cache = {}  # (name, step) -> (color float32, centroid_y)
    sel = {name: [] for name in runs}
    yc_min, yc_max = np.inf, -np.inf
    for name, r in runs.items():
        for tt in timeline:
            step = r["steps"][int(np.argmin(np.abs(r["ttr_steps"] - tt)))]
            if (name, step) not in cache:
                c = color_at(r, step).astype(np.float32)
                cy = float((c * r["ycc"][:, None]).sum() / c.sum())
                cache[(name, step)] = (c, cy)
            sel[name].append(step)
            yc_min = min(yc_min, cache[(name, step)][1])
            yc_max = max(yc_max, cache[(name, step)][1])

    half_x = 1.6
    y_lo, y_hi = yc_min - 1.25, yc_max + 1.25

    # One reusable figure (fixed size -> identical even frame dims for libx264/yuv420p).
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.8), dpi=100, constrained_layout=True, sharex=True, sharey=True)
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(0.0, 1.0))
    cbar = fig.colorbar(sm, ax=axes, shrink=0.85, pad=0.02)
    cbar.set_label("color function  (0 = ambient, 1 = drop;  0.5 = interface)")

    os.makedirs(ANIM, exist_ok=True)
    out = os.path.join(ANIM, "case1_fig5_tiled.mp4")
    with imageio.get_writer(out, fps=fps, codec="libx264", pixelformat="yuv420p", macro_block_size=1, ffmpeg_log_level="error") as writer:
        for fi, tt in enumerate(timeline):
            for name, row, col, case_lbl in PANELS:
                r = runs[name]
                step = sel[name][fi]
                c, cy = cache[(name, step)]
                ax = axes[row, col]
                ax.cla()
                ax.pcolormesh(r["xcc"], r["ycc"], c, cmap="RdBu_r", vmin=0.0, vmax=1.0, shading="auto")
                ax.contour(r["xcc"], r["ycc"], c, levels=[0.5], colors="k", linewidths=1.1)
                ax.plot(0.0, cy, "k+", ms=8, mew=1.4)
                ax.set_aspect("equal")
                ax.set_xlim(-half_x, half_x)
                ax.set_ylim(y_lo, y_hi)
                ax.set_title(rf"{r['cells_per_D']:.1f} cells/$D$", fontsize=10)
                if col == 0:
                    ax.set_ylabel(f"{case_lbl}\n$y$", fontsize=10)
                if row == 1:
                    ax.set_xlabel("$x$", fontsize=10)
            fig.suptitle(rf"2D thermocapillary rise — interface field, six fig5 runs   $t/t_r = {tt:5.2f}$", fontsize=13)
            fig.canvas.draw()
            frame = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
            writer.append_data(np.ascontiguousarray(frame))
            if fi % 20 == 0 or fi == n_frames - 1:
                print(f"  frame {fi + 1}/{n_frames}  (t/t_r={tt:.2f})")
    plt.close(fig)
    print(f"matched span t/t_r = 0..{target_max:.2f}  ({n_frames} frames @ {fps} fps = {n_frames / fps:.1f}s)")
    print(f"wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "anim":
        animate()
    else:
        main()
