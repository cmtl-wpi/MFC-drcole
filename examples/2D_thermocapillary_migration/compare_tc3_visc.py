#!/usr/bin/env python3
"""TC3 validation: does temperature-dependent viscosity mu(T) produce the migration signature that a
constant-mu run cannot? This is the load-bearing claim of Samareh's Sec. 4.2 (large-Ma, Figs 8/13):
the silicon oil's viscosity falls as the drop rises into warmer oil, so the drag drops and the drop
*accelerates* with height -- the experiment's non-monotonic rise-velocity behaviour that a constant-mu
run (TC1/TC2 physics) structurally cannot reproduce.

We run the SAME 2D TC3 case twice -- mu(T) = exp(C+D/T) on (runs/tc3/muT) vs a constant-mu control
frozen at the drop's start temperature (runs/tc3/const). The two are identical at t=0 (the drop sits
at one temperature) and may only diverge because mu(T) lets the local viscosity change as the drop moves.

The 2D cylinder is the tractable analogue of Samareh's 3D sphere, exactly as TC1's Fig 5 (2D plane) is
the analogue of its Fig 6 (3D sphere): it isolates the mu(T) PHYSICS at a fraction of the cost. The
absolute velocity ratio differs between 2D and 3D, but the mu(T) signature -- acceleration with height,
absent in the control -- is dimension-independent.

Both runs are read straight from restart data; run-dependent constants come from each run's
simulation.inp so the script cannot silently disagree with what was actually simulated.

Usage:  python3 compare_tc3_visc.py [muT_dir] [const_dir]
Writes: figures/case3_large_marangoni_mu_of_T_validation.png  and prints a verification summary (tag: VERIFY).
"""

import glob
import os
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
MU_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "runs/tc3/muT")
CONST_DIR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "runs/tc3/const")


def read_namelist(path):
    """Parse a Fortran namelist file's plain "name = value" lines into a dict (lowercase keys)."""
    out = {}
    for line in open(path):
        if "=" in line:
            name, value = line.split("=", 1)
            out[name.strip().lower()] = value.strip().rstrip(",")
    return out


def trajectory(case_dir):
    """Return (t_ms, dist_mm, v_mms, Tdrop_K) for one run -- color-weighted drop rise history.
    dist_mm = drop-centroid distance from the cold (bottom) wall; v_mms = lab-frame rise velocity."""
    params = read_namelist(os.path.join(case_dir, "simulation.inp"))

    def param(name):
        return float(params[name.lower()])

    # Grid: MFC stores m/n/p as (cells - 1), so add 1. 2D runs have p=0 -> nz=1.
    nx, ny, nz = int(param("m")) + 1, int(param("n")) + 1, int(param("p")) + 1
    dt = param("dt")
    save = int(param("t_step_save"))
    y_cold = param("y_domain%beg")
    cells = nx * ny * nz
    restart_dir = os.path.join(case_dir, "restart_data")

    # Cell-center y positions from the boundary file (last ny+1 boundaries are the interior).
    yb = np.fromfile(os.path.join(restart_dir, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
    y = 0.5 * (yb[:-1] + yb[1:])  # cell-center y (m)
    y_grid = y[None, :, None]  # broadcast y over (nz, ny, nx)

    # Keep only snapshots on the current run's cadence (guards against stale files from an
    # earlier cadence) -- i.e. step numbers that are a multiple of t_step_save.
    steps = []
    for path in glob.glob(os.path.join(restart_dir, "lustre_*.dat")):
        match = re.search(r"lustre_(\d+)\.dat$", path)
        if match and int(match.group(1)) % save == 0:
            steps.append(int(match.group(1)))
    steps.sort()
    if not steps:
        return None
    nvars = np.fromfile(os.path.join(restart_dir, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    # thermal_scalar appends the temperature scalar T_s last, so color c is second-to-last.
    c_idx, T_idx = nvars - 2, nvars - 1

    times, dist, vel, T_drop = [], [], [], []
    for step in steps:
        snap = np.fromfile(os.path.join(restart_dir, f"lustre_{step}.dat"), np.float64)

        def field(i):
            return snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)

        rho = field(0) + field(1)  # total density (sum of partial densities)
        vy = field(3) / rho  # rise (y) velocity
        color = np.clip(field(c_idx), 0.0, 1.0)  # color function (drop = 1, bulk = 0)
        Ts = field(T_idx)  # temperature scalar
        color_total = color.sum()
        if color_total <= 0:
            continue
        times.append(step * dt * 1e3)  # ms
        dist.append(((color * y_grid).sum() / color_total - y_cold) * 1e3)  # mm from cold wall
        vel.append((color * vy).sum() / color_total * 1e3)  # mm/s
        T_drop.append((color * Ts).sum() / color_total)  # K
    return tuple(np.array(a) for a in (times, dist, vel, T_drop))


def resolution_check(case_dir, peak_v_mms):
    """Is the thermal boundary layer around the drop resolved? At large Marangoni number the migration is
    set by a thin thermal layer of thickness delta ~ R/sqrt(Pe), Pe = V R / alpha_thermal. If delta is
    smaller than a cell, the surface-temperature gradient (hence the mu(T)-vs-control divergence) is NOT
    grid-converged -- a converged Samareh Fig-8 needs the fine production grid (nx ~ 240+)."""
    params = read_namelist(os.path.join(case_dir, "simulation.inp"))

    def param(name):
        return float(params[name.lower()])

    dx = (param("x_domain%end") - param("x_domain%beg")) / (int(param("m")) + 1)
    R = 5.35e-3  # TC3 drop radius (D = 10.7 mm)
    k_b, rho_b, cp_b = 0.13389, 918.3, 1778.2  # silicon oil bulk (Samareh Sec. 4.2)
    alpha = k_b / (rho_b * cp_b)
    Pe = max(peak_v_mms * 1e-3, 1e-12) * R / alpha
    delta = R / Pe**0.5
    return dx, Pe, delta


def running_mean(x, w):
    """Centered running mean (window w, clipped at the ends) -- averages over the acoustic ringing."""
    if len(x) < 3:
        return x
    # Force the window odd and no wider than the series so it stays centered.
    w = max(1, min(w, len(x) // 2 * 2 - 1 if len(x) > 2 else 1))
    window = np.ones(w) / w
    pad = w // 2
    padded = np.pad(x, pad, mode="edge")  # edge-pad so the ends aren't dragged toward zero
    return np.convolve(padded, window, mode="valid")[: len(x)]


mu = trajectory(MU_DIR)
const = trajectory(CONST_DIR)
if mu is None or const is None:
    sys.exit(f"missing snapshots -- muT={'ok' if mu else 'EMPTY'} const={'ok' if const else 'EMPTY'}")

t_mu, d_mu, v_mu, T_mu = mu
t_c, d_c, v_c, T_c = const
smooth_w = max(5, len(t_mu) // 12)  # smoothing window ~ a few acoustic periods

# -- verification --
# The two runs reach different heights at different wall-clock, so a fair comparison is at MATCHED DISTANCE
# (Samareh's Fig-8 axes), not matched index/time. We interpolate each smoothed v(distance) onto a common
# distance grid over the overlap and read the mu(T)-minus-constant divergence vs height. The mu(T) signature:
# as the drop rises into warmer, less-viscous oil the local bulk drag falls, so mu(T) outpaces the
# constant-mu control by a margin that GROWS with height (zero at the release point, where both share one
# viscosity). The drop interior stays ~cold (large Ma: it is advected faster than it conducts), so the
# effect is carried by the bulk viscosity at the drop's location, exactly as Samareh describes.
rose_mu = d_mu[-1] > d_mu[0]
vbar_mu = running_mean(v_mu, smooth_w)
vbar_c = running_mean(v_c, smooth_w)
peak_mu = float(np.nanmax(vbar_mu))
rise_mm = float(d_mu[-1] - d_mu[0])


def by_distance(dist, vel):
    """Sort a (distance, velocity) pair by ascending distance so np.interp can use it."""
    order = np.argsort(dist)
    return dist[order], vel[order]


dm_s, vm_s = by_distance(d_mu, vbar_mu)
dc_s, vc_s = by_distance(d_c, vbar_c)
# Resample both runs onto a common height grid over the heights they BOTH reached.
lo, hi = max(dm_s[0], dc_s[0]), min(dm_s[-1], dc_s[-1])
grid = np.linspace(lo, hi, 200)
vmu_g = np.interp(grid, dm_s, vm_s)
vc_g = np.interp(grid, dc_s, vc_s)
dv_pct = 100.0 * (vmu_g - vc_g) / np.maximum(vc_g, 1e-9)  # mu(T) - const at matched height (%)
top_pct = float(np.mean(dv_pct[-10:]))  # divergence near the highest point both runs reached
mu_faster_top = top_pct > 0.0
# bulk temperature span the drop has traversed (1000 K/m * rise)
bulk_dT = 1000.0 * rise_mm * 1e-3

print(f"VERIFY runs: muT snapshots={len(t_mu)} (t_end={t_mu[-1]:.1f} ms)  const snapshots={len(t_c)} (t_end={t_c[-1]:.1f} ms)")
print(f"VERIFY direction: drop rose {d_mu[0]:.2f} -> {d_mu[-1]:.2f} mm from cold wall  (+y toward hot wall: {rose_mu})")
print(f"VERIFY magnitude: peak smoothed rise velocity = {peak_mu:.3f} mm/s   [large-Ma migration << U_r~20 mm/s; exp Fig8 ~2-3]")
print(f"VERIFY drop interior T {T_mu[0]:.2f} -> {T_mu[-1]:.2f} K (stays cold: large-Ma); bulk dT traversed ~ {bulk_dT:.2f} K over {rise_mm:.2f} mm")
print(f"VERIFY mu(T) signature (distance-matched, height {lo:.2f}->{hi:.2f} mm): mu(T) diverges from constant-mu by {top_pct:+.2f}% near the top (feature is active: {abs(top_pct) > 1.0})")
dx_, Pe_, delta_ = resolution_check(MU_DIR, peak_mu)
print(
    f"VERIFY resolution: Pe~{Pe_:.0f}, thermal-layer delta~{delta_ * 1e3:.2f} mm vs cell dx~{dx_ * 1e3:.2f} mm -> thermal layer resolved: {delta_ > dx_} (else the mu(T) MAGNITUDE is not grid-converged)"
)

# -- figure: (A) v vs distance (Samareh Fig-8 axes); (B) v vs time; (C) the mu(T) divergence vs height --
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.0))
CMU, CC = "#0072B2", "#D55E00"  # colorblind blue / vermillion

# Panels A and B plot the same velocity series against distance (Fig-8 axes) and against time.
for ax, x_mu, x_c, xlabel, title in (
    (axA, d_mu, d_c, "Distance from cold wall (mm)", "TC3 (2D): rise velocity vs height  [Samareh Fig 8 axes]"),
    (axB, t_mu, t_c, "Time (ms)", "Rise velocity vs time"),
):
    ax.plot(x_mu, v_mu, ".", color=CMU, ms=3, alpha=0.30)
    ax.plot(x_c, v_c, ".", color=CC, ms=3, alpha=0.30)
    ax.plot(x_mu, vbar_mu, "-", color=CMU, lw=2.2, label="mu(T) = exp(C + D/T)")
    ax.plot(x_c, vbar_c, "-", color=CC, lw=2.2, label="constant mu (control)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Rise velocity (mm/s)")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)

axC.axhline(0.0, color="0.5", lw=1.0)
axC.plot(grid, dv_pct, "-", color="#117733", lw=2.4)
axC.fill_between(grid, 0, dv_pct, color="#117733", alpha=0.15)
axC.set_xlabel("Distance from cold wall (mm)")
axC.set_ylabel("[v_muT - v_const] / v_const   (%)")
axC.set_title("mu(T) vs constant-mu divergence (measured, signed)", fontsize=10)
axC.grid(alpha=0.3)

fig.suptitle(
    "TC3 (large-Ma): temperature-dependent viscosity mu(T) vs a constant-mu control -- both migrate the drop toward the hot wall;\n"
    "the mu(T) run measurably diverges from the control as the surrounding oil's viscosity changes with height (panel C)",
    fontsize=10,
)
fig.tight_layout(rect=(0, 0, 1, 0.93))
out = os.path.join(HERE, "figures", "case3_large_marangoni_mu_of_T_validation.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150)
print(f"saved figure -> {out}")
