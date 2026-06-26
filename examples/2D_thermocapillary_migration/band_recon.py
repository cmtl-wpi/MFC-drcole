#!/usr/bin/env python3
"""Interface-band diagnostic + overlay figure for the MUSCL vs MUSCL+THINC (int_comp) TC2 sweep.

Answers the THINC experiment's three questions quantitatively, for each variant in
runs/recon_tc2/{muscl,muscl_thinc}:

  (a) does int_comp=T keep the VOLUME-FRACTION band (drop alpha, the field THINC compresses)
      sharper over time?
  (b) does it change the migration velocity U*(t*)?
  (c) crucially -- does the COLOR-function band (which THINC never touches, and which drives the
      surface-tension CSF force) stay diffused even if alpha sharpens?

Band-thickness recipe (the droop study's co-area estimator), per snapshot, per field f in [0,1]:
    band_area = count(0.1 < f < 0.9) * dx*dy          (mixed-cell area)
    length    = sum(|grad f|) * dx*dy                 (interface length; co-area: int |grad f| dA)
    thickness = band_area / length / dx               (in CELL widths, dx)
Applied to BOTH the color function (var nvars-1) AND the drop volume fraction (the adv field
that goes 0->1; auto-detected as the [0,1]-bounded non-color field most correlated with color at
t=0). U*(t*) reuses measure.py's color-weighted rise velocity / U_r.

Usage:  python3 band_recon.py [out_png]
Default out_png = runs/recon_tc2/recon_thinc_compare.png. Prints a quantitative summary.
"""

import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNROOT = os.path.join(HERE, "runs", "recon_tc2")
VARIANTS = ["muscl", "muscl_thinc"]
LABEL = {"muscl": "MUSCL (int_comp=F)", "muscl_thinc": "MUSCL+THINC (int_comp=T)"}
COLOR = {"muscl": "C0", "muscl_thinc": "C3"}


def read_namelist(path):
    out = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().lower()] = v.strip().rstrip(",")
    return out


def reference_scales(P, Ly, R):
    """Marangoni U_r, t_r from the run's constants (mirrors measure.py fig7)."""
    mu = 1.0 / float(P["fluid_pp(1)%re(1)"])
    sigma_T = float(P["sigma_dtdt"])
    gradT = abs(float(P["bc_y%twall_out"]) - float(P["bc_y%twall_in"])) / Ly
    G = abs(sigma_T * gradT)
    return mu, sigma_T, gradT, G, G * R / mu, mu / G  # mu, sigma_T, gradT, G, U_r, t_r


def load_variant(variant, R=0.5):
    """Walk a variant's restart_data; return time series of t*, U*, and band thickness (cells)
    for the color function and the drop volume fraction."""
    d = os.path.join(RUNROOT, variant)
    P = read_namelist(os.path.join(d, "simulation.inp"))
    nx, ny, nz = int(P["m"]) + 1, int(P["n"]) + 1, int(P["p"]) + 1
    cells = nx * ny * nz
    dt = float(P["dt"])
    Wx = float(P["x_domain%end"]) - float(P["x_domain%beg"])
    Ly = float(P["y_domain%end"]) - float(P["y_domain%beg"])
    dx, dy = Wx / nx, Ly / ny
    mu, sigma_T, gradT, G, U_r, t_r = reference_scales(P, Ly, R)

    rd = os.path.join(d, "restart_data")
    yb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), np.float64)[-(ny + 1) :]
    y = 0.5 * (yb[:-1] + yb[1:])
    steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (m := re.search(r"lustre_(\d+)\.dat$", f)))
    if not steps:
        sys.exit(f"no snapshots in {rd}")

    snap0 = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64)
    nvars = snap0.size // cells
    c_idx = nvars - 1

    def field(snap, i):
        return snap[i * cells : (i + 1) * cells].reshape(nz, ny, nx)

    # Auto-detect the drop volume fraction: among non-color fields bounded in [0,1], the one most
    # positively correlated with the color field at t=0 (alpha2 = eta == color IC). alpha1 = 1-eta
    # is anti-correlated, so this lands on the drop's volume fraction.
    c0 = field(snap0, c_idx).ravel()
    best_i, best_corr = None, -2.0
    for i in range(nvars):
        if i == c_idx:
            continue
        f = field(snap0, i).ravel()
        if f.min() >= -1e-6 and f.max() <= 1.0 + 1e-4 and f.std() > 0:
            corr = float(np.corrcoef(f, c0)[0, 1])
            if corr > best_corr:
                best_corr, best_i = corr, i
    alpha_idx = best_i

    def band_thickness(f2d):
        """f2d: (ny,nx). Returns interface band thickness in cell widths (dx)."""
        f = f2d[0]  # 2D: squeeze z
        gy, gx = np.gradient(f, dy, dx)
        gmag = np.sqrt(gx**2 + gy**2)
        band_area = float(((f > 0.1) & (f < 0.9)).sum()) * dx * dy
        length = float(gmag.sum()) * dx * dy
        return band_area / length / dx if length > 0 else 0.0

    times, U, th_color, th_alpha, ycen = [], [], [], [], []
    yb3 = y[None, :, None]
    for s in steps:
        snap = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        rho = field(snap, 0) + field(snap, 1)
        vy = field(snap, 3) / rho
        c = np.clip(field(snap, c_idx), 0.0, 1.0)
        csum = c.sum()
        times.append(s * dt)
        U.append((c * vy).sum() / csum)
        ycen.append((c * yb3).sum() / csum)
        th_color.append(band_thickness(field(snap, c_idx)))
        th_alpha.append(band_thickness(field(snap, alpha_idx)))
    times = np.array(times)
    out = {
        "variant": variant,
        "nx": nx,
        "ny": ny,
        "nvars": nvars,
        "c_idx": c_idx,
        "alpha_idx": alpha_idx,
        "alpha_corr": best_corr,
        "dx": dx,
        "t_r": t_r,
        "U_r": U_r,
        "tstar": times / t_r,
        "Ustar": np.array(U) / U_r,
        "th_color": np.array(th_color),
        "th_alpha": np.array(th_alpha),
        "rises": float(ycen[-1]) > float(ycen[0]),
    }
    return out


def summarize(r):
    Us = r["Ustar"]
    pk = int(np.argmax(Us))
    tail = r["tstar"] >= r["tstar"][-1] - 1.0  # final t_r window
    return {
        "peak_Ustar": float(Us[pk]),
        "t_peak": float(r["tstar"][pk]),
        "terminal_Ustar": float(Us[tail].mean()),
        "th_color_0": float(r["th_color"][0]),
        "th_color_end": float(r["th_color"][-1]),
        "th_alpha_0": float(r["th_alpha"][0]),
        "th_alpha_end": float(r["th_alpha"][-1]),
        "t_end_tr": float(r["tstar"][-1]),
    }


def main():
    out_png = sys.argv[1] if len(sys.argv) > 1 else os.path.join(RUNROOT, "recon_thinc_compare.png")
    data, summ = {}, {}
    for v in VARIANTS:
        if not os.path.isdir(os.path.join(RUNROOT, v, "restart_data")):
            print(f"  SKIP {v}: no restart_data (run sweep_recon.py first)")
            continue
        data[v] = load_variant(v)
        summ[v] = summarize(data[v])

    if not data:
        sys.exit("no variants found under runs/recon_tc2/")

    # Report
    print("\n=== THINC (int_comp) on TC2 Ma=20: MUSCL vs MUSCL+THINC ===")
    for v in VARIANTS:
        if v not in data:
            continue
        r, s = data[v], summ[v]
        print(f"\n[{v}]  {LABEL[v]}   nvars={r['nvars']} color_idx={r['c_idx']} alpha_idx={r['alpha_idx']} (corr {r['alpha_corr']:.3f})  run={s['t_end_tr']:.1f} t_r")
        print(f"  U*: peak={s['peak_Ustar']:.4f} @ t*={s['t_peak']:.2f}   terminal={s['terminal_Ustar']:.4f}   rises={r['rises']}")
        print(f"  alpha band thickness (cells): {s['th_alpha_0']:.2f} -> {s['th_alpha_end']:.2f}  (x{s['th_alpha_end'] / s['th_alpha_0']:.2f})")
        print(f"  color band thickness (cells): {s['th_color_0']:.2f} -> {s['th_color_end']:.2f}  (x{s['th_color_end'] / s['th_color_0']:.2f})")

    if len(data) == 2:
        a, b = summ["muscl"], summ["muscl_thinc"]
        print("\n--- THINC effect (muscl_thinc vs muscl) ---")
        print(f"  alpha band end:  {a['th_alpha_end']:.2f} -> {b['th_alpha_end']:.2f} cells  ({100 * (b['th_alpha_end'] - a['th_alpha_end']) / a['th_alpha_end']:+.1f}%)")
        print(f"  color band end:  {a['th_color_end']:.2f} -> {b['th_color_end']:.2f} cells  ({100 * (b['th_color_end'] - a['th_color_end']) / a['th_color_end']:+.1f}%)")
        print(f"  peak U*:         {a['peak_Ustar']:.4f} -> {b['peak_Ustar']:.4f}  ({100 * (b['peak_Ustar'] - a['peak_Ustar']) / a['peak_Ustar']:+.1f}%)")
        print(f"  terminal U*:     {a['terminal_Ustar']:.4f} -> {b['terminal_Ustar']:.4f}  ({100 * (b['terminal_Ustar'] - a['terminal_Ustar']) / a['terminal_Ustar']:+.1f}%)")

    # Overlay figure: U*(t*) on top, band thickness(t*) below.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.8, 8.0), sharex=True)
    for v in VARIANTS:
        if v not in data:
            continue
        r = data[v]
        ax1.plot(r["tstar"], r["Ustar"], "-", color=COLOR[v], lw=1.6, label=LABEL[v])
    ax1.axhline(0.13, ls="--", color="0.5", lw=1.2, label="Nas & Tryggvason peak ~ 0.13")
    ax1.set_ylabel(r"$U^* = U/U_r$")
    ax1.set_title("TC2 Ma=20: THINC interface compression (int_comp) effect")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="upper right", fontsize=8.5)

    for v in VARIANTS:
        if v not in data:
            continue
        r = data[v]
        ax2.plot(r["tstar"], r["th_alpha"], "-", color=COLOR[v], lw=1.7, label=f"{LABEL[v]} -- alpha (vol frac)")
        ax2.plot(r["tstar"], r["th_color"], "--", color=COLOR[v], lw=1.3, label=f"{LABEL[v]} -- color (CSF)")
    ax2.set_xlabel(r"$t^* = t/t_r$")
    ax2.set_ylabel("interface band thickness (cells, $dx$)")
    ax2.set_xlim(left=0.0)
    ax2.grid(alpha=0.3)
    ax2.legend(loc="upper left", fontsize=8.0, ncol=1)
    ax2.text(0.98, 0.03, "solid = volume fraction (THINC-compressed)\ndashed = color function (CSF, never compressed)", transform=ax2.transAxes, ha="right", va="bottom", fontsize=8, color="0.3")

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"\nsaved overlay -> {out_png}")


if __name__ == "__main__":
    main()
