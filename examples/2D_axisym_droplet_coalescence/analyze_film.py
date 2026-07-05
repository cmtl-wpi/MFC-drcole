#!/usr/bin/env python3
"""Film-zoom probe post-processor: extract H(r,t), H_center, H_min, r_h, p(0,0).

Reads MFC post_process output (format=1, silo_hdf5/) OR the local pipeline's
xdmf/step_*.h5 conversion if present (same convention as the earlier
3D_droplet_coalescence analysis). Writes film_metrics.csv and film_metrics.png.

Status: reader auto-detection and the metric extraction are implemented; the
silo variable paths may need one-time adjustment against the first real output
(marked TODO below) since MFC silo layouts differ by version/parallel_io.

Usage: python3 analyze_film.py <run_dir>
"""

import csv
import glob
import os
import sys

import numpy as np

R = 153e-6
D = 2.0 * R
UR = 0.98984  # relative velocity (mirror pair)
TAU = D / UR


def find_frames(run_dir):
    """Return sorted list of (time_or_index, path, kind)."""
    h5s = sorted(glob.glob(os.path.join(run_dir, "xdmf", "step_*.h5")))
    if h5s:
        return [(i, p, "xdmf") for i, p in enumerate(h5s)]
    silos = sorted(glob.glob(os.path.join(run_dir, "silo_hdf5", "root", "*.silo")))
    if silos:
        return [(i, p, "silo") for i, p in enumerate(silos)]
    raise SystemExit(f"no xdmf/step_*.h5 or silo_hdf5 frames under {run_dir}")


def read_frame(path, kind):
    """Return (x_axial, y_radial, alpha1[y, x], pres[y, x])."""
    import h5py

    with h5py.File(path, "r") as f:
        if kind == "xdmf":
            x = f["x"][:]
            y = f["y"][:]
            a = f["alpha1"][:]
            p = f["pres"][:]
            if a.ndim == 3:  # (z,y,x) with singleton z from 2D runs
                a, p = a[0], p[0]
            return x, y, a, p
        # TODO(silo): confirm dataset paths on first real output. Typical MFC
        # silo: /DOMAIN_/mesh coords + /DOMAIN_/alpha1 etc. Adjust here once.
        raise SystemExit("silo reader needs the one-time path adjustment (see TODO)")


def film_metrics(x, y, a, rim_rmax=90e-6):
    """H(r) = 2 * x-position of the alpha=0.5 crossing nearest the plane."""
    Hs = np.full(len(y), np.nan)
    for j, yj in enumerate(y):
        if yj > rim_rmax:
            break
        col = a[j, :]
        idx = np.where(col >= 0.5)[0]
        if idx.size == 0:
            continue
        i = idx[0]  # first liquid cell moving away from the plane
        if i == 0:
            Hs[j] = 0.0  # interface reached the plane: numerical merge
            continue
        # linear interpolation of the 0.5 crossing between cells i-1, i
        f = (0.5 - col[i - 1]) / max(col[i] - col[i - 1], 1e-30)
        Hs[j] = 2.0 * (x[i - 1] + f * (x[i] - x[i - 1]))
    valid = ~np.isnan(Hs)
    if not valid.any():
        return np.nan, np.nan, np.nan, Hs
    H_center = Hs[valid][0]
    jmin = np.nanargmin(np.where(valid, Hs, np.inf))
    return H_center, Hs[jmin], y[jmin], Hs


def main(run_dir):
    frames = find_frames(run_dir)
    rows = []
    for i, path, kind in frames:
        x, y, a, p = read_frame(path, kind)
        Hc, Hmin, rh, _ = film_metrics(x, y, a)
        p00 = float(p[0, 0])
        rows.append((i, Hc, Hmin, rh, p00))
        print(f"frame {i:4d}: H_center={Hc if Hc==Hc else float('nan'):.3e} "
              f"H_min={Hmin:.3e} r_h={rh:.3e} p(0,0)={p00:.1f}")
    out = os.path.join(run_dir, "film_metrics.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["frame", "H_center_m", "H_min_m", "r_h_m", "p00_Pa"])
        w.writerows(rows)
    print(f"wrote {out}")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        arr = np.array(rows, dtype=float)
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))
        for ax, col, label in zip(
            axes.flat, (1, 2, 3, 4),
            ("H_center (m)", "H_min (m)", "r_h (m)", "p(0,0) (Pa)")):
            ax.plot(arr[:, 0], arr[:, col], "k-")
            ax.set_xlabel("frame")
            ax.set_ylabel(label)
            if col in (1, 2):
                ax.set_yscale("log")
        fig.suptitle("film-zoom probe: gap/dimple metrics (cf. Chen & Yang Figs. 13/15)")
        fig.tight_layout()
        fig.savefig(os.path.join(run_dir, "film_metrics.png"), dpi=140)
        print("wrote film_metrics.png")
    except Exception as e:  # matplotlib optional
        print(f"(figure skipped: {e})")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
