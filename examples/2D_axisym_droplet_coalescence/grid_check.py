#!/usr/bin/env python3
"""Grid designer/verifier for the axisymmetric film-zoom probe.

Replicates MFC's hyperbolic-tangent grid stretching EXACTLY
(src/pre_process/m_grid.f90, s_generate_serial_grid):

  pre-stretch (x):  dx0 = (x_end - x_beg)/(m+1); boundaries x_cb[i] = x_beg + dx0*i
  pre-stretch (y, axisym half-cell offset when cyl_coord and y_beg == 0):
                    dy0 = (y_end - y_beg)/(2n+1); y_cb[-1]=0, y_cb[i-1] = dy0*(2i-1)
  stretch (per loop, coordinates normalized by total length L; x_a, x_b normalized once):
      x <- (x/a) * ( a + log(cosh(a*(x - x_a))) + log(cosh(a*(x - x_b)))
                       - 2*log(cosh(a*(x_b - x_a)/2)) )
  then de-normalize by L. Note: x = 0 is a fixed point and the map slope there is
  exactly 1 when x_a = -x_b, so the CELL SIZE AT THE PLANE equals the pre-stretch
  spacing; stretching only expands the far field.

Run:  python3 grid_check.py     (checks the pilot and production parameter sets)
"""

import json
import math
import sys


def logcosh(z):
    # numerically safe log(cosh(z)) for large |z|
    az = abs(z)
    if az > 30.0:
        return az - math.log(2.0)
    return math.log(math.cosh(z))


def stretch(cb, L, a, xa, xb, loops):
    """Apply MFC's tanh stretch to boundary coords cb (list, physical units)."""
    cbn = [c / L for c in cb]
    xan, xbn = xa / L, xb / L
    c3 = logcosh(a * (xbn - xan) / 2.0)
    for _ in range(loops):
        cbn = [
            (c / a) * (a + logcosh(a * (c - xan)) + logcosh(a * (c - xbn)) - 2.0 * c3)
            for c in cbn
        ]
    return [c * L for c in cbn]


def compute_grid(p):
    """p: dict with x_beg,x_end,m,a_x,x_a,x_b,loops_x, y_end,n,a_y,y_a,y_b,loops_y.
    Returns dict of boundary arrays + diagnostics. Assumes x_beg=0 half-domain,
    y axisym half-cell offset (cyl_coord, y_beg=0)."""
    m, n = p["m"], p["n"]
    # x boundaries: i = -1..m -> m+2 boundaries, uniform
    dx0 = (p["x_end"] - p["x_beg"]) / (m + 1)
    xcb = [p["x_beg"] + dx0 * i for i in range(0, m + 2)]  # index shift: [i-1] i=0..m+1
    xcb = stretch(xcb, abs(xcb[-1] - xcb[0]), p["a_x"], p["x_a"], p["x_b"], p["loops_x"])
    # y boundaries: axisym offset grid
    dy0 = (p["y_end"] - 0.0) / (2 * n + 1)
    ycb = [0.0] + [dy0 * (2 * i - 1) for i in range(1, n + 1)] + [p["y_end"]]
    ycb = stretch(ycb, abs(ycb[-1] - ycb[0]), p["a_y"], p["y_a"], p["y_b"], p["loops_y"])
    dxs = [xcb[i + 1] - xcb[i] for i in range(len(xcb) - 1)]
    dys = [ycb[i + 1] - ycb[i] for i in range(len(ycb) - 1)]
    return {"xcb": xcb, "ycb": ycb, "dx": dxs, "dy": dys}


def band_max(cb, d, lo, hi):
    vals = [d[i] for i in range(len(d)) if cb[i] >= lo and cb[i + 1] <= hi]
    return max(vals) if vals else float("nan")


def growth_max(cb, d, size_below, skip_first=False):
    """max cell-to-cell growth ratio among cells smaller than size_below.
    skip_first exempts the axis half-cell (dy0 -> 2*dy0 is by construction)."""
    g = 1.0
    for i in range(1 if skip_first else 0, len(d) - 1):
        if d[i] < size_below and d[i + 1] < size_below and d[i] > 0:
            r = d[i + 1] / d[i]
            g = max(g, r, 1.0 / r if r > 0 else 1.0)
    return g


def check(name, p, dx_fine):
    D = 306e-6
    g = compute_grid(p)
    dx, dy, xcb, ycb = g["dx"], g["dy"], g["xcb"], g["ycb"]
    cells = (p["m"] + 1) * (p["n"] + 1)
    rep = {
        "name": name,
        "cells": cells,
        "dx_plane": dx[0],
        "dx_fine_band_max(0..4um)": band_max(xcb, dx, 0.0, 4e-6),
        "x_growth_max(dx<1um)": growth_max(xcb, dx, 1e-6),
        "dx_max_over_drop(x<350um)": band_max(xcb, dx, 0.0, 350e-6),
        "Lx_phys": xcb[-1],
        "dy_axis": dy[0],
        "dy_rim_max(0..70um)": band_max(ycb, dy, 0.0, 70e-6),
        "y_growth_max(dy<1um)": growth_max(ycb, dy, 1e-6, skip_first=True),
        "dy_max_over_drop(y<350um)": band_max(ycb, dy, 0.0, 350e-6),
        "Ly_phys": ycb[-1],
        "film_aspect(dy0/dx0)": dy[0] / dx[0],
    }
    crit = {
        "a_dx_plane": rep["dx_plane"] <= 1.10 * dx_fine,
        "b_fine_band": rep["dx_fine_band_max(0..4um)"] <= 1.2 * dx_fine,
        "c_x_growth": rep["x_growth_max(dx<1um)"] <= 1.06,
        "d_drop_resolved_x": rep["dx_max_over_drop(x<350um)"] <= 4e-6,
        "d_extent_x": rep["Lx_phys"] >= 1.5 * D,
        "e_rim_dy": rep["dy_rim_max(0..70um)"] <= 0.6e-6,
        "e_y_growth": rep["y_growth_max(dy<1um)"] <= 1.06,
        "e_drop_resolved_y": rep["dy_max_over_drop(y<350um)"] <= 4e-6,
        "e_extent_y": rep["Ly_phys"] >= 1.5 * D,
        "f_cells": cells <= (320_000 if name == "pilot" else 650_000),
        "f_aspect": rep["film_aspect(dy0/dx0)"] <= 7.0,
    }
    # cost table (artificial c = 100 m/s; |u| ~ 1)
    for cfl in (0.1, 0.2, 0.3):
        dt = cfl * min(dx[0], dy[0]) / 101.0
        steps = 4.0e-4 / dt
        rep[f"dt_cfl{cfl}"] = dt
        rep[f"steps_cfl{cfl}"] = int(steps)
        rep[f"hours_cfl{cfl}_at_0.01s/step"] = steps * 0.01 / 3600.0
    ok = all(crit.values())
    return ok, rep, crit


# Parameter sets (must match case.py's emitted values for the same --dx-fine)
def params(dx_fine):
    # LOCKED parameters - must match case.py exactly (cross-checked in __main__).
    # Axial: dx(plane) = Lx_pre/(m+1) = dx_fine (map slope 1 at x=0); fine band
    # +/-4 um; two loops expand the far field to ~1.5 mm. Radial: axisym offset
    # grid, interior dy = 2*Ly_pre/(2n+1) = 0.50 um uniform through the rim zone
    # (band +/-90 um), one loop to the far field.
    Lx_pre = 80e-6
    m = int(round(Lx_pre / dx_fine)) - 1
    Ly_pre = 260e-6
    n = 520
    return {
        "x_beg": 0.0, "x_end": Lx_pre, "m": m,
        "a_x": 10.0, "x_a": -4e-6, "x_b": 4e-6, "loops_x": 2,
        "y_end": Ly_pre, "n": n,
        "a_y": 8.0, "y_a": -90e-6, "y_b": 90e-6, "loops_y": 1,
    }


if __name__ == "__main__":
    fail = False
    for name, dxf in (("pilot", 150e-9), ("production", 75e-9)):
        p = params(dxf)
        ok, rep, crit = check(name, p, dxf)
        print(f"\n=== {name} (dx_fine={dxf*1e9:.0f} nm)  m={p['m']} n={p['n']} ===")
        for k, v in rep.items():
            if isinstance(v, float):
                print(f"  {k:38s} {v:.6g}")
            else:
                print(f"  {k:38s} {v}")
        for k, v in crit.items():
            print(f"  [{'PASS' if v else 'FAIL'}] {k}")
        fail |= not ok
    # cross-check: case.py must emit the same grid parameters
    import subprocess
    for dxf in (150e-9, 75e-9):
        out = subprocess.run(
            [sys.executable, "case.py", "--mode", "probe", "--dx-fine", str(dxf)],
            capture_output=True, text=True, cwd=__import__("os").path.dirname(
                __import__("os").path.abspath(__file__)))
        emitted = json.loads(out.stdout)
        p = params(dxf)
        agree = (emitted["m"] == p["m"] and emitted["n"] == p["n"]
                 and emitted["a_x"] == p["a_x"] and emitted["x_b"] == p["x_b"]
                 and emitted["a_y"] == p["a_y"] and emitted["y_b"] == p["y_b"]
                 and emitted["loops_y"] == p["loops_y"]
                 and emitted["x_domain%end"] == p["x_end"]
                 and emitted["y_domain%end"] == p["y_end"])
        print(f"case.py cross-check dx_fine={dxf*1e9:.0f}nm: {'PASS' if agree else 'FAIL'}")
        fail |= not agree
    sys.exit(1 if fail else 0)
