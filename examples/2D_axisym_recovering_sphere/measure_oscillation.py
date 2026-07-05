#!/usr/bin/env python3
"""Measure the recovering-sphere validation results.

Oscillate mode: extract the l=2 oscillation period from the drop-centre
pressure probe (D/probe1_prim.dat) and compare against the sphere and
planar-drop (cylinder) Rayleigh-Lamb targets. The probe signal mixes fast
acoustic ringing with the slow capillary mode, so it is smoothed with a
moving-average window much longer than the acoustic period and much shorter
than the capillary period, then the period is read from zero crossings of the
mean-subtracted signal (median crossing-to-crossing interval x 2).

Static mode: time-average the drop-centre minus far-gas pressure over the
second half of the record and compare to the Laplace jump 2*sigma/R.

Usage: python3 measure_oscillation.py <run_dir> [--mode oscillate|static]
"""

import argparse
import sys

import numpy as np

RHO_L = 763.0
RHO_G = RHO_L / 666.0
SIGMA = 0.0266
R = 50e-6

T_SPHERE = 2.0 * np.pi / np.sqrt(24.0 * SIGMA / ((3 * RHO_L + 2 * RHO_G) * R**3))
T_CYL = 2.0 * np.pi / np.sqrt(6.0 * SIGMA / ((RHO_L + RHO_G) * R**3))
JUMP = 2.0 * SIGMA / R


def smooth(x, w):
    k = np.ones(w) / w
    return np.convolve(x, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--mode", default="oscillate", choices=["oscillate", "static"])
    a = ap.parse_args()

    d1 = np.loadtxt(f"{a.run_dir}/D/probe1_prim.dat")  # t, rho, u, p (drop centre)
    t, p1 = d1[:, 0], d1[:, 3]

    if a.mode == "static":
        d2 = np.loadtxt(f"{a.run_dir}/D/probe2_prim.dat")
        p2 = d2[:, 3]
        half = len(t) // 2
        jump = np.mean(p1[half:]) - np.mean(p2[half:])
        err = jump / JUMP - 1
        print(f"time-averaged Laplace jump = {jump:.1f} Pa  target 2*sigma/R = {JUMP:.1f} Pa")
        print(f"ratio = {jump/JUMP:.4f}   ({err*100:+.1f}%)   "
              f"[unfixed axisym held ~sigma/R = {JUMP/2:.0f} Pa or aborted]")
        print(f"VERDICT: {'PASS' if abs(err) < 0.10 else 'FAIL'} (+/-10%)")
        sys.exit(0 if abs(err) < 0.10 else 1)

    # oscillate: smooth away acoustics (window ~ 2 us of samples)
    dt_s = np.median(np.diff(t))
    w = max(3, int(2e-6 / dt_s))
    ps = smooth(p1, w)
    # drop the transient first 20 us and the edge-contaminated tail
    i0 = np.searchsorted(t, 2e-5)
    i1 = len(t) - w
    sig = ps[i0:i1] - smooth(ps, 5 * w)[i0:i1]  # detrend against slow damping
    ts = t[i0:i1]
    # zero crossings
    s = np.sign(sig)
    cross = np.where(np.diff(s) != 0)[0]
    if len(cross) < 3:
        print("too few crossings - inspect signal manually")
        sys.exit(2)
    half_periods = np.diff(ts[cross])
    T_meas = 2.0 * np.median(half_periods)
    print(f"samples={len(t)}  smooth window={w}  crossings={len(cross)}")
    print(f"measured period  T = {T_meas*1e6:.2f} us")
    print(f"sphere target      = {T_SPHERE*1e6:.2f} us  (ratio {T_meas/T_SPHERE:.4f})")
    print(f"planar (wrong)     = {T_CYL*1e6:.2f} us  (ratio {T_meas/T_CYL:.4f})")
    is_sphere = abs(T_meas / T_SPHERE - 1) < 0.07
    is_cyl = abs(T_meas / T_CYL - 1) < 0.07
    verdict = "SPHERE (fix validated)" if is_sphere and not is_cyl else \
              "PLANAR (azimuthal curvature still missing!)" if is_cyl and not is_sphere else \
              "AMBIGUOUS - inspect"
    print(f"VERDICT: {verdict}")
    sys.exit(0 if (is_sphere and not is_cyl) else 1)


if __name__ == "__main__":
    main()
