#!/usr/bin/env python3
"""
Independent comparison of mu(T) vs constant-mu thermocapillary drop migration.
Written from scratch; does not consult compare_tc3_visc.py.

For each saved snapshot in both runs, compute:
  - color-weighted drop rise velocity  v = sum(c * v_y) / sum(c)   [mm/s]
  - color-weighted drop height (centroid y measured from y_domain%beg) [mm]
Then compare the two runs at matched step numbers.
"""
import os
import re
import glob
import numpy as np

RUN_MUT = "runs/tc3_2d_muT_fine"
RUN_CONST = "runs/tc3_2d_const_fine"

# variable indices within each 11-variable snapshot
IDX_RHO1 = 0
IDX_RHO2 = 1
IDX_MOMY = 3
IDX_C = 9    # nvars-2 color function (drop=1)
IDX_TS = 10  # nvars-1 temperature scalar


def read_namelist(path):
    """Parse plain 'name = value' Fortran namelist lines we care about."""
    vals = {}
    keymap = {
        "m": "m", "n": "n", "p": "p", "dt": "dt",
        "t_step_save": "t_step_save",
        "x_domain%beg": "xbeg", "x_domain%end": "xend",
        "y_domain%beg": "ybeg", "y_domain%end": "yend",
    }
    with open(path) as fh:
        for line in fh:
            mo = re.match(r"\s*([A-Za-z_%]+)\s*=\s*([-+0-9.eEdD]+)", line)
            if not mo:
                continue
            name, raw = mo.group(1), mo.group(2)
            if name in keymap:
                v = float(raw.replace("D", "e").replace("d", "e"))
                vals[keymap[name]] = v
    return vals


def y_centers(rundir, ny):
    """Cell-center y coords from the last ny+1 cell-boundary entries."""
    yb = np.fromfile(os.path.join(rundir, "restart_data", "lustre_y_cb.dat"),
                     dtype="<f8")
    yb = yb[-(ny + 1):]
    return 0.5 * (yb[:-1] + yb[1:])


def list_saved_steps(rundir, t_step_save):
    steps = []
    for f in glob.glob(os.path.join(rundir, "restart_data", "lustre_*.dat")):
        base = os.path.basename(f)
        mo = re.match(r"lustre_(\d+)\.dat$", base)
        if not mo:
            continue   # skips lustre_y_cb.dat and friends
        s = int(mo.group(1))
        if s % t_step_save == 0:    # only true save-cadence snapshots
            steps.append(s)
    return sorted(steps)


def analyze(rundir):
    nml = read_namelist(os.path.join(rundir, "simulation.inp"))
    nx = int(nml["m"]) + 1
    ny = int(nml["n"]) + 1
    nz = int(nml["p"]) + 1
    cells = nx * ny * nz
    dt = nml["dt"]
    tss = int(nml["t_step_save"])
    ybeg = nml["ybeg"]

    yc = y_centers(rundir, ny)                 # shape (ny,)
    # broadcast y over a (nz, ny, nx) field: y varies along axis 1
    ygrid = yc[None, :, None]                  # (1, ny, 1)

    steps = list_saved_steps(rundir, tss)
    out = {"step": [], "time": [], "v_mm_s": [], "h_mm": []}
    for s in steps:
        path = os.path.join(rundir, "restart_data", f"lustre_{s}.dat")
        flat = np.fromfile(path, dtype="<f8")
        nvars = flat.size // cells
        assert nvars == 11, f"unexpected nvars={nvars} in {path}"

        def var(i):
            return flat[i * cells:(i + 1) * cells].reshape(nz, ny, nx)

        rho = var(IDX_RHO1) + var(IDX_RHO2)
        vy = var(IDX_MOMY) / rho
        c = var(IDX_C)

        csum = c.sum()
        v = (c * vy).sum() / csum                      # m/s
        ycen = (c * ygrid).sum() / csum                # m (absolute)
        h = (ycen - ybeg)                              # height above floor, m

        out["step"].append(s)
        out["time"].append(s * dt)
        out["v_mm_s"].append(v * 1e3)
        out["h_mm"].append(h * 1e3)

    for k in out:
        out[k] = np.asarray(out[k])
    return out, nml


def smooth_trend(y, win=7):
    """Centered moving-average trend (odd window)."""
    if win % 2 == 0:
        win += 1
    half = win // 2
    n = len(y)
    tr = np.empty(n)
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        tr[i] = y[a:b].mean()
    return tr


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    os.chdir(base)

    mut, nml_m = analyze(RUN_MUT)
    con, nml_c = analyze(RUN_CONST)

    # match on common step numbers
    common = np.intersect1d(mut["step"], con["step"])
    im = {s: k for k, s in enumerate(mut["step"])}
    ic = {s: k for k, s in enumerate(con["step"])}

    steps = common
    t = np.array([mut["time"][im[s]] for s in steps])
    vm = np.array([mut["v_mm_s"][im[s]] for s in steps])
    vc = np.array([con["v_mm_s"][ic[s]] for s in steps])
    hm = np.array([mut["h_mm"][im[s]] for s in steps])
    hc = np.array([con["h_mm"][ic[s]] for s in steps])

    print(f"grid nx={int(nml_m['m'])+1} ny={int(nml_m['n'])+1} "
          f"nz={int(nml_m['p'])+1}  dt={nml_m['dt']:.4e}  "
          f"t_step_save={int(nml_m['t_step_save'])}")
    print(f"matched snapshots: {len(steps)}  "
          f"steps {steps[0]}..{steps[-1]}  "
          f"t {t[0]*1e3:.3f}..{t[-1]*1e3:.3f} ms\n")

    # divergence ratio (guard tiny denominators)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(np.abs(vc) > 1e-12, (vm - vc) / vc * 100.0, np.nan)

    # report at ~50%, ~75%, 100% through the run (by index)
    n = len(steps)
    idxs = {"~50%": n // 2, "~75%": (3 * n) // 4, "final": n - 1}
    print("frac    step      t(ms)   h_muT   h_const   v_muT     v_const    (vm-vc)/vc")
    for label, i in idxs.items():
        print(f"{label:6s} {steps[i]:8d} {t[i]*1e3:7.2f} "
              f"{hm[i]:7.3f} {hc[i]:8.3f}  "
              f"{vm[i]:8.4f}  {vc[i]:8.4f}   {ratio[i]:+8.2f}%")

    # noise estimate: residual of v about a smooth trend (both runs)
    trend_m = smooth_trend(vm)
    trend_c = smooth_trend(vc)
    res_m = vm - trend_m
    res_c = vc - trend_c
    noise_m = res_m.std()
    noise_c = res_c.std()
    print(f"\nacoustic ripple (std of v about 7-pt trend): "
          f"muT={noise_m:.4f} mm/s   const={noise_c:.4f} mm/s")

    # signed difference vs noise, using smoothed velocities to kill ripple
    diff_raw = vm - vc                       # raw signed diff (mm/s)
    diff_sm = trend_m - trend_c              # smoothed signed diff (mm/s)
    # combined noise floor on the difference of two independent ripples
    noise_diff = np.hypot(noise_m, noise_c)
    print(f"noise floor on (v_muT - v_const): ~{noise_diff:.4f} mm/s")

    # final-portion mean difference (last quarter) for robustness
    q = slice(3 * n // 4, n)
    print(f"\nmean signed diff (smoothed) over last quarter: "
          f"{diff_sm[q].mean():+.4f} mm/s "
          f"(raw {diff_raw[q].mean():+.4f} mm/s, "
          f"noise floor {noise_diff:.4f})")
    print(f"mean %% divergence over last quarter: "
          f"{np.nanmean(ratio[q]):+.2f}%")

    # does divergence grow with height? correlate smoothed signed diff vs height
    havg = 0.5 * (hm + hc)
    # use second half where drop is clearly rising and signal is established
    half = slice(n // 2, n)
    cc = np.corrcoef(havg[half], diff_sm[half])[0, 1]
    # linear slope of signed diff vs height (mm/s per mm risen)
    A = np.vstack([havg[half], np.ones(half.stop - half.start
                                       if half.stop else n - n // 2)]).T
    slope, intercept = np.linalg.lstsq(A, diff_sm[half], rcond=None)[0]
    print(f"\ndivergence-vs-height (2nd half): corr={cc:+.3f}  "
          f"slope={slope:+.5f} mm/s per mm risen")
    print(f"height span 2nd half: {havg[n//2]:.3f} -> {havg[-1]:.3f} mm")

    # dump full table for transparency
    print("\nfull matched series (step, t_ms, h_avg_mm, v_muT, v_const, "
          "diff_sm, %%):")
    for i in range(n):
        print(f"  {steps[i]:8d} {t[i]*1e3:7.2f} {havg[i]:7.3f} "
              f"{vm[i]:8.4f} {vc[i]:8.4f} {diff_sm[i]:+8.4f} "
              f"{ratio[i]:+7.2f}")


if __name__ == "__main__":
    main()
