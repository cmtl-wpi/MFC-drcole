#!/usr/bin/env python3
# 2D analytic validation of the bulk-conduction operator. A single 2D Fourier mode of the independent
# temperature scalar T_s decays on a doubly-periodic square. T_s is imposed directly (decoupled from
# density, so u stays 0) and therefore obeys dT/dt = alpha*lap(T) exactly, with closed-form solution
#     T_s(x,y,t) = T0 * (1 + eps * sin(kx x) * sin(ky y) * exp(-alpha*(kx^2 + ky^2)*t)).
# Different x/y wavenumbers (nx=1, ny=2) make kx != ky, so the test exercises the x- and y-conduction
# fluxes separately -- a genuine 2D check, not a 1D test in disguise. This is the 2D counterpart of
# verify_1d_thermal_scalar.py: the cell-resolved field must match the closed form and the mode must
# decay at alpha*(kx^2 + ky^2).
#
#   ./mfc.sh run examples/2D_thermocapillary_migration/verify_2d_conduction.py -n 4 -t pre_process simulation
#   python3 examples/2D_thermocapillary_migration/verify_2d_conduction.py --artifact <rundir>
import json
import math
import sys

# Shared physical constants (module level so the generator and the measurer agree).
Lx = Ly = 1.0
mx = my = 127  # 128 x 128 cells
nx_mode, ny_mode = 1, 2  # Fourier mode numbers -> kx != ky exercises x- and y-conduction separately
gam, p_inf, cv, p0, T0 = 2.0, 100.0, 12.5, 25.0, 10.0
eps = 0.05  # mode amplitude (small; density stays exactly uniform regardless since T_s is decoupled)
alpha = 0.05  # target thermal diffusivity
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * T0)  # = 1.0 (uniform -> uniform alpha = k/(rho*cp))
k_therm = alpha * rho0 * cv * gam  # alpha = k/(rho*cp), cp = gam*cv
kx = nx_mode * 2.0 * math.pi / Lx
ky = ny_mode * 2.0 * math.pi / Ly
rate_an = alpha * (kx**2 + ky**2)  # analytic decay rate of the mode

# Conserved-var layout, model_eqns=2 / num_fluids=1 / 2D / thermal_scalar:
#   0 rho, 1 mom_x, 2 mom_y, 3 E, 4 alpha, 5 T_s
TS_IDX0 = 5


def make_case():
    c0 = (gam * (p0 + p_inf)) ** 0.5
    dx = Lx / (mx + 1)
    dt = min(0.4 * dx / c0, 0.35 * dx**2 / (2.0 * 2 * alpha))  # 2 = spatial dims in the diffusion limit
    t_end = 2.0 / rate_an  # ~2 e-folds of decay
    t_step_stop = int(round(t_end / dt))
    ts_expr = f"{T0}*(1.0 + {eps}*sin({kx:.9f}*x)*sin({ky:.9f}*y))"  # T_s imposed directly (no 'e' literal)
    data = {
        "run_time_info": "T",
        "x_domain%beg": 0.0,
        "x_domain%end": Lx,
        "y_domain%beg": 0.0,
        "y_domain%end": Ly,
        "m": mx,
        "n": my,
        "p": 0,
        "cyl_coord": "F",
        "dt": dt,
        "t_step_start": 0,
        "t_step_stop": t_step_stop,
        "t_step_save": max(1, t_step_stop // 40),
        "model_eqns": 2,
        "num_fluids": 1,
        "mpp_lim": "F",
        "mixture_err": "T",
        "time_stepper": 3,
        "weno_order": 5,
        "weno_eps": 1e-16,
        "mapped_weno": "T",
        "null_weights": "F",
        "mp_weno": "T",
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "bc_x%beg": -1,  # doubly periodic
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "num_patches": 1,
        "thermal_scalar": "T",
        "thermal_conduction": "T",
        "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
        "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
        "fluid_pp(1)%cv": cv,
        "fluid_pp(1)%k_therm": k_therm,
        "patch_icpp(1)%geometry": 3,  # 2D rectangle spanning the domain
        "patch_icpp(1)%x_centroid": Lx / 2,
        "patch_icpp(1)%y_centroid": Ly / 2,
        "patch_icpp(1)%length_x": Lx,
        "patch_icpp(1)%length_y": Ly,
        "patch_icpp(1)%vel(1)": 0.0,
        "patch_icpp(1)%vel(2)": 0.0,
        "patch_icpp(1)%pres": p0,
        "patch_icpp(1)%alpha_rho(1)": rho0,  # uniform density
        "patch_icpp(1)%alpha(1)": 1.0,
        "patch_icpp(1)%T_temp_val": ts_expr,  # T_s imposed directly
        "format": 1,
        "precision": 2,
        "prim_vars_wrt": "T",
        "cons_vars_wrt": "T",
        "parallel_io": "T",
    }
    print(json.dumps(data))


def artifact(rundir):
    import glob
    import os
    import re

    import numpy as np
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ncx, ncy = mx + 1, my + 1
    rd = os.path.join(rundir, "restart_data")
    inp = {}
    for line in open(os.path.join(rundir, "simulation.inp")):
        if "=" in line:
            kk, vv = line.split("=", 1)
            inp[kk.strip().lower()] = vv.strip().rstrip(",")
    dt = float(inp["dt"])
    steps = sorted(int(g.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (g := re.search(r"lustre_(\d+)\.dat$", os.path.basename(f))))
    xc = (np.arange(ncx) + 0.5) * (Lx / ncx)
    yc = (np.arange(ncy) + 0.5) * (Ly / ncy)
    X, Y = np.meshgrid(xc, yc)  # (ncy, ncx): x varies along axis 1, y along axis 0
    basis = np.sin(kx * X) * np.sin(ky * Y)
    norm = (basis**2).sum()

    def record(s):
        a = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64)
        return a.reshape(a.size // (ncx * ncy), ncy, ncx)  # MFC stores x innermost -> [var, y, x]

    show = len(steps) // 2  # an intermediate output where the mode is still clearly visible
    t, amp, ferr, umax = [], [], [], []
    field_mfc = field_an = None
    t_show = 0.0
    for i, s in enumerate(steps):
        rec = record(s)
        ts, rho, mx_, my_ = rec[TS_IDX0], rec[0], rec[1], rec[2]
        ts_t = s * dt
        an = T0 * (1.0 + eps * np.exp(-rate_an * ts_t) * basis)
        amp.append((ts * basis).sum() / norm)
        ferr.append(np.abs(ts - an).max() / (eps * T0) * 100.0)
        umax.append(float(np.max(np.sqrt(mx_**2 + my_**2) / rho)))
        t.append(ts_t)
        if i == show:
            field_mfc, field_an, t_show = ts, an, ts_t
    t, amp, ferr = np.array(t), np.array(amp), np.array(ferr)
    rate = -np.polyfit(t, np.log(np.abs(amp / amp[0])), 1)[0]
    rate_err = 100 * abs(rate / rate_an - 1)
    max_ferr, max_u = float(ferr.max()), float(max(umax))

    fig = plt.figure(figsize=(13, 7))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1.0])
    vmin, vmax = float(field_an.min()), float(field_an.max())
    ax0, ax1, ax2 = (fig.add_subplot(gs[0, c]) for c in range(3))
    ax0.pcolormesh(xc, yc, field_mfc, vmin=vmin, vmax=vmax, cmap="inferno", shading="auto")
    ax0.set_title(f"MFC  T_s(x,y),  t={t_show:.3f}")
    im1 = ax1.pcolormesh(xc, yc, field_an, vmin=vmin, vmax=vmax, cmap="inferno", shading="auto")
    ax1.set_title("analytic  T0(1+ε·sin·sin·e$^{-α(k_x^2+k_y^2)t}$)")
    fig.colorbar(im1, ax=[ax0, ax1], shrink=0.85, label="T_s")
    diff = field_mfc - field_an
    dm = float(np.abs(diff).max()) or 1e-12
    im2 = ax2.pcolormesh(xc, yc, diff, vmin=-dm, vmax=dm, cmap="coolwarm", shading="auto")
    ax2.set_title(f"MFC − analytic  (max {dm:.2e})")
    fig.colorbar(im2, ax=ax2, shrink=0.85)
    for a in (ax0, ax1, ax2):
        a.set_xlabel("x")
        a.set_ylabel("y")
        a.set_aspect("equal")

    axd = fig.add_subplot(gs[1, :2])
    axd.plot(t, amp, "-", color="C0", lw=2.0, label="MFC (mode projection)")
    axd.plot(t, amp[0] * np.exp(-rate_an * t), "--", color="k", lw=1.5, label=f"analytic  A0·e$^{{-{rate_an:.2f}\\,t}}$")
    axd.set(xlabel="t", ylabel="mode amplitude A(t)", title=f"mode decay: MFC {rate:.3f} vs analytic {rate_an:.3f}  ({rate_err:.1f}%)")
    axd.legend(fontsize=9)
    axe = fig.add_subplot(gs[1, 2])
    axe.plot(t, ferr, "C3-", lw=1.6)
    axe.set(xlabel="t", ylabel="max|T$_{sim}$−T$_{an}$| (% of εT0)", title=f"pointwise error (peak {max_ferr:.2f}%)")
    fig.suptitle(f"verify_2d_conduction: 2D scalar mode (nx={nx_mode}, ny={ny_mode}) vs analytic  —  max|u|={max_u:.1e}", fontweight="bold")
    fig.tight_layout()

    art = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "verify")
    os.makedirs(art, exist_ok=True)
    png = os.path.join(art, "conduction_2d.png")
    fig.savefig(png, dpi=130)
    plt.close(fig)
    sp = os.path.join(art, "summary.json")
    summary = json.load(open(sp)) if os.path.exists(sp) else {}
    summary["conduction_2d"] = {"measured_rate": rate, "analytic_rate": rate_an, "rate_error_pct": rate_err, "field_err_peak_pct": max_ferr, "max_u": max_u, "kx": kx, "ky": ky, "n_outputs": len(steps), "png": "conduction_2d.png"}
    json.dump(summary, open(sp, "w"), indent=2)
    print(f"  field err at t=0 (reshape sanity) = {ferr[0]:.3f}%   (must be ~0 if x/y order is right)")
    print(f"  measured rate {rate:.4f} vs analytic {rate_an:.4f} ({rate_err:.1f}%)   field err peak {max_ferr:.2f}%   max|u| {max_u:.2e}")
    print(f"  wrote {png}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--artifact":
        artifact(sys.argv[2])
    else:
        make_case()
