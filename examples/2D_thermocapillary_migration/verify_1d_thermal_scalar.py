#!/usr/bin/env python3
# Analytic verification of the independent temperature scalar (thermal_scalar): decay of a
# sinusoidal T_s mode on a PERIODIC 1D domain at UNIFORM density. This is the cleanest possible
# demonstration that T_s is decoupled from the EOS -- unlike verify_1d_diffusion.py, the
# temperature perturbation is imposed directly on the advected scalar (patch_icpp%T_temp_val),
# NOT through a density gradient, so the density and pressure stay exactly uniform and there is no
# flow to advect or contaminate the signal.
#
# T_s(x,0) = T0*(1 + eps*sin(2*pi*x/L)). With thermal_conduction also enabled, T_s diffuses at the
# thermal diffusivity alpha = k/(rho*cp), cp = gamma*cv, so the fundamental mode amplitude obeys
# A(t) = A(0)*exp(-alpha*kappa^2*t), kappa = 2*pi/L. Reading the T_s field directly (no EOS
# inversion), projecting onto sin(kappa x), and fitting ln(A) vs t must recover alpha*kappa^2.
#
# Usage:
#   ./mfc.sh run examples/2D_thermocapillary_migration/verify_1d_thermal_scalar.py -n 1 -t pre_process simulation
#   python3 examples/2D_thermocapillary_migration/verify_1d_thermal_scalar.py --measure <rundir>
import json
import sys

# Shared physical constants (module level so both the generator and the measurer agree)
L = 1.0
m = 399  # 400 cells (raised from 200 for a finer field-vs-analytic curve fit)
gam, p_inf, cv, p0, T0 = 2.0, 100.0, 12.5, 25.0, 10.0
eps = 0.05  # T_s perturbation amplitude (small, but density stays exactly uniform regardless)
alpha = 0.05  # target thermal diffusivity
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * T0)  # = 1.0 (uniform)
k_therm = alpha * rho0 * cv * gam  # alpha = k/(rho*cp), cp = gam*cv

# Conserved-variable layout for model_eqns=2, num_fluids=1, 1D with thermal_scalar:
#   1: rho, 2: mom_x, 3: E, 4: alpha (volume fraction), 5: T_s
TS_IDX0 = 4  # zero-based index of T_s in the restart record


def make_case():
    import math

    kappa = 2.0 * math.pi / L
    c0 = (gam * (p0 + p_inf)) ** 0.5
    dx = L / (m + 1)
    dt = min(0.4 * dx / c0, 0.35 * dx**2 / (2.0 * alpha))
    t_end = 0.6 / (alpha * kappa**2)  # ~0.6 e-folding times of decay
    t_step_stop = int(round(t_end / dt))
    # Temperature imposed DIRECTLY on the scalar (no 'e' literal -> avoids the Euler footgun):
    ts_expr = f"{T0}*(1.0 + {eps}*sin({kappa:.9f}*x))"
    data = {
        "run_time_info": "T",
        "x_domain%beg": 0.0,
        "x_domain%end": L,
        "m": m,
        "n": 0,
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
        "bc_x%beg": -1,  # periodic
        "bc_x%end": -1,
        "num_patches": 1,
        "thermal_scalar": "T",
        "thermal_conduction": "T",
        "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
        "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
        "fluid_pp(1)%cv": cv,
        "fluid_pp(1)%k_therm": k_therm,
        "patch_icpp(1)%geometry": 1,
        "patch_icpp(1)%x_centroid": L / 2,
        "patch_icpp(1)%length_x": L,
        "patch_icpp(1)%vel(1)": 0.0,
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


def measure(rundir):
    import glob
    import math
    import os
    import re

    import numpy as np

    cells = m + 1
    rd = os.path.join(rundir, "restart_data")
    inp = {}
    for line in open(os.path.join(rundir, "simulation.inp")):
        if "=" in line:
            kk, vv = line.split("=", 1)
            inp[kk.strip().lower()] = vv.strip().rstrip(",")
    dt = float(inp["dt"])
    steps = sorted(int(mm.group(1)) for f in glob.glob(os.path.join(rd, "lustre_*.dat")) if (mm := re.search(r"lustre_(\d+)\.dat$", f)))
    nv = np.fromfile(os.path.join(rd, f"lustre_{steps[0]}.dat"), np.float64).size // cells
    x = (np.arange(cells) + 0.5) * (L / cells)
    kappa = 2.0 * math.pi / L
    basis = np.sin(kappa * x)
    norm = (basis**2).sum()
    t, amp, umax = [], [], []
    for s in steps:
        d = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64).reshape(nv, cells)
        rho, mom = d[0], d[1]
        T_s = d[TS_IDX0]  # the scalar, read directly -- no EOS inversion
        a = (T_s * basis).sum() / norm  # projection onto sin mode
        t.append(s * dt)
        amp.append(a)
        umax.append(np.max(np.abs(mom / rho)))
    t, amp = np.array(t), np.array(amp)
    a0 = amp[0]
    good = np.abs(amp) > 1e-12
    rate = -np.polyfit(t[good], np.log(np.abs(amp[good] / a0)), 1)[0]
    analytic = alpha * kappa**2
    print(f"  initial mode amplitude A0 = {a0:.4f}  (eps*T0 = {eps * T0:.4f})")
    print(f"  measured decay rate    = {rate:.4f}")
    print(f"  analytic alpha*kappa^2 = {analytic:.4f}  (alpha={alpha}, kappa=2pi/L)")
    print(f"  ratio measured/analytic = {rate / analytic:.4f}  -> error {100 * abs(rate / analytic - 1):.1f}%")
    print(f"  max |u| over the run    = {max(umax):.3e}  (should stay ~0: density is uniform, T_s is decoupled)")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--measure":
        measure(sys.argv[2])
    else:
        make_case()
