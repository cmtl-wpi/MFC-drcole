#!/usr/bin/env python3
# Quantitative analytic verification of the bulk-conduction flux: decay of a sinusoidal thermal
# mode on a PERIODIC 1D domain (no temperature-BC needed -- the cleanest possible conduction test).
#
# A small temperature perturbation T(x,0) = T0*(1 + eps*sin(2*pi*x/L)), imposed through the density
# at uniform pressure, decays under conduction. In the low-Mach (isobaric) limit the fundamental
# mode amplitude obeys A(t) = A(0)*exp(-alpha*kappa^2*t) with kappa = 2*pi/L and alpha = k/(rho*cp),
# cp = gamma*cv. Projecting the temperature field onto sin(2*pi*x/L) at each output and fitting
# ln(A) vs t recovers the decay rate; it must match alpha*kappa^2 to a few percent. Acoustic
# oscillations (period ~ L/c0, far faster than the conduction e-fold) average out of the slow
# envelope. This isolates the conduction OPERATOR from the droplet/Marangoni coupling.
#
# Usage:
#   ./mfc.sh run examples/2D_thermocapillary_migration/verify_1d_diffusion.py -n 1 -t pre_process simulation
#   python3 examples/2D_thermocapillary_migration/verify_1d_diffusion.py --measure <rundir>
import json
import sys

# Shared physical constants (module level so both the generator and the measurer agree)
L = 1.0
m = 199
gam, p_inf, cv, p0, T0 = 2.0, 100.0, 12.5, 25.0, 10.0
eps = 0.03  # perturbation amplitude (small -> low-Mach, isobaric)
alpha = 0.05  # target thermal diffusivity
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * T0)  # = 1.0
k_therm = alpha * rho0 * cv * gam  # alpha = k/(rho*cp), cp = gam*cv


def make_case():
    import math

    kappa = 2.0 * math.pi / L
    c0 = (gam * (p0 + p_inf)) ** 0.5
    dx = L / (m + 1)
    dt = min(0.4 * dx / c0, 0.35 * dx**2 / (2.0 * alpha))
    t_end = 0.6 / (alpha * kappa**2)  # ~0.6 e-folding times of decay
    t_step_stop = int(round(t_end / dt))
    # density realizing T(x) = T0*(1+eps*sin(kappa x)) at uniform p0:
    #   rho(x) = (p0+p_inf)/((gam-1)*cv*T(x))  ->  rho0 / (1 + eps*sin(kappa x))
    rho_expr = f"{rho0:.9f}/(1.0 + {eps}*sin({kappa:.9f}*x))"
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
        "patch_icpp(1)%alpha_rho(1)": rho_expr,
        "patch_icpp(1)%alpha(1)": 1.0,
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
    t, amp = [], []
    for s in steps:
        d = np.fromfile(os.path.join(rd, f"lustre_{s}.dat"), np.float64).reshape(nv, cells)
        rho, mom, E = d[0], d[1], d[2]
        u = mom / rho
        e = E / rho - 0.5 * u**2
        pres = (gam - 1.0) * rho * e - gam * p_inf  # stiffened-gas pressure
        T = (pres + p_inf) / ((gam - 1.0) * rho * cv)
        a = (T * basis).sum() / norm  # projection onto sin mode
        t.append(s * dt)
        amp.append(a)
    t, amp = np.array(t), np.array(amp)
    a0 = amp[0]
    # fit ln|A/A0| = -rate * t over the decaying window
    good = np.abs(amp) > 1e-12
    rate = -np.polyfit(t[good], np.log(np.abs(amp[good] / a0)), 1)[0]
    analytic = alpha * kappa**2
    print(f"  initial mode amplitude A0 = {a0:.4f}  (eps*T0 = {eps * T0:.4f})")
    print(f"  measured decay rate  = {rate:.4f}")
    print(f"  analytic alpha*kappa^2 = {analytic:.4f}  (alpha={alpha}, kappa=2pi/L)")
    print(f"  ratio measured/analytic = {rate / analytic:.4f}  -> error {100 * abs(rate / analytic - 1):.1f}%")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--measure":
        measure(sys.argv[2])
    else:
        make_case()
