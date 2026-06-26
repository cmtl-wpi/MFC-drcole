#!/usr/bin/env python3
# 3D thermocapillary YGB migration -- STRETCHED-GRID variant for high droplet resolution.
#
# The uniform YGB sweep showed v_t/u_YGB is limited only by grid resolution (cells per drop
# diameter), climbing 0.676 -> 0.757 over cpd=6.4 -> 9.6 and extrapolating to ~0.92 as dx->0.
# Reaching cpd~20 with a UNIFORM cube over W=10D would be ~200^3 = 8M cells. Instead we put a fine
# UNIFORM CORE around the drop and STRETCH outward to an unconfined far field: cpd~20 on the drop at
# ~0.5M cells.
#
# MFC grid stretching is a DOMAIN-EXPANDING sponge (docs/documentation/case.md): the input
# [*_domain%beg,*_domain%end] is the ~uniform fine core (core dx = (end-beg)/(m+1)); cells grow
# outward and the actual boundary lands well beyond the input end. So here the input domain is the
# fine core (+/- YGB_COREHALF in D), and stretching pushes the walls out to ~YGB_WEFF in D.
#
# THE Y-AXIS TRAP: stretching y pushes the isothermal gradient walls from +/-core_half to the ACTUAL
# expanded boundary. The wall temperatures must be pinned to T(y)=T0+gradT*y at those TRUE positions
# (computed below by replicating MFC's stretch transform), else a temperature discontinuity at the
# wall corrupts the imposed field and the Marangoni forcing. The analytic rho(y)/pres/cf ICs need no
# change -- MFC evaluates them on the actual stretched cell centers.
#
# Physics block (EOS, YGB scales, conduction, analytic patch, surface tension) is copied verbatim
# from case_ygb.py (cube geom). Only the grid (stretched) and the wall temps differ.

import json
import math
import os
import sys

# -- Stretched-grid knobs (env) --
cpd = float(os.environ.get("YGB_CPD", "20"))  # core cells per D == drop resolution
core_half = float(os.environ.get("YGB_COREHALF", "2.0"))  # fine-core half-extent in D
W_eff = float(os.environ.get("YGB_WEFF", "10.0"))  # target effective FULL width in D (far field)
frac = float(os.environ.get("YGB_FRAC", "0.75"))  # x_a/x_b as fraction of core_half (uniform band)
loops = int(os.environ.get("YGB_LOOPS", "3"))  # stretch recursion count
Ma = float(os.environ.get("YGB_MA", "0.5"))  # thermal Marangoni number (conduction strength)
n_tr = float(os.environ.get("YGB_TR", "3.0"))  # run length in capillary-thermal times t_r
assert Ma > 0, "3D needs conduction (YGB_MA > 0): a frozen-T 3D rise has no plateau to validate."


def _stretch_cb(beg, end, ncells, a, xa, xb, loops):
    """Replicate src/pre_process/m_grid.f90 grid stretching on one axis. Returns the stretched
    cell-boundary array (length ncells+1), matching MFC's lustre_*_cb.dat. ncells = m+1."""
    m = ncells - 1
    dx = (end - beg) / (m + 1)
    cb = [beg + dx * k for k in range(ncells + 1)]
    cb[-1] = end
    length = abs(cb[-1] - cb[0])
    cb = [v / length for v in cb]
    A, B = xa / length, xb / length
    for _ in range(loops):
        cb = [c / a * (a + math.log(math.cosh(a * (c - A))) + math.log(math.cosh(a * (c - B))) - 2.0 * math.log(math.cosh(a * (B - A) / 2.0))) for c in cb]
    return [v * length for v in cb]


def _expanded_half(beg, end, ncells, a, xa, xb, loops):
    return _stretch_cb(beg, end, ncells, a, xa, xb, loops)[-1]


def _solve_a(beg, end, ncells, xa, xb, loops, target_half, lo=0.2, hi=12.0):
    """Bisection for the stretch rate a giving expanded half-width == target_half."""

    def f(a):
        return _expanded_half(beg, end, ncells, a, xa, xb, loops) - target_half

    if f(lo) * f(hi) > 0:
        raise ValueError(f"W_eff={2 * target_half}D unreachable at loops={loops}, frac={frac}, core_half={core_half}; raise YGB_LOOPS or move x_a/x_b inward (lower YGB_FRAC).")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# -- Geometry: cube, drop centered. Input domain is the fine core; stretching expands it. --
D = 1.0
r = D / 2.0
y_drop = 0.0
dx = D / cpd  # core cell size == drop resolution
beg, end = -core_half, core_half  # input domain (the uniform fine core), in D
m = round(2.0 * core_half / dx) - 1  # core cells - 1; same on all axes (cube)
n = m
p = m
ncells = m + 1
xa = -frac * core_half  # stretch onset (negative); [xa,xb] stays uniform
xb = frac * core_half

# -- Solve the stretch rate for the requested effective width; report the achieved grid. --
a = _solve_a(beg, end, ncells, xa, xb, loops, W_eff / 2.0)
ycb = _stretch_cb(beg, end, ncells, a, xa, xb, loops)  # y identical to x,z (cube)
y_lo, y_hi = ycb[0], ycb[-1]  # ACTUAL expanded boundaries (where the walls really are)
full_eff = y_hi - y_lo
widths = [ycb[i + 1] - ycb[i] for i in range(len(ycb) - 1)]
dx_core_actual, dx_edge = min(widths), max(widths)

# -- Equation of state (two IDENTICAL stiffened-gas fluids, gamma=2; mu*=1, k*=1) -- (copied) --
gam = 2.0
p_inf, p0 = 32.0, 8.0
mu = 0.1
T0 = 10.0
gradT = 2.0 / 15.0
sigma0 = 0.1
sigma_T = -0.1
rho_b = 0.2
rho_coeff = rho_b * T0  # = 2.0
cv_b = (p0 + p_inf) / ((gam - 1.0) * rho_coeff)  # = 20
T_ref = T0 + gradT * y_drop  # sigma = sigma0 at the drop's initial center
eps = 1.0e-9

# -- YGB terminal speed and time scales -- (copied) --
v_YGB = (2.0 / 15.0) * (-sigma_T) * gradT * r / mu  # = 8.889e-3
U_r = (-sigma_T) * gradT * r / mu
t_r = mu / abs(sigma_T * gradT)  # = 7.5

# -- Bulk thermal conduction: alpha_T from the thermal Marangoni number Ma -- (copied) --
alpha_T = U_r * r / Ma
cp_b = gam * cv_b
k_therm = alpha_T * rho_b * cp_b

# -- One analytic patch: smooth sphere eta(x,y,z) drives color, Laplace jump, and density -- (copied) --
xc_d, yc_d, zc_d, r_d = 0.0, y_drop, 0.0, r
w_if = 0.75 * dx  # interface half-width (~3-cell transition, in the fine core)
laplace = sigma0 / r
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2 + (z - ({zc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
cf_expr = f"({eta})"
pres_expr = f"{p0:.9f} + {laplace:.9f}*({eta})"
rho_num = (1.0 - eps) * rho_coeff
rho_expr = f"{rho_num:.9f}*({p0 + p_inf:.9f} + {laplace:.9f}*({eta}))/({p0 + p_inf:.9f}*({T0} + {gradT:.9f}*y))"

# -- Time stepping: min(acoustic CFL, 3D explicit-diffusion). Min cell == core dx; rho_min/c_max use
#    the EXPANDED hot wall (hottest, lowest-density cell after stretching). --
rho_min = rho_coeff / (T0 + gradT * y_hi)
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5
mydt = 0.35 * dx / c_max
mydt = min(mydt, 0.35 * dx**2 / (6.0 * alpha_T))
t_step_stop = int(round(n_tr * t_r / mydt))
t_step_save = max(1, t_step_stop // 60)

print(
    f"[case_ygb_stretch] cpd={cpd:g} core=+/-{core_half:g}D m=n=p={m} dx_core={dx:.5g} "
    f"a={a:.5g} loops={loops} -> W_eff={full_eff:.4g}D (+/-{y_hi:.4g}) edge/core={dx_edge / dx_core_actual:.2f} "
    f"cells={ncells**3} Twall_in/out={T0 + gradT * y_lo:.6g}/{T0 + gradT * y_hi:.6g} "
    f"dt={mydt:.4g} steps={t_step_stop}",
    file=sys.stderr,
)

data = {
    "run_time_info": "T",
    # Input domain = the fine CORE (+/- core_half); stretching expands the actual walls to +/- y_hi.
    "x_domain%beg": beg,
    "x_domain%end": end,
    "y_domain%beg": beg,
    "y_domain%end": end,
    "z_domain%beg": beg,
    "z_domain%end": end,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
    # Grid stretching on all three axes: uniform fine core in [xa,xb], coarsening outward.
    "stretch_x": "T",
    "a_x": a,
    "x_a": xa,
    "x_b": xb,
    "loops_x": loops,
    "stretch_y": "T",
    "a_y": a,
    "y_a": xa,
    "y_b": xb,
    "loops_y": loops,
    "stretch_z": "T",
    "a_z": a,
    "z_a": xa,
    "z_b": xb,
    "loops_z": loops,
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # Simulation algorithm (6-equation model; proven WENO5/HLLC settings)
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "F",
    "time_stepper": 3,
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Boundaries: slip walls (-2) on all six faces; isothermal gradient walls on y set below
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "bc_z%beg": -2,
    "bc_z%end": -2,
    "num_patches": 1,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + temperature-dependent surface tension sigma(T)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    # Continuous phase (fluid 1)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    "fluid_pp(1)%k_therm": k_therm,
    # Second phase (fluid 2) -- identical properties (mu* = 1, k* = 1)
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_b,
    "fluid_pp(2)%Re(1)": 1.0 / mu,
    "fluid_pp(2)%k_therm": k_therm,
    # Single analytic patch (3D cuboid spanning the EXPANDED domain, so every stretched cell gets the
    # analytic field). The drop lives entirely in eta(x,y,z): cf, pres, and alpha_rho all share it.
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": full_eff,
    "patch_icpp(1)%length_y": full_eff,
    "patch_icpp(1)%length_z": full_eff,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": rho_expr,
    "patch_icpp(1)%alpha_rho(2)": eps,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": cf_expr,
    # Isothermal Dirichlet gradient walls -- pinned to T(y) at the ACTUAL expanded boundaries.
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T0 + gradT * y_lo,  # cold floor at the expanded -y wall
    "bc_y%Twall_out": T0 + gradT * y_hi,  # hot ceiling at the expanded +y wall
}

print(json.dumps(data))
