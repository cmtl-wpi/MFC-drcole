#!/usr/bin/env python3
"""
Phase 1 -- static-drop Laplace seam test (AMR + surface tension containment experiment).

One static liquid drop (radius R) in gas, pressure initialized with the 2D Laplace
jump so the IC is in mechanical equilibrium. Any residual velocity is a parasitic
(spurious) current from the surface-tension discretization. We compare three runs at
matched physical time:

    coarse : uniform R/50   (400x400)          baseline parasitic current
    fine   : uniform R/100  (800x800)          what R/100 resolution produces
    amr    : R/50 base + static 2:1 block      R/100 near interface, coarse elsewhere

The AMR block (+-1.75R) sits 0.75R (~37 coarse cells) OUTSIDE the interface (at R), so
c is exactly {0,1} across the coarse/fine seam and |grad c| < capillary_cutoff there --
the inconsistent-normal seam mechanism (PR 1628 gate) has nothing to act on. The test:
does the contained AMR run stay at the parasitic baseline (PASS) or grow a 27-540x seam
current (FAIL)?

Physics (EOS, densities, sigma, sound speed) are taken verbatim from case f of
../../case.py (the cavitation-fix c_l=100 minimizes step count and keeps Ma~0.005).
Numerics deviate from case.py in two required/deliberate ways:
  * recon_type: WENO5 (AMR requires WENO; case.py uses MUSCL). ST-safe WENO settings
    are copied from the shipped examples/2D_laplace_pressure_jump case.
  * Laplace jump: sigma/R (correct 2D value; case.py's 2*sigma/R is the 3D-sphere
    formula and would leave a static drop unbalanced).

Usage:  case_laplace.py --variant {coarse,fine,amr} [--n-periods N] [--n-start K]
"""

import argparse
import json
import math
import sys

p = argparse.ArgumentParser()
p.add_argument("--variant", choices=["coarse", "fine", "amr"], required=True)
p.add_argument("--n-periods", type=float, default=10.0,
               help="run length in capillary periods tau = sqrt(rho_l R^3 / sigma)")
p.add_argument("--n-saves", type=int, default=200, help="number of checkpoint saves")
p.add_argument("--n-start", type=int, default=0, help="restart checkpoint index (cfl_dt)")
p.add_argument("--amr-max-blocks", type=int, default=1, help="AMR block slots")
p.add_argument("--no-st", action="store_true", help="disable surface_tension (bug isolation)")
# feasibility knobs (defaults = handoff spec). Lowering c_l cuts step count (dt ~ 1/c);
# the seam mechanism is sound-speed-independent, so this is safe for the mechanism test.
p.add_argument("--c-l", type=float, default=100.0, help="liquid sound speed [m/s]")
p.add_argument("--domain-R", type=float, default=4.0, help="domain half-width in units of R")
p.add_argument("--block-R", type=float, default=1.75, help="AMR block half-width in units of R")
# MFC passes --mfc '{...}' when it execs the case; accept and ignore it.
p.add_argument("--mfc", type=json.loads, default={})
args = p.parse_args()

# ---- physics from case f (../../case.py) --------------------------------------
rho_l = 763.0
rho_g = rho_l / 666.0
sigma = 0.0266
gamma_l = 3.7
gamma_g = 1.4
R_um = 159.0                      # case f radius
R = R_um * 1e-6
c_l_phys = args.c_l               # 100 = cavitation-fix speed; lowered for np=1 feasibility

p_gas = c_l_phys**2 * rho_g / gamma_g            # c_g == c_l (no acoustic dt penalty)
pi_inf_l = rho_l * c_l_phys**2 / gamma_l - p_gas
pi_inf_g = 0.0
p_liq = p_gas + sigma / R                        # 2D Laplace jump (one curvature)

# viscosity from case f's We/Re (absolute mu, independent of the static setup)
We_f, Re_f = 32.8, 210.8
D = 2.0 * R
Ur_f = math.sqrt(We_f * sigma / (rho_l * D))
mu_l = rho_l * Ur_f * D / Re_f
mu_g = mu_l / 119.0

# ---- grid ---------------------------------------------------------------------
# coarse dx = R/50 over the +-domain_R box; fine = R/100 (2x cells).
half = args.domain_R * R
ncell_coarse = round(2.0 * args.domain_R * 50)     # cells at R/50
ncell = 2 * ncell_coarse if args.variant == "fine" else ncell_coarse
Nx = Ny = ncell - 1
dx = 2.0 * half / ncell          # == R/50 (coarse) or R/100 (fine)

# ---- capillary time & schedule ------------------------------------------------
tau = math.sqrt(rho_l * R**3 / sigma)
t_stop = args.n_periods * tau
t_save = t_stop / args.n_saves

eps = 1e-9

data = {
    "run_time_info": "T",
    "x_domain%beg": -half,
    "x_domain%end": half,
    "y_domain%beg": -half,
    "y_domain%end": half,
    "m": Nx,
    "n": Ny,
    "p": 0,
    "cyl_coord": "F",
    "n_start": args.n_start,
    "cfl_adap_dt": "T",
    "cfl_target": 0.1,
    "t_stop": t_stop,
    "t_save": t_save,
    # ---- algorithm (5-eq, WENO5, RK3, HLLC) -----------------------------------
    "model_eqns": 2,
    "num_fluids": 2,
    "mpp_lim": "T",              # required by AMR for num_fluids>1
    "mixture_err": "T",
    "alt_soundspeed": "F",
    "low_Mach": 0,
    "time_stepper": 3,
    "recon_type": 1,            # WENO (AMR requirement)
    "weno_order": 5,
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "weno_Re_flux": "F",
    "weno_avg": "T",
    "avg_state": 2,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "viscous": "T",
    "surface_tension": "F" if args.no_st else "T",
    "sigma": sigma,
    # ---- non-reflecting outflow on all sides (as case f) ----------------------
    "bc_x%beg": -8, "bc_x%end": -8,
    "bc_y%beg": -8, "bc_y%end": -8,
    # ---- output ---------------------------------------------------------------
    "format": 1,
    "precision": 2,
    "parallel_io": "T",
    "alpha_wrt(1)": "T",
    "cf_wrt": "T",
    "vel_wrt(1)": "T",
    "vel_wrt(2)": "T",
    "pres_wrt": "T",
    # ---- fluids (case f EOS + viscosity) --------------------------------------
    "fluid_pp(1)%gamma": 1.0 / (gamma_l - 1.0),
    "fluid_pp(1)%pi_inf": gamma_l * pi_inf_l / (gamma_l - 1.0),
    "fluid_pp(1)%Re(1)": 1.0 / mu_l,
    "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "fluid_pp(2)%Re(1)": 1.0 / mu_g,
    # ---- patches: gas background + liquid drop --------------------------------
    "num_patches": 2,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": 32.0 * R,      # oversize; clipped to grid
    "patch_icpp(1)%length_y": 32.0 * R,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p_gas,
    "patch_icpp(1)%alpha_rho(1)": eps * rho_l,
    "patch_icpp(1)%alpha_rho(2)": (1.0 - eps) * rho_g,
    "patch_icpp(1)%alpha(1)": eps,
    "patch_icpp(1)%alpha(2)": 1.0 - eps,
    "patch_icpp(1)%cf_val": 0,
    "patch_icpp(2)%geometry": 2,
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%radius": R,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": 0.99,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%pres": p_liq,
    "patch_icpp(2)%alpha_rho(1)": (1.0 - eps) * rho_l,
    "patch_icpp(2)%alpha_rho(2)": eps * rho_g,
    "patch_icpp(2)%alpha(1)": 1.0 - eps,
    "patch_icpp(2)%alpha(2)": eps,
    "patch_icpp(2)%cf_val": 1,
}

# ---- static 2:1 AMR block, contained around the interface ---------------------
if args.variant == "amr":
    # block spans +-block_R*R -> level-0 cell indices; interface at R stays
    # (block_R-1)*50 coarse cells inside, so the c/f seam sees c in {0,1} exactly.
    beg = int(math.floor((-args.block_R * R + half) / dx))
    end = int(math.ceil((args.block_R * R + half) / dx))
    width = end - beg + 1
    assert 2 * width - 1 <= Nx, f"block width {width} violates 2w-1<=m_glb ({Nx})"
    data.update({
        "amr": "T",
        "amr_block_beg(1)": beg,
        "amr_block_beg(2)": beg,
        "amr_block_end(1)": end,
        "amr_block_end(2)": end,
        "amr_regrid_int": 0,     # static
        "amr_max_blocks": args.amr_max_blocks,
    })
    print(f"[case_laplace] amr block cells {beg}-{end} (w={width}); "
          f"interface at cell {int(round((R+half)/dx))}", file=sys.stderr)

if args.no_st:                     # bug-isolation: strip ST-only knobs
    for k in ("cf_wrt", "sigma"):
        data.pop(k, None)

print(f"[case_laplace] variant={args.variant} ncell={ncell} dx=R/{R/dx:.0f} "
      f"tau={tau:.4e}s t_stop={t_stop:.4e}s ({args.n_periods} periods) "
      f"p_gas={p_gas:.1f} p_liq={p_liq:.1f} dP={sigma/R:.1f}Pa", file=sys.stderr)

print(json.dumps(data))
