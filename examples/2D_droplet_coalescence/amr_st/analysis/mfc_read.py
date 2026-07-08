"""Shared MFC lustre reader for the AMR+ST experiment.

Reads a run's simulation.inp and raw restart_data/lustre_<step>.dat checkpoints.
Field layout is fixed by the model (model_eqns=2, num_fluids=2, surface_tension=T,
2D) -- verified against an archived case-f run:

    0: alpha_rho_1 (liquid)   4: E (total energy)
    1: alpha_rho_2 (gas)      5: alpha_1  (liquid volume fraction)
    2: rho*u                  6: alpha_2  (gas volume fraction)
    3: rho*v                  7: cf       (color function, drives surface tension)

Flat float64, C-reshape to (nvars, n+1, m+1); x is the innermost (fastest) axis.
For an AMR run the lustre_<step>.dat holds the level-0 (coarse) field -- exactly the
grid on which the coarse/fine seam current would appear.

Never hardcode grid/EOS: everything comes from simulation.inp (RESEARCH_WORKFLOWS s7).
"""
import glob
import os
import re

import numpy as np

NVARS = 8
IDX = dict(arho1=0, arho2=1, mom_x=2, mom_y=3, E=4, alpha1=5, alpha2=6, cf=7)


def read_inp(run_dir, fname="simulation.inp"):
    """Parse an MFC .inp (Fortran namelist) into a dict of python scalars.

    simulation.inp holds grid/EOS/algorithm/AMR; pre_process.inp additionally holds
    the patch geometry/density (radius, alpha_rho) that simulation.inp omits.
    """
    txt = open(os.path.join(run_dir, fname)).read()
    d = {}
    for m in re.finditer(r"^\s*([A-Za-z0-9_%()]+)\s*=\s*([^\n]+?)\s*$", txt, re.M):
        k, v = m.group(1), m.group(2).rstrip("/").strip()
        if v in ("T", "F"):
            d[k] = (v == "T")
        else:
            try:
                d[k] = int(v) if re.fullmatch(r"[-+]?\d+", v) else float(v)
            except ValueError:
                d[k] = v
    return d


def grid(run_dir):
    """Cell-center coordinates (x, y) from the lustre cell-boundary files."""
    rd = os.path.join(run_dir, "restart_data")
    xcb = np.fromfile(os.path.join(rd, "lustre_x_cb.dat"), dtype=np.float64)
    ycb = np.fromfile(os.path.join(rd, "lustre_y_cb.dat"), dtype=np.float64)
    x = 0.5 * (xcb[:-1] + xcb[1:])
    y = 0.5 * (ycb[:-1] + ycb[1:])
    return x, y


def list_steps(run_dir):
    rd = os.path.join(run_dir, "restart_data")
    steps = []
    for f in glob.glob(os.path.join(rd, "lustre_*.dat")):
        m = re.search(r"lustre_(\d+)\.dat$", f)
        if m:
            steps.append(int(m.group(1)))
    return sorted(steps)


def read_step(run_dir, step, inp, x, y):
    """Return primitive fields for one checkpoint as a dict of (ny, nx) arrays."""
    m, n = inp["m"], inp["n"]
    nx, ny = m + 1, n + 1
    f = os.path.join(run_dir, "restart_data", f"lustre_{step}.dat")
    q = np.fromfile(f, dtype=np.float64).reshape(NVARS, ny, nx)

    arho1, arho2 = q[IDX["arho1"]], q[IDX["arho2"]]
    rho = arho1 + arho2
    u = q[IDX["mom_x"]] / rho
    v = q[IDX["mom_y"]] / rho
    speed = np.sqrt(u * u + v * v)
    alpha1 = q[IDX["alpha1"]]
    cf = q[IDX["cf"]]

    # mixture stiffened-gas pressure: MFC stores gamma_mfc=1/(g-1), pi_inf_mfc=g*pi/(g-1)
    g1, g2 = inp["fluid_pp(1)%gamma"], inp["fluid_pp(2)%gamma"]
    pi1, pi2 = inp["fluid_pp(1)%pi_inf"], inp["fluid_pp(2)%pi_inf"]
    alpha2 = q[IDX["alpha2"]]
    gamma_mix = alpha1 * g1 + alpha2 * g2
    pi_mix = alpha1 * pi1 + alpha2 * pi2
    pres = (q[IDX["E"]] - 0.5 * rho * speed**2 - pi_mix) / gamma_mix

    return dict(rho=rho, u=u, v=v, speed=speed, alpha1=alpha1, cf=cf, pres=pres)


def buff_size(inp):
    """Replicate s_configure_coordinate_bounds for the seam-band half-width."""
    if inp.get("recon_type", 1) == 1:            # WENO
        weno_polyn = (inp.get("weno_order", 5) - 1) // 2
        return 2 * weno_polyn + 2 if inp.get("viscous", False) else weno_polyn + 2
    return inp.get("muscl_order", 2) + 2          # MUSCL


def amr_block(inp):
    """(beg, end) level-0 cell indices of the static block, or None if not an AMR run."""
    if not inp.get("amr", False):
        return None
    return (int(inp["amr_block_beg(1)"]), int(inp["amr_block_end(1)"]),
            int(inp["amr_block_beg(2)"]), int(inp["amr_block_end(2)"]))
