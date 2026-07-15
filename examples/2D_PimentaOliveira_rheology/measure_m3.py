#!/usr/bin/env python3
# Pimenta & Oliveira (2021) rheology diagnostic (guide Table 3.3). Bulk interfacial contribution to the
# suspension stress = volume average of the CSF capillary stress tensor T_ij = sigma*(delta_ij - n_i n_j)*|grad c|,
# with n = grad c/|grad c| and sigma(Gamma) the Langmuir EOS (Gamma = surf/|grad c|, matching the solver).
# Decompose sigma = <sigma>_band + sigma' : the mean gives the CAPILLARY stress [eta_c], the fluctuation the
# MARANGONI stress [eta_m]; by linearity [eta] = [eta_c] + [eta_m] exactly. Reports intrinsic viscosities
# [eta*] = Sigma_xy/(mu*gdot*phi) and first normal-stress difference N1 = Sigma_xx - Sigma_yy. Gates:
# clean drop (X=0) -> [eta_m]~0; surfactant -> [eta_m]>0 grows with coverage; N1>0. Argv[1]=case dir.
import glob
import json
import os
import re
import sys

import numpy as np

C = sys.argv[1]


def inp(n):
    m = re.search(rf"^\s*{re.escape(n)}\s*=\s*([^\s,]+)", open(os.path.join(C, "simulation.inp")).read(), re.M)
    return float(m.group(1).replace("d", "e").replace("D", "e"))


nx, ny = int(inp("m")) + 1, int(inp("n")) + 1
x0, y0 = inp("x_domain%beg"), inp("y_domain%beg")
Wd, Hd = inp("x_domain%end") - x0, inp("y_domain%end") - y0
dx, dy = Wd / nx, Hd / ny
mu = 1.0 / inp("fluid_pp(1)%Re(1)")  # matrix viscosity
gdot = inp("bc_y%ve1") / (Hd / 2)  # shear rate from moving-wall velocity
sig0, E_el, surf_max = inp("sigma"), inp("sigma_El"), inp("surf_max")
R = 1.0
phi = np.pi * R**2 / (Wd * Hd)  # drop area fraction

X = np.tile(x0 + (np.arange(nx) + 0.5) * dx, ny).reshape(ny, nx)
Y = np.repeat(y0 + (np.arange(ny) + 0.5) * dy, nx).reshape(ny, nx)
fs = sorted(
    [f for f in glob.glob(C + "/restart_data/lustre_*.dat") if re.fullmatch(r"\d+", os.path.basename(f)[7:-4])],
    key=lambda f: int(os.path.basename(f)[7:-4]),
)
ss = os.path.getsize(fs[0]) // 8 // (nx * ny)
cutoff = 1e-6


def frame(f):
    a = np.fromfile(f, "<f8").reshape(ss, ny, nx)
    c, surf = a[ss - 2], a[ss - 1]
    # deformation (inertia tensor)
    w = np.clip(c, 0, 1)
    wn = w / w.sum()
    mx, my = (X * wn).sum(), (Y * wn).sum()
    Ixx = ((X - mx) ** 2 * wn).sum()
    Iyy = ((Y - my) ** 2 * wn).sum()
    Ixy = ((X - mx) * (Y - my) * wn).sum()
    tr = Ixx + Iyy
    disc = np.sqrt(max(tr * tr / 4 - (Ixx * Iyy - Ixy**2), 0.0))
    aa, bb = np.sqrt(tr / 2 + disc), np.sqrt(max(tr / 2 - disc, 0.0))
    D = (aa - bb) / (aa + bb) if (aa + bb) > 0 else 0.0
    theta = 0.5 * np.degrees(np.arctan2(2 * Ixy, Ixx - Iyy))
    # CSF capillary stress tensor, volume-averaged over the interface band
    gx = np.gradient(c, dx, axis=1)
    gy = np.gradient(c, dy, axis=0)
    normc = np.sqrt(gx**2 + gy**2)
    band = normc > cutoff
    gam = np.where(band, surf / np.where(band, normc, 1.0), 0.0)  # Gamma = surf/|grad c|
    arg = np.clip(1.0 - gam / surf_max, cutoff, None)
    sig = sig0 * (1.0 + E_el * np.log(arg))  # Langmuir sigma(Gamma), matches solver
    sig = np.maximum(sig, 1e-3 * sig0)  # solver floors c_sigma at 1e-3*sigma
    sig = np.where(band, sig, 0.0)
    # <sigma> = |grad c|-weighted band mean; split into mean (capillary) + fluctuation (Marangoni)
    wsum = normc[band].sum()
    sig_mean = (sig[band] * normc[band]).sum() / wsum if wsum > 0 else 0.0
    fac = dx * dy / (Wd * Hd)
    nxn = np.where(band, gx / np.where(band, normc, 1.0), 0.0)
    nyn = np.where(band, gy / np.where(band, normc, 1.0), 0.0)

    # T_ij = sigma*(delta_ij - n_i n_j)*|grad c|;  Sigma_ij = fac * sum(T_ij)
    def stress(sfield):
        Sxy = fac * np.sum(sfield * (-nxn * nyn) * normc)
        Sxx = fac * np.sum(sfield * (1.0 - nxn * nxn) * normc)
        Syy = fac * np.sum(sfield * (1.0 - nyn * nyn) * normc)
        return Sxy, Sxx, Syy

    Sxy, Sxx, Syy = stress(sig)  # total
    Cxy, Cxx, Cyy = stress(np.where(band, sig_mean, 0.0))  # capillary (mean sigma)
    norm = mu * gdot * phi
    eta = Sxy / norm
    eta_c = Cxy / norm
    eta_m = eta - eta_c
    N1 = (Sxx - Syy) / norm
    return D, theta, eta, eta_c, eta_m, N1, sig_mean, surf.sum()


# Select the latest window where surfactant mass is still conserved (<2% drift) -- excludes any late-time
# tip-instability onset (gentle Ca=0.1 rarely triggers it, but keep the guard consistent with M2).
allf = [(f, frame(f)) for f in fs]
m0 = allf[0][1][7]
valid = [(f, v) for f, v in allf if m0 <= 0 or v[7] / m0 < 1.02]  # m0=0 at X=0 (no surfactant)
use = valid[-4:] if len(valid) >= 4 else (valid or allf[-4:])
v = np.array([vv for _, vv in use]).mean(axis=0)
D, theta, eta, eta_c, eta_m, N1, sig_mean, mass = v
print(
    json.dumps(
        {
            "D": round(float(D), 5),
            "theta_deg": round(float(theta), 2),
            "eta_intrinsic": round(float(eta), 5),
            "eta_capillary": round(float(eta_c), 5),
            "eta_marangoni": round(float(eta_m), 5),
            "N1_star": round(float(N1), 5),
            "sigma_mean_band": round(float(sig_mean), 5),
            "phi": round(float(phi), 5),
        }
    )
)
