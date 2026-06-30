#!/usr/bin/env python3
# Parametrized 2D thermocapillary migration of a single drop/bubble, for the terminal-velocity TREND
# sweeps (sweep.py) that validate the single-particle claims in the Introduction of Nas & Tryggvason,
# Int. J. Multiphase Flow 29 (2003) 1117-1135 (a summary of prior single-particle literature, not a
# figure in that paper):
#   * gas bubble V_t: decreases rapidly with Ca, increases very weakly with Re, decreases with Ma
#   * single drop V_t: decreases with Ma, reaches a minimum, then increases (U-shape)
#   * drops deform oblate/prolate depending on the density ratio
#
# This is the case_Ma_20.py construction (distinct-fluid temperature via a per-fluid 1/T density proxy
# + bulk conduction + sigma(T); see that file's header for the full rationale) generalized to (a)
# arbitrary Re/Ma/Ca and (b) INDEPENDENT property ratios so a low-ratio gas bubble and an O(1) drop
# share one case. A "bubble" is all ratios ~1/25; a "drop" is ratios ~0.5. The box is 2D wide x 4D
# tall (Nas-Tryggvason Fig 3 geometry) with the particle 1D above the cold floor, giving vertical room
# to reach a quasi-steady terminal velocity before it nears the hot wall. Geometry per their Fig 1:
# periodic in x, no-slip isothermal walls on the gradient axis (+y).
#
# The block below is rewritten line-by-line by sweep.py (regex on `^name = ...`); the committed
# defaults reproduce the Fig 3 point (Re=5, Ma=20, Ca=0.01666, ratios 0.5). Keep `e`-notation out of
# the analytic strings (MFC's IC parser substitutes bare e/r/eps).

import json

# --- sweep knobs (rewritten by sweep.py) -------------------------------------------------------
Re = 5.0
Ma = 20.0
Ca = 0.01666
rho_ratio = 0.5  # drop / bulk density ratio
mu_ratio = 0.5  # drop / bulk dynamic-viscosity ratio
cv_ratio = 0.5  # drop / bulk heat-capacity ratio
k_ratio = 0.5  # drop / bulk conductivity ratio
Nx = 32  # cells across the box width (16 cells/D at Wx = 2D)
n_tau = 1.5  # run length in VISCOUS times tau = rho_b*r^2/mu_b. The thermocapillary overshoot peak
# (the trend metric) occurs at ~1 tau independent of Re, whereas t_r = mu_b/G scales as 1/Re; pinning
# the run to tau keeps the step count bounded and the drop's migration inside the box across the sweep.
# ----------------------------------------------------------------------------------------------

# Geometry: D = 1 particle 1D above the cold floor of a 2D-wide x 4D-tall box; gradient axis = y
D = 1.0
r = D / 2.0
Wx = 2.0 * D
Hy = 4.0 * D
y_bottom = -Hy / 2.0
y_drop = y_bottom + 1.0 * D

dx = Wx / Nx
Ny = round(Hy / dx)

# Bulk reference scales (arbitrary; only Re, Ma, Ca and the ratios are physical)
rho_b = 1.0
mu_b = 0.02
gam = 2.0

# Invert the non-dimensional numbers for the physical surface-tension and conduction properties
G = Re * mu_b**2 / (rho_b * r**2)  # |sigma_T*gradT|
gradT = 1.0 / Hy
sigma_T = -G / gradT  # dsigma/dT < 0
sigma0 = G * r / Ca  # surface tension at T_ref
alpha_b = G * r**2 / (mu_b * Ma)  # bulk thermal diffusivity (definition of Ma)

T_base = 1.0  # T(y) = T_base + gradT*(y - y_bottom); offset keeps rho > 0


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # sigma = sigma0 at the particle's initial position

# Stiffened-gas EOS, bulk fluid. Two competing requirements set the bulk stiffening p_inf_b:
#   (i)  low migration Mach -> a baseline sound speed c >= ~15  (p0 + p_inf_b >= gam/2 * 15^2 * rho_b);
#   (ii) p_inf_d >= 0 for the PARTICLE. With the per-fluid 1/T proxy, p0 + p_inf_d closes to
#        rho_ratio*cv_ratio*(p0 + p_inf_b), so a low-ratio gas bubble collapses p_inf_d negative unless
#        p0 + p_inf_b is large enough. We size p0 + p_inf_b to satisfy BOTH, then derive c_ref. p0 is
#        floored well above the Laplace jump sigma0/r so the pressure stays positive (the largest jump
#        is at high Re / low Ca).
laplace = sigma0 / r
p0 = max(0.5, 5.0 * laplace)
margin = 1.05  # 5% headroom so p_inf_d stays strictly positive
pb_total = max(gam / 2.0 * 15.0**2 * rho_b, p0 * margin / (rho_ratio * cv_ratio))  # = p0 + p_inf_b
p_inf_b = pb_total - p0
c_ref = (gam * pb_total / rho_b) ** 0.5
cv_b = pb_total / ((gam - 1.0) * rho_b * T_ref)
cp_b = gam * cv_b
k_b = alpha_b * rho_b * cp_b

# Particle = (per-property ratio) * bulk. The EOS stiffening p_inf_d closes the drop density to
# rho_d at (p0, T_ref); the per-fluid 1/T proxy then keeps the density ratio height-independent.
rho_d = rho_ratio * rho_b
mu_d = mu_ratio * mu_b
cv_d = cv_ratio * cv_b
k_d = k_ratio * k_b
cp_d = gam * cv_d
p_inf_d = (gam - 1.0) * rho_d * cv_d * T_ref - p0  # = rho_ratio*cv_ratio*pb_total - p0 >= 0 by construction

U_r = G * r / mu_b  # Marangoni velocity scale
t_r = mu_b / G  # capillary-thermal time
tau = rho_b * r**2 / mu_b  # viscous time (fixed across the sweep; sets the run length)

# Time stepping: min(acoustic CFL, explicit diffusion limit), each taken over BOTH fluids. Density is
# lowest (sound speed highest) at the hot wall; the per-fluid 1/T proxy gives rho_i_hot = ratio*T_ref/T_hot.
T_hot = T_of_y(Hy / 2.0)
rho_b_hot = rho_b * T_ref / T_hot
rho_d_hot = rho_d * T_ref / T_hot
c_b = (gam * (p0 + p_inf_b) / rho_b_hot) ** 0.5
c_d = (gam * (p0 + p_inf_d) / rho_d_hot) ** 0.5
c_max = max(c_b, c_d)
diff_max = max(alpha_b, k_d / (rho_d * cp_d), mu_b / rho_b, mu_d / rho_d)  # thermal + viscous diffusivities
mydt = min(0.35 * dx / c_max, 0.35 * dx**2 / (4.0 * diff_max))
t_step_stop = round(n_tau * tau / mydt)
t_step_save = max(1, t_step_stop // 120)

# One analytic patch over the whole box: the smooth circle eta(x,y) ~ 1 inside / 0 outside drives the
# color, the volume-fraction split, the Laplace pressure jump, and BOTH per-fluid densities together,
# so the mixture EOS recovers the linear T(y) exactly (particle included) while the fluid-to-fluid
# density ratio is preserved. Hardcode the center/radius as decimals.
xc_d, yc_d, r_d = 0.0, y_drop, r
w_if = 0.75 * dx
dist = f"sqrt((x - ({xc_d:.9f}))**2 + (y - ({yc_d:.9f}))**2)"
eta = f"0.5*(1.0 - tanh(({dist} - {r_d:.9f})/{w_if:.9f}))"
Texpr = f"({T_base:.9f} + {gradT:.9f}*(y - ({y_bottom:.9f})))"
pres_expr = f"({p0:.9f} + {laplace:.9f}*({eta}))"
alpha2_expr = f"({eta})"
alpha1_expr = f"(1.0 - ({eta}))"
cb = (gam - 1.0) * cv_b
cd = (gam - 1.0) * cv_d
arho1_expr = f"(1.0 - ({eta}))*({pres_expr} + {p_inf_b:.9f})/({cb:.9f}*{Texpr})"
arho2_expr = f"({eta})*({pres_expr} + {p_inf_d:.9f})/({cd:.9f}*{Texpr})"
cf_expr = f"({eta})"

data = {
    "run_time_info": "T",
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Hy / 2,
    "y_domain%end": Hy / 2,
    "m": Nx - 1,
    "n": Ny - 1,
    "p": 0,
    "cyl_coord": "F",
    "dt": mydt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    "model_eqns": 3,
    "alt_soundspeed": "F",
    "mixture_err": "T",
    "mpp_lim": "T",
    "time_stepper": 3,
    "weno_order": 5,  # weno5 needs >=25 cells/dir; the 32-cell width (16/D) cannot fit weno7
    "weno_eps": 1e-16,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",  # monotonicity-preserving: required for stability at 16/D (3-cell interface + sigma)
    "weno_avg": "T",
    "weno_Re_flux": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # Slip-wall box (Samareh / case_Ma_20 proven config; isothermal gradient walls set below). No-slip
    # (-16) goes unstable at finite Re on this coarse 16/D grid -- fine for the creeping case_NT_fig2.py
    # but not here -- so the sweep uses reflective walls, the validated finite-Re geometry. Confinement
    # is identical across each sweep, so the trend (the quantity of interest) is unaffected.
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 1,
    "num_fluids": 2,
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_conduction": "T",
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "parallel_io": "T",
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf_b / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    "fluid_pp(1)%k_therm": k_b,
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf_d / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    "fluid_pp(2)%k_therm": k_d,
    "patch_icpp(1)%geometry": 3,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Hy,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": pres_expr,
    "patch_icpp(1)%alpha_rho(1)": arho1_expr,
    "patch_icpp(1)%alpha_rho(2)": arho2_expr,
    "patch_icpp(1)%alpha(1)": alpha1_expr,
    "patch_icpp(1)%alpha(2)": alpha2_expr,
    "patch_icpp(1)%cf_val": cf_expr,
    "bc_y%isothermal_in": "T",
    "bc_y%isothermal_out": "T",
    "bc_y%Twall_in": T_of_y(-Hy / 2),
    "bc_y%Twall_out": T_of_y(Hy / 2),
}

print(json.dumps(data))
