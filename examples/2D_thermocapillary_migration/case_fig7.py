#!/usr/bin/env python3
# Thermocapillary droplet migration -- MFC validation against Samareh, Mostaghimi & Moreau,
# "Thermocapillary migration of a deformable droplet", Int. J. Heat Mass Transfer 73 (2014) 616-626.
#
# WHAT THIS REPRODUCES
# Samareh's FINITE-Marangoni validation case in 2D (their Section 4.1.2 / Fig 7, "Thermocapillary
# migration of a droplet with low Marangoni number"), the test originally devised by Nas & Tryggvason
# [Samareh ref. 19; Int. J. Multiphase Flow 29 (2003) 1117-1135]. Unlike the Fig 5 case (case.py,
# Ma = 0, a FROZEN imposed temperature field), here the Marangoni number is finite, so the ENERGY
# equation is coupled: the droplet's motion distorts the temperature field, a thermal boundary layer
# forms at the interface, and surface tension responds to the evolving local T. This is the case that
# actually exercises bulk Fourier CONDUCTION -- it is the paper's validation of the heat-transfer
# module, and the right home for MFC's thermal_conduction + thermal_scalar features.
#
# Samareh's Fig 7 setup (their Sec. 4.1.2, verbatim):
#   "a two-dimensional droplet is placed inside a rectangular box of size 2 x 4 droplet diameters.
#    The center of the droplet is initially located 1 droplet diameter away from the bottom wall.
#    The ratio of the droplet material properties to that of the bulk fluid is set equal to 0.5 and
#    the non-dimensional numbers are selected as Re = 5, Ma = 20, and Ca = 0.01666."
#   Result: U* = U/U_r peaks at ~0.13 (near t* = t/t_r ~ 5), then relaxes toward a terminal value;
#   the fine-grid (n_x = 128) terminal velocity matches Nas & Tryggvason to within 1.7%, and the
#   64 vs 128 grids differ by 1.2%. (Note: Samareh report U* = U/U_r here, NOT v/v_YGB -- the
#   v/v_YGB normalization is used only for the Ma = 0 Fig 5 case in case.py.)
#
# WHY thermal_scalar (decoupled T), NOT the density-proxy
# The Fig 5 case (case.py) fakes the imposed T as a density gradient (rho = rho_coeff/T), which only
# works because there both fluids are IDENTICAL. Fig 7's droplet is a genuinely DIFFERENT fluid (all
# material properties at ratio 0.5), so density is set by the two-fluid composition and CANNOT also
# encode T. Temperature must therefore be carried as MFC's INDEPENDENT advected+diffused scalar
# (thermal_scalar -> eqn_idx%T_s), decoupled from the EOS, conducting at alpha = k/(rho cp). sigma(T)
# reads T_s. This is exactly the configuration thermal_scalar was built for.
#
# NON-DIMENSIONAL INVERSION (Samareh Eqs. 13-15; bulk phase = subscript b)
#   length scale r_0 (drop radius), velocity scale U_r = |sigma_T gradT| r_0 / mu_b,
#   time scale t_r = mu_b / |sigma_T gradT|, temperature scale T_r = |gradT| r_0.
#   Ca = |sigma_T gradT| r_0 / sigma_0 ,  Re = |sigma_T gradT| rho_b r_0^2 / mu_b^2 ,
#   Ma = U_r r_0 / alpha_b = |sigma_T gradT| r_0^2 / (mu_b alpha_b) ,  alpha_b = k_b/(rho_b cp_b).
# Writing G = |sigma_T gradT| (the Marangoni stress magnitude), and choosing r_0, rho_b, mu_b freely:
#   G = Re mu_b^2 / (rho_b r_0^2);  sigma_0 = G r_0 / Ca;  alpha_b = G r_0^2 / (mu_b Ma).
# The split G = |sigma_T| * gradT is free; we fix gradT so T runs 0..1 across the box (Samareh's
# Fig 5 convention), then sigma_T = -G/gradT (sigma falls with T, so the drop rises toward hot).
# rho_b = 1, mu_b = 0.02 are picked so the Marangoni velocity scale U_r and the sound speed give a
# deeply incompressible state (Mach = U_r/c ~ 0.02); the absolute scales are arbitrary, only the
# three non-dimensional numbers are physical.
#
# ONE ASSUMPTION TO VERIFY: "material properties ratio = 0.5" is taken to mean the droplet's density,
# viscosity, conductivity AND specific heat are each 0.5x the bulk's (the conventional reading). With
# all four at 0.5 the droplet thermal diffusivity is alpha_d = k_d/(rho_d cp_d) = 2 alpha_b. If Nas &
# Tryggvason 2003 instead held the volumetric heat capacity (rho cp) ratio fixed (so alpha_d = alpha_b),
# flip PROP_RATIO handling for cv below -- it is a single knob (prop_ratio) applied to rho/mu/k/cv.
#
# GOTCHA (analytic-IC parser): MFC expands a bare `e` to Euler's number -- even inside a literal like
# 1e-9. Keep `e`-notation OUT of every analytic patch string (see eps / T_expr below).

import json
import os

# -- Variant selection (env vars; defaults = Samareh Fig 7 at the coarse grid) --
width_cells = int(os.environ.get("FIG7_NX", "64"))  # cells across the 2D box WIDTH (Samareh: 64, 128)
n_tr = float(os.environ.get("FIG7_TR", "15"))  # run length in capillary-thermal times t_r (Fig 7 x-axis ~ 0..20)
prop_ratio = float(os.environ.get("FIG7_RATIO", "0.5"))  # droplet/bulk material-property ratio (Samareh: 0.5)
conduction = os.environ.get("FIG7_COND", "1") == "1"  # diagnostic: FIG7_COND=0 -> advection of T_s only (no Fourier)
adiabatic = os.environ.get("FIG7_ADIABATIC", "0") == "1"  # diagnostic: adiabatic walls instead of the faithful isothermal Dirichlet
# Acoustic-ring control. The box is compressible and closed (reflective -2 walls); a curved interface
# at uniform pressure leaves the Laplace jump sigma/r unbalanced at t=0, which launches an acoustic
# standing wave that the walls trap (a domain-filling vertical mode, fundamental f = c/2Ly), riding on
# the slow migration as a ~few-% ripple. We remove it at the source by initializing the droplet at its
# Laplace overpressure (see patch_icpp(2)%pres below): an ablation showed this cuts the ripple by ~97%,
# down to the diffuse-interface parasitic-current floor -- bulk viscosity / EOS stiffening then add
# nothing. FIG7_UNBALANCED=1 restores Samareh's bare uniform-pressure IC (rings) for that diagnostic.
unbalanced_ic = os.environ.get("FIG7_UNBALANCED", "0") == "1"

# -- Samareh's Fig 7 non-dimensional targets --
Re = 5.0
Ma = 20.0
Ca = 0.01666

# -- Geometry (2D droplet in a 2D-wide x 4D-tall box; rise/gradient axis = y) --
D = 1.0  # droplet diameter
r0 = D / 2.0  # droplet radius = 0.5 (the length scale)
Wx = 2.0 * D  # box width (short axis, x)
Hy = 4.0 * D  # box height (the gradient / rise axis, y)
y_bottom = -Hy / 2.0  # cold floor at y = -2
y_drop = y_bottom + 1.0 * D  # drop center 1D above the cold floor -> y = -1 (Samareh)

dx = Wx / width_cells  # isotropic cell size set by the box-width resolution
height_cells = round(Hy / dx)  # = 2 * width_cells
m = width_cells - 1  # x (short axis)        MFC index = cells - 1
n = height_cells - 1  # y (long / rise axis)
p = 0  # 2D

# Diffuse color interface ~2 cells wide at every resolution (smooth_coeff = dx/w, w = 2*dx -> 0.5).
cf_smooth_coeff = 0.5

# -- Bulk-phase reference scales (arbitrary; only Re, Ma, Ca are physical) --
rho_b = 1.0  # bulk density
mu_b = 0.02  # bulk dynamic viscosity (small -> deeply incompressible Marangoni flow)
gam = 2.0  # stiffened-gas exponent (same for both fluids: equal compressibility)
cv_b = 1.0  # bulk specific heat at constant volume (free; only enters conduction via cp = gam*cv)

# -- Invert the non-dimensional numbers for the physical (sigma, gradT, alpha, k) --
G = Re * mu_b**2 / (rho_b * r0**2)  # |sigma_T * gradT| (Marangoni stress magnitude) = 0.008
gradT = 1.0 / Hy  # |dT/dy|: T runs 0 (cold floor) .. 1 (hot ceiling) over the 4D box = 0.25
sigma_T = -G / gradT  # dsigma/dT < 0 (sigma falls with T) = -0.032
sigma0 = G * r0 / Ca  # surface tension at T_ref (Ca closure) ~ 0.240
alpha_b = G * r0**2 / (mu_b * Ma)  # bulk thermal diffusivity = 0.005
cp_b = gam * cv_b  # bulk specific heat at constant pressure = 2.0
k_b = alpha_b * rho_b * cp_b  # bulk conductivity (so alpha_b = k_b/(rho_b cp_b)) = 0.01

# -- Droplet properties = prop_ratio x bulk (the Samareh 0.5 material-property ratio) --
rho_d = prop_ratio * rho_b  # = 0.5
mu_d = prop_ratio * mu_b  # = 0.01
cv_d = prop_ratio * cv_b  # = 0.5  (-> cp_d = gam*cv_d = 1.0)
k_d = prop_ratio * k_b  # = 0.005 (-> alpha_d = k_d/(rho_d cp_d) = 2 alpha_b; see header note)

# -- Diagnostic scales (Samareh's normalization of Fig 7) --
U_r = G * r0 / mu_b  # Marangoni velocity scale |sigma_T gradT| r_0/mu_b = 0.2
t_r = mu_b / G  # capillary-thermal time scale mu_b/|sigma_T gradT| = 2.5

# -- Equation of state (stiffened gas; background pressure well above the Laplace jump sigma0/r0) --
# The flow is low-Mach; the sound speed (set by p0, p_inf) fixes only the acoustic CFL, not U_r.
p0, p_inf = 5.0, 20.0  # Laplace jump sigma0/r0 ~ 0.48 << p0; c stays moderate
rho_min = rho_d  # lowest density is the droplet -> sets the max sound speed / acoustic CFL
c_max = (gam * (p0 + p_inf) / rho_min) ** 0.5  # ~ 10

# -- Imposed temperature field T(y) = T_base + gradT*(y - y_bottom): T runs T_base .. T_base+1 --
# A uniform baseline T_base is added because the isothermal-BC validator requires Twall > 0 (an EOS
# thermodynamics check). For the DECOUPLED scalar T_s the offset is physically inert: sigma depends
# only on (T - T_ref), and both advection and conduction of T_s are shift-invariant. So T runs
# 1 (cold floor) .. 2 (hot ceiling); only the gradient gradT and slope sigma_T are physical.
T_base = 1.0  # positive baseline so Twall > 0 (inert for the decoupled T_s; cf. case.py's T0 shift)


def T_of_y(y):
    return T_base + gradT * (y - y_bottom)


T_ref = T_of_y(y_drop)  # T at the drop's initial position -> sigma = sigma0 there (= 1.25)
GRAD = "y"  # gradient / rise axis
# Plain decimals only (no `e` notation): MFC would expand a bare `e` to Euler's number in the IC parser.
T_expr = f"{T_base} + {gradT}*({GRAD} - ({y_bottom}))"  # = 1.0 + 0.25*(y - (-2.0)) ; T runs 1..2
eps = 1.0e-9  # trace volume fraction of the "other" fluid in each patch (kept out of analytic strings)

# -- Time stepping: min(acoustic CFL, explicit-conduction limit) --
mydt = 0.35 * dx / c_max  # acoustic CFL (RK3 + WENO5, ICFL 0.35)
alpha_max = max(alpha_b, k_d / (rho_d * gam * cv_d))  # fastest-diffusing phase sets the conduction dt
mydt = min(mydt, 0.35 * dx**2 / (4.0 * alpha_max))  # explicit diffusion number, d = 2 (2D)
t_end = n_tr * t_r  # default 15 t_r = 37.5 (covers the Fig 7 ramp, overshoot, and terminal approach)
t_step_stop = int(round(t_end / mydt))
t_step_save = max(1, t_step_stop // 100)  # ~100 snapshots to resolve the U*(t*) curve

data = {
    # Logistics
    "run_time_info": "T",
    # Computational domain: rise (gradient) axis y in [-2, 2]; short axis x in [-1, 1]
    "x_domain%beg": -Wx / 2,
    "x_domain%end": Wx / 2,
    "y_domain%beg": -Hy / 2,
    "y_domain%end": Hy / 2,
    "m": m,
    "n": n,
    "p": p,
    "cyl_coord": "F",
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
    # Boundaries: closed slip-wall box (-2 on all sides), Samareh's / Nas & Tryggvason's geometry.
    "bc_x%beg": -2,
    "bc_x%end": -2,
    "bc_y%beg": -2,
    "bc_y%end": -2,
    "num_patches": 2,
    "num_fluids": 2,
    # Physics: viscosity + bulk conduction + temperature-dependent surface tension sigma(T)
    "viscous": "T",
    "surface_tension": "T",
    "sigma": sigma0,
    "sigma_model": 1,
    "sigma_T_ref": T_ref,
    "sigma_dTdT": sigma_T,
    "thermal_scalar": "T",  # carry T as an independent scalar T_s, decoupled from the EOS/density
    # Output
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T",
    "cf_wrt": "T",
    "T_s_wrt": "T",  # expose the independent temperature scalar for visualization
    "parallel_io": "T",
    # Continuous phase (fluid 1 = bulk)
    "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(1)%cv": cv_b,
    "fluid_pp(1)%Re(1)": 1.0 / mu_b,
    # Dispersed phase (fluid 2 = droplet), all material properties at prop_ratio x bulk
    "fluid_pp(2)%gamma": 1.0 / (gam - 1.0),
    "fluid_pp(2)%pi_inf": gam * p_inf / (gam - 1.0),
    "fluid_pp(2)%cv": cv_d,
    "fluid_pp(2)%Re(1)": 1.0 / mu_d,
    # Patch 1 -- bulk medium spanning the domain: density rho_b, ambient color c = 0, T_s(y) linear.
    "patch_icpp(1)%geometry": 3,  # 2D rectangle
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%length_x": Wx,
    "patch_icpp(1)%length_y": Hy,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": (1.0 - eps) * rho_b,
    "patch_icpp(1)%alpha_rho(2)": eps * rho_d,
    "patch_icpp(1)%alpha(1)": 1.0 - eps,
    "patch_icpp(1)%alpha(2)": eps,
    "patch_icpp(1)%cf_val": 0.0,
    "patch_icpp(1)%T_temp_val": T_expr,
    # Patch 2 -- droplet (2D circle): fluid 2 at rho_d, color c = 1, smeared over ~2 cells. SAME
    # imposed T_s(y) as the bulk (temperature is continuous across the interface; only the color and
    # the material properties jump). The smoothed alter_patch sets the real two-fluid composition.
    "patch_icpp(2)%geometry": 2,  # 2D circle
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": y_drop,
    "patch_icpp(2)%radius": r0,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%smoothen": "T",
    "patch_icpp(2)%smooth_patch_id": 1,
    "patch_icpp(2)%smooth_coeff": cf_smooth_coeff,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    # Initialize the droplet interior at the Laplace overpressure p_out + sigma/r (2D circle: single
    # curvature 1/r) so the t=0 interface is in mechanical equilibrium -- the smoothed patch tanh-blends
    # this into a balanced pressure profile, removing the unbalanced CSF kick that otherwise launches
    # the acoustic standing wave. FIG7_UNBALANCED=1 reverts to Samareh's bare uniform-pressure IC.
    "patch_icpp(2)%pres": p0 if unbalanced_ic else p0 + sigma0 / r0,
    "patch_icpp(2)%alpha_rho(1)": eps * rho_b,
    "patch_icpp(2)%alpha_rho(2)": (1.0 - eps) * rho_d,
    "patch_icpp(2)%alpha(1)": eps,
    "patch_icpp(2)%alpha(2)": 1.0 - eps,
    "patch_icpp(2)%cf_val": 1.0,
    "patch_icpp(2)%T_temp_val": T_expr,
}

if conduction:
    # Finite Ma: solve the energy transport (Fourier -k grad T) for the decoupled scalar T_s, at
    # alpha = k/(rho cp). Bulk conduction time domain^2/alpha ~ 3200 >> t_end ~ 37.5, so the far-field
    # gradient is effectively held by the IC; the drop-local thermal boundary layer is the physics.
    data.update(
        {
            "thermal_conduction": "T",
            "fluid_pp(1)%k_therm": k_b,  # bulk conductivity (alpha_b = k_b/(rho_b cp_b))
            "fluid_pp(2)%k_therm": k_d,  # droplet conductivity (prop_ratio x bulk)
        }
    )
    if not adiabatic:
        # Faithful Samareh/Nas & Tryggvason setup (DEFAULT): isothermal Dirichlet on the top/bottom
        # walls pins T to the imposed gradient (T = T_base cold floor, T_base + 1 hot ceiling), so the
        # walls SINK the drop's thermal wake and the migration reaches a clean plateau. The side walls
        # (bc_x) are left non-isothermal -> adiabatic. NOTE: this requires the MPI rank-guard fix in
        # s_apply_thermal_conduction_bc -- before it, the isothermal BC overwrote interior ranks' halo
        # cells, scrambling T into rank-boundary bands and REVERSING the drop (this, not a density
        # proxy or advective throughflow, was the CONDUCTION_REVERSAL_SAGA root cause). FIG7_ADIABATIC=1
        # leaves the walls adiabatic instead (gradient IC-sustained over the window; the wake is not
        # sunk, so the plateau over-declines) -- a diagnostic, not the faithful run.
        data.update(
            {
                "bc_y%isothermal_in": "T",
                "bc_y%isothermal_out": "T",
                "bc_y%Twall_in": T_of_y(-Hy / 2),  # cold floor (y%beg)
                "bc_y%Twall_out": T_of_y(Hy / 2),  # hot ceiling (y%end)
            }
        )

print(json.dumps(data))
