#!/usr/bin/env python3
"""2D thermal flat-plate boundary layer -- non-chemistry bulk Fourier conduction.

This case is the thermal_conduction counterpart of 2D_Thermal_Flatplate.
Instead of chemistry + chem_params%diffusion for heat conduction, it uses the
standalone thermal_conduction feature (Fourier's law -k*grad(T)) with the
mixture stiffened-gas EOS temperature. All other physics -- domain, grid,
boundary conditions, viscous, and initial state -- are kept identical so the
two cases can be directly compared.

Physical setup:
  - 2D square domain, 0.05 m x 0.05 m
  - Single fluid: ideal-gas air (gamma = 1.4, pi_inf = 0)
  - Initial state: quiescent, uniform at 1125 K, 1 atm
  - Left boundary: subsonic characteristic inflow
  - Right boundary: ghost-cell extrapolation outflow
  - Bottom boundary: no-slip isothermal wall at 600 K (the flat plate)
  - Top boundary: ghost-cell extrapolation outflow
  - Viscous: Re = 100,000
  - Bulk thermal conduction: k = 0.07 W/(m*K)  (air at ~1000 K)

The 525 K temperature difference between the free stream (1125 K) and the wall
(600 K) drives a thermal boundary layer through Fourier conduction. The
isothermal Dirichlet wall BC pins the ghost-cell temperature so that the face
flux through the wall is consistent with the prescribed 600 K wall temperature.
Boundaries without an isothermal flag are adiabatic (zero-gradient) by default.
"""

import json

Lx = 0.05  # m
Ly = 0.05  # m

# Air properties (ideal gas, gamma = 1.4)
c_gamma = 1.4  # specific heat ratio
R_air = 287.0  # J/(kg*K)
cv = R_air / (c_gamma - 1.0)  # ≈ 717.5 J/(kg*K)
p0 = 101325.0  # Pa (1 atm)
T0 = 1125.0  # K (free-stream temperature)

# Thermal conductivity of air near 1000 K
k_air = 0.07  # W/(m*K)

# Grid: 499 x 499 cells (isotropic dx ~ 1.00e-4 m)
m_cells = 499
n_cells = 499

# Acoustic CFL limit: c = sqrt(gamma*R*T) ~ 672 m/s
c_sound = (c_gamma * R_air * T0) ** 0.5  # ~ 672 m/s
dx = Lx / m_cells  # ~ 1.00e-4 m
dt_acoustic = 0.3 * dx / c_sound  # ~ 4.47e-8 s

# Diffusion CFL limit: alpha = k/(rho*cp), rho = p/(R*T)
rho0 = p0 / (R_air * T0)  # ~ 0.314 kg/m^3
cp = c_gamma * R_air / (c_gamma - 1.0)  # ~ 1004.5 J/(kg*K)
alpha_T = k_air / (rho0 * cp)  # ~ 2.22e-4 m^2/s
dt_diff = 0.3 * dx**2 / (4.0 * alpha_T)  # ~ 3.39e-6 s  (not limiting)

dt = min(dt_acoustic, dt_diff)  # acoustic-CFL limited

t_end = 5.0e-3  # 5 ms simulated time
t_step_stop = int(round(t_end / dt))
t_step_save = max(1, t_step_stop // 10)  # ~10 snapshots

case = {
    # -- Run info --
    "run_time_info": "T",
    # -- Domain --
    "x_domain%beg": 0.0,
    "x_domain%end": Lx,
    "y_domain%beg": 0.0,
    "y_domain%end": Ly,
    "m": m_cells,
    "n": n_cells,
    "p": 0,
    # -- Time stepping --
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": t_step_stop,
    "t_step_save": t_step_save,
    # -- Model --
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "num_fluids": 1,
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": 3,
    # -- Spatial scheme --
    "mp_weno": "F",
    "weno_order": 5,
    "weno_eps": 1e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 2,
    # -- Boundary conditions --
    "bc_x%beg": -7,  # subsonic characteristic inflow
    "bc_x%end": -3,  # ghost-cell extrapolation outflow
    "bc_y%beg": -16,  # no-slip isothermal wall (flat plate)
    "bc_y%end": -3,  # ghost-cell extrapolation outflow
    "bc_y%isothermal_in": "T",
    "bc_y%Twall_in": 600.0,
    # -- Output --
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    # -- Physics: viscous --
    "viscous": "T",
    "fluid_pp(1)%Re(1)": 100000,
    # -- Physics: bulk thermal conduction (non-chemistry) --
    "thermal_conduction": "T",
    "chemistry": "F",
    # -- Fluid properties --
    # MFC stores gamma as 1/(gamma-1) in the stiffened-gas EOS
    "fluid_pp(1)%gamma": 1.0 / (c_gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(1)%cv": cv,
    "fluid_pp(1)%k_therm": k_air,
    # -- Initial condition: single uniform patch --
    # Quiescent, uniform free stream at T0 = 1125 K, 1 atm. The partial
    # density must be rho0 = p0/(R*T0) so the EOS-derived temperature is
    # actually 1125 K (T = (p + pi_inf)/((gamma-1)*rho*cv)).
    "num_patches": 1,
    "patch_icpp(1)%geometry": 3,  # 2D rectangle
    "patch_icpp(1)%x_centroid": Lx / 2,
    "patch_icpp(1)%y_centroid": Ly / 2,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%vel(1)": 0,
    "patch_icpp(1)%vel(2)": 0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha(1)": 1,
}

if __name__ == "__main__":
    print(json.dumps(case))
