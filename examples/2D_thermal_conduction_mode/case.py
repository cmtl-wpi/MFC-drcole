#!/usr/bin/env python3
# 2D heat equation, decaying Fourier mode on a doubly-periodic square.
# T = T0 + A*sin(k*x)*sin(k*y), k = 2*pi/L; decays in place as exp(-alpha*2*k^2*t).
# A periodic, smooth mode (no wall step) stays quiescent, so it validates the
# operator cleanly in the compressible energy path -- unlike a strong-Dirichlet
# plate, whose hot wall drives a sustained near-sonic flow. Temperature is
# recovered from the stiffened-gas EOS: the mode is imposed through the density
# at uniform pressure, rho = rho0*T0/T.
import json
import math

L = 1.0
gam, p_inf, cv = 2.0, 100.0, 12.5
p0, T0, A = 25.0, 10.0, 3.0
alpha = 0.05
rho0 = (p0 + p_inf) / ((gam - 1.0) * cv * T0)
k_therm = alpha * rho0 * cv * gam
kw = 2.0 * math.pi / L

Nx = Ny = 127
dx = L / (Nx + 1)
c0 = (gam * (p0 + p_inf)) ** 0.5
# stability needs the cv-based diffusivity k/(rho*cv) = gam*alpha, not the cp-based alpha that sets the decay rate
dt = min(0.4 * dx / c0, 0.35 * dx**2 / (4.0 * gam * alpha))
t_end = 1.0 / (alpha * 2.0 * kw**2)  # one e-fold
Nt = int(round(t_end / dt))

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": 0.0,
            "x_domain%end": L,
            "y_domain%beg": 0.0,
            "y_domain%end": L,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": Nt,
            "t_step_save": max(1, Nt // 30),
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "model_eqns": 2,
            "num_fluids": 1,
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "parallel_io": "T",
            # Thermal conduction (EOS temperature)
            "thermal_conduction": "T",
            # Fluids Physical Parameters
            "fluid_pp(1)%gamma": 1.0 / (gam - 1.0),
            "fluid_pp(1)%pi_inf": gam * p_inf / (gam - 1.0),
            "fluid_pp(1)%cv": cv,
            "fluid_pp(1)%k_therm": k_therm,
            # Patch 1: full square, periodic sine mode imposed through density
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": L / 2,
            "patch_icpp(1)%y_centroid": L / 2,
            "patch_icpp(1)%length_x": L,
            "patch_icpp(1)%length_y": L,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": f"{rho0 * T0:.12f}/({T0} + {A}*sin({kw:.12f}*x)*sin({kw:.12f}*y))",
            "patch_icpp(1)%alpha(1)": 1.0,
        }
    )
)
