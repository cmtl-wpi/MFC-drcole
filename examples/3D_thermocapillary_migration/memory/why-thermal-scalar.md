---
name: why-thermal-scalar
description: "The u_YGB validation (case_ygb.py) uses thermal_scalar=T, NOT the density-proxy temperature of case.py — the proxy advects with the drop's flow and reverses the gradient"
metadata:
  type: decision
---

`case_ygb.py` carries temperature as an independent advected+diffused scalar `T_s`
(`thermal_scalar=T`), with uniform density. The legacy `case.py` instead fakes temperature through
density (`rho = rho_coeff/T(y)`), with conduction acting on the EOS temperature.

**Why the scalar:** the density proxy is a *transported* field, so the drop's own recirculating flow
advects the temperature it is meant to hold fixed — the local interfacial gradient collapses and
reverses, and the rise velocity decays (measured in the 2D frozen-T case; see the project memory
`frozen-t-proxy-advects-not-frozen`). With `thermal_scalar=T`, the surface-tension closure reads
`T_s` directly (`src/simulation/m_surface_tension.fpp`, the `if (thermal_scalar)` branch), decoupled
from density. Both fluids are identical (μ*=k*=1) and density is uniform, so the ONLY thing driving
the drop is the σ(T) gradient — the clean YGB setup, and a *proper* validation of variable surface
tension rather than an EOS artifact.

**How to apply:** for any thermocapillary validation here, use `case_ygb.py` (scalar), not `case.py`
(proxy). `measure.py`/`fields_ygb.py` already handle `ts_mode` (color at `nvars-2`, `T_s` at
`nvars-1`). Verified 2026-06-19 on a cube/W5/Nx40 smoke run: `T_s` initializes exactly linear
(slope 0.13332 vs imposed 0.13333) and wall-pinned, and the drop rises toward the hot wall (u/v_YGB
ramps 0.34→0.48 over the first 0.03 t_r). See [[confinement-to-one]].
