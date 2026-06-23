---
name: thermocap-nan-pressure-relaxation
description: case_Ma_20 NaN was a 0/0 in the 6-eq pressure relaxation when the bulk phase vanishes inside the drop; fixed by clamping rho_K_s. NOT a GPU bug, NOT a wall corner.
metadata:
  type: decision
---

**Root cause (found & fixed).** `case_Ma_20.py` (TC2) aborted with "NaN(s) in timestep
output" partway through (originally seen at the save step 15147; true ignition **t_step
14870**). The blowup is a **0/0 in the six-equation phase-pressure relaxation**
(`model_eqns=3`), `src/simulation/m_pressure_relaxation.fpp`, `s_equilibrate_pressure`.

Mechanism: the rising drop is nearly pure fluid 2, so the bulk phase drains out of the drop
interior — `alpha_rho_1` (fluid-1 partial density) decays to **exactly 0.0** while `alpha_1`
is still just above `sgm_eps` (1e-16). The Newton loop forms `rho_K_s(i) =
alpha_rho_i/max(alpha_i,sgm_eps)*(p-ratio)^(1/gs_min)`, which is `0` when `alpha_rho_i=0`;
then `alpha_rho_i/rho_K_s(i)` (lines 195, 196, 207) is `0/0 = NaN`. The only guard was
`if (alpha_i > sgm_eps)` — it protects a vanishing *volume fraction*, not a vanishing
*partial density*. `mpp_lim` doesn't help (it only zeros *negative* values). Pressure/sound
speed stay healthy (p≈5.85), so it is NOT negative-pressure; it's a pure FP-invalid event,
which is why it's a single-step NaN from a smooth state and **platform-independent**
(CPU≡GPU). Confirmed by an `nvfortran -Ktrap=fp` debug build trapping at the exact line at
t_step 14870, and by the source audit.

**Corrects the earlier note** (was `gpu-nan-thermal-marangoni.md`): it is **not** a GPU bug,
**not** the top-right wall corner (the abort print's trailing `63 127 0` is `m,n,p` grid
dims, not a cell index — the real first-NaN cell is reported as `j,k,l,i`), and the cause was
NOT in the new σ(T)/conduction kernels. The new physics only *triggers* it: σ(T) Marangoni +
conduction drive the migration that drains the phase. With `sigma_model=0` there's no
migration, the phase never fully drains, so it never hit the latent divide.

**Fix:** clamp `rho_K_s(i) = max(rho_K_s(i), sgm_eps)` right after its assignment (covers all
three divides). `rho_K_s` is an intrinsic density (O(1)) except in the exact-zero degenerate
cell, where the massless phase then contributes 0 — physically correct, negligible elsewhere.
The divide is pre-existing (file ≈ unchanged from master), so this hardens *all* `model_eqns=3`
cases. **Verified:** GPU rebuild completes the full run to t_step 45999 (was dying at 15147)
with zero NaN; a per-step NaN-checked CPU debug rerun is clean from 14868 through 15366
(covering the 14870 ignition and the old 15147 abort). (The pre-fix root cause was pinned by an
`nvfortran -Ktrap=fp` build that trapped exactly at the 0/0 line at t_step 14870.)

Complementary case-quality note: `case_Ma_20.py` uses `weno_eps=1e-16` (vs the usual 1e-6),
which sharpens interface overshoot that pushes `alpha_rho_1`→0. The code clamp is the real
fix; a saner `weno_eps`/volume-fraction floor would keep the case better-posed. Build/run env
on this workstation: [[gpu-build-on-nighthawk]].
