# Temperature-dependent surface tension — σ(T), the thermal Marangoni effect

Branch `pr/variable-sigma`, 6 commits on top of the thermal-conduction base
(`origin/pr/thermal-conduction..origin/pr/variable-sigma`, `8e06deac`..`63260266`). Line numbers
cite that branch tip.

Adds a temperature-dependent surface tension. The elegant part: it introduces **no explicit
Marangoni source term**. Making σ a spatially varying field is *sufficient* — the tangential
Marangoni stress falls out of the divergence of the existing capillary stress tensor.

## 1. The physics

Surface tension σ is the energy per unit interfacial area. When σ is **uniform**, its only effect
is the normal Laplace pressure σκ. When σ **varies along the interface**, the tangential gradient
∇ₛσ is an unbalanced tangential force that drags interfacial fluid from low- to high-tension
regions — the **Marangoni stress**. For a freely suspended drop the reaction propels the drop
bodily (**thermocapillary migration**). Since σ usually falls with T (`dσ/dT < 0`), a drop in a
temperature gradient has lower tension on its hot side and migrates **toward the hot** region.

The feature implements a single linear model, `sigma_model = 1`:

```
σ(T) = sigma + sigma_dTdT · (T − sigma_T_ref)
```

with exact parameter names:
- `sigma` — reference surface tension (the existing constant-σ coefficient).
- `sigma_dTdT` — the slope `dσ/dT` (signed, typically negative).
- `sigma_T_ref` — the reference temperature at which σ = `sigma`.

Fortran evaluation, `src/simulation/m_surface_tension.fpp:320`:

```fortran
c_sigma(j, k, l) = sigma + sigma_dTdT*(T_cell - sigma_T_ref)
```

Documented validity: the closure is physical only while σ(T) > 0, i.e. T within
`T_ref ± sigma/|dσ/dT|`.

## 2. The injection point

MFC models surface tension with the **Continuum-Surface-Force (CSF)** formulation (Schmidmayer et
al. 2017): a capillary stress tensor built from the color-function gradient **w** = ∇c,

```
Ω = σ·( |∇c|·I − ∇c⊗∇c/|∇c| )
```

whose divergence ∇·Ω is added to the momentum flux (and a companion σ(∇·n) term to the energy
flux), per direction in `s_capilliary`. **With constant σ, ∇·Ω is purely the normal capillary
force.** The moment σ becomes a spatially varying field σ(x), the divergence picks up an extra
`(∇σ)·(…)` piece whose tangential component *is* the Marangoni stress. Nothing else is added.

Mechanically the feature:
1. Adds a cell-centered field `c_sigma(j,k,l)` (`m_surface_tension.fpp:43`), filled from σ(T) over
   the full buffer range when `sigma_model == 1` (`m_surface_tension.fpp:313`), reusing
   `f_compute_mixture_temperature` from the conduction base (this range does **not** touch
   `m_sim_helpers.fpp`).
2. At each face forms a face-local σ by averaging the two adjacent cells, e.g. x-direction
   (`m_surface_tension.fpp:116`):
   ```fortran
   sigma_face = sigma
   if (sigma_model == 1) sigma_face = (c_sigma(j,k,l) + c_sigma(j+1,k,l))/2._wp
   ```
3. Passes `sigma_face` into the stress-tensor macro and the energy flux. When `sigma_model == 0`,
   `sigma_face` stays equal to the constant `sigma` and the path is **bit-for-bit unchanged**.

`c_sigma` is allocated **unconditionally** whenever surface tension is active (not gated on
`sigma_model==1`) so the device descriptor is always valid in the capillary kernel; it is only
written/read when `sigma_model==1`. Full-buffer computation reuses already-populated `q_prim_vf`
ghost cells, so no extra halo exchange.

## 3. How it threads the architecture

- **Params** (`toolchain/mfc/params/definitions.py:627`): `sigma_model` (INT), `sigma_T_ref`,
  `sigma_dTdT` (REAL), all plain scalars in the namelist registry → Fortran declarations
  **auto-generate** (none appear in `m_derived_types.fpp`). Defaults
  (`m_global_parameters_common.fpp:447`): `sigma_model=0`, `sigma_T_ref=dflt_real` (sentinel),
  `sigma_dTdT=0._wp` (a genuine no-op slope). That default asymmetry is what makes the
  `966dc60d` validator fix necessary.
- **The computation** — all in `m_surface_tension.fpp` (the CSF module), above. The base commit
  also parametrized the stress-tensor macro to take a `sig` argument (see the parenthesization
  gotcha), and fixed a pre-existing `c_divs` finalize deallocate off-by-one
  (`do j=1,num_dims` → `do j=1,num_dims+1`).
- **GPU** (`d79874fc`): the three scalars are read inside device kernels, so they were added to
  the same `GPU_DECLARE` and `GPU_UPDATE(device=...)` lines as `sigma`
  (`m_global_parameters_common.fpp:96`, `m_start_up.fpp:1042`). Without the device push, the
  device copy holds the default `sigma_model=0`, silently disabling the closure on GPU. Same class
  of bug as the thermal-conduction GPU-declare fix. `sigma_face` was added to the `private=[...]`
  clause of all three directional capillary loops.
- **Validator** (`case_validator.py:755`): `sigma_model ∈ {0,1}`; model 1 / `sigma_dTdT` /
  `sigma_T_ref` require `surface_tension`; model 1 requires `sigma_dTdT` and (after `966dc60d`)
  `sigma_T_ref`; `cv > 0` for both fluids (temperature recovery divides by `mCP`).

## 4. Gotchas and subtleties

- **The macro-argument parenthesization fix** (`d79874fc`). The base commit parametrized the Fypp
  macro `compute_capillary_stress_tensor` with a **bare** expansion `-${sig}$*(...)`. Fypp does
  *textual* substitution. With today's callers (`sig = sigma_face`, a bare identifier) it is
  correct, so the bug is **latent**. But if a future call site passed an *expression* — say
  `a + b` — the bare expansion would become:
  ```fortran
  Omega(1,1) = -a + b*(w2*w2 + w3*w3)/normW      ! parses as (-a) + (b*…) — WRONG
  ```
  because unary minus and `*` bind tighter than the injected `+`. The intended `-(a+b)*(…)`
  silently becomes `-a + b*(…)`. Wrapping every expansion as `-(${sig}$)*(…)` makes the macro
  correct for its whole contract, not just today's caller.
- **The orphan-parameter validator fix** (`966dc60d`). `sigma_T_ref` defaults to the `dflt_real`
  sentinel, so an unset ref would silently poison σ(T) = σ + slope·(T − garbage) — hence the new
  hard requirement. Separately, the early-return guard was widened so a case that sets
  `sigma_dTdT`/`sigma_T_ref` *without* `sigma_model` (a likely typo) no longer returns early and
  skips the "require surface_tension" check.
- **Precision**: `c_sigma`, `sigma_face`, `T_cell` are `wp`. `c_sigma` is a transient computed
  field, not a stored/I-O field, so it correctly uses `wp` not `stp`. No mixing.

## 5. Validation

**Example — `2D_thermocapillary_migration`.** A neutrally-buoyant 2D drop in an imposed linear
temperature field migrating toward the hot wall, validated against Samareh et al. (IJHMT 73, 2014)
for the two σ(T)-only cases. The analytic reference is the **Young-Goldstein-Block** terminal
speed (`measure.py:153`):

```python
v_YGB = (2.0/15.0) * (-dsigma_dT) * gradT * R / mu     # μ*=k*=1
```

Measured plateau `v_t/v_YGB ≈ 0.80` (finite slip-wall box; the unbounded-cylinder limit is
15/16 = 0.938). Setup notes: σ params derived from target Marangoni/Capillary numbers; the linear
`T(y)` is imposed through the **density IC** at uniform pressure, painted as a **single full-box
analytic patch** so T stays continuous across the interface (a two-patch drop IC injects a
spurious force exactly where the Marangoni force lives, ~22% error — quantified by a companion
`_2patch.py` case). Uses `thermal_conduction=T` + isothermal walls to sustain the gradient — so it
builds on the conduction base.

**Golden test — `E39EFF77`** (`toolchain/mfc/test/cases.py`, 2D only): `sigma_model=1`,
`sigma_T_ref=2.0`, `sigma_dTdT=-0.1`, with **differing per-fluid cv (1.0 vs 2.0)** to exercise the
mixture-`mCP` temperature recovery, layered on the existing two-interface capillary case (50 steps).

## 6. Commit-by-commit

1. `8e06deac` Add temperature-dependent surface tension σ(T) — params, defaults, parametrized
   stress-tensor macro, `c_sigma` field + linear closure, face-averaged `sigma_face` in the
   momentum/energy fluxes, first-cut validator; fixes a pre-existing `c_divs` deallocate off-by-one.
2. `b1756d75` Add the `2D_thermocapillary_migration` example (vs Samareh 2014 / YGB).
3. `966dc60d` case_validator: require `sigma_T_ref`; fix the orphan-parameter early-return.
4. `d79874fc` device-declare the σ(T) scalars; parenthesize the capillary macro argument.
5. `56ba5644` align example README/scripts with committed results (docs/tooling only).
6. `63260266` add the 2D `sigma_model=1` golden test (`E39EFF77`).

## Key files

- `src/simulation/m_surface_tension.fpp` — `c_sigma` (43), unconditional alloc (69), `sigma_face`
  (116/163/210), σ(T) fill loop (313-323)
- `src/simulation/include/inline_capillary.fpp` — the parametrized stress-tensor macro + paren fix
- `src/simulation/m_sim_helpers.fpp:48` — `f_compute_mixture_temperature` (reused, unchanged here)
- `src/common/m_global_parameters_common.fpp` — GPU_DECLARE (96), defaults (447)
- `src/simulation/m_start_up.fpp:1042` — GPU_UPDATE
- `toolchain/mfc/params/definitions.py:627`; `case_validator.py:755`
- `examples/2D_thermocapillary_migration/` — README, `case_Ma_20.py`, `measure.py:153`
- `toolchain/mfc/test/cases.py`; golden files under `tests/E39EFF77/`
