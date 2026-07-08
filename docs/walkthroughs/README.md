# Thermal-Marangoni stack — implementation walkthroughs

Study notes for the three stacked features on the `pr/*` branches, written so a maintainer
can explain the code in review as if they wrote it. These are *implementation notes*, not
user-facing feature docs (those live in `docs/documentation/case.md`). Line numbers cite the
tip of each feature's own branch.

- [Thermal conduction](thermal-conduction.md) — `pr/thermal-conduction` (the 15-commit base)
- [Variable viscosity μ(T)](variable-viscosity.md) — `pr/variable-viscosity` (+7 commits)
- [Variable surface tension σ(T)](variable-sigma.md) — `pr/variable-sigma` (+6 commits)

## The one idea that unifies all three

MFC's state vector has **no temperature**. All three features hang off a single helper the
thermal-conduction base introduced — `f_compute_mixture_temperature`
(`src/simulation/m_sim_helpers.fpp`), which inverts the mixture stiffened-gas EOS:

```
T = ((Γ_mix + 1)·p + Π∞_mix) / mCP,    mCP = Σ αᵢρᵢ·cvᵢ·γᵢ   (= ρ·c_p)
```

Once you can recover `T` cheaply, per cell, anywhere, the recipe for each feature is identical:

> **Find the one place a constant physical coefficient is read. Replace it with `coeff(T)`
> evaluated from the local temperature. Change nothing downstream.**

That is why each diff is small, and why "feature off" is bit-for-bit identical to upstream:

| Feature | Constant coefficient replaced | Where |
|---|---|---|
| Thermal conduction | *(new term entirely)* `k∇T` in the energy flux | new `m_thermal_conduction` module |
| Viscosity | shear viscosity `μ` (as `1/Re`) | HLLC Reynolds-number accumulation |
| Surface tension | capillary coefficient `σ` | CSF stress tensor, face-averaged |

## The template every feature repeats (the 7 layers)

Each feature touches the same seven layers in the same order. Recite this and you can
reconstruct any of them:

| Layer | Where | Conduction | μ(T) | σ(T) |
|---|---|---|---|---|
| 1. Param registration | `params/definitions.py` | `thermal_conduction`, `k_therm` | `visc_model/c/d` | `sigma_model`, `sigma_T_ref`, `sigma_dTdT` |
| 2. Derived-type member | `m_derived_types.fpp` | `k_therm` | `visc_c/d/model` | *(none — plain scalars, auto-gen)* |
| 3. Flat device array | `m_variables_conversion.fpp` | `kappas` | `visc_cs/ds` + `viscous_T_dependent` | `c_sigma` field |
| 4. **The injection point** | `src/simulation/` | new conduction flux | HLLC shear-Re override | capillary `sigma_face` |
| 5. GPU mapping | `GPU_DECLARE`/`GPU_UPDATE` | declare the flag | declare arrays | declare the scalars |
| 6. Validator + checker | `case_validator.py` + `m_checker.fpp` | conduction constraints | Arrhenius constraints | σ(T) constraints |
| 7. Example + golden test | `examples/`, `test/cases.py` | 1D conduction `B38D8D17` | Couette `6D29F444` | capillary `E39EFF77` |

Layer 4 is where the physics lives; the rest is the *same plumbing* dressed differently. The
`kappas`-mirrors-`cvs` pattern from conduction is the template the other two copy.

## Recurring gotchas (the reviewer's greatest hits)

These four show up in all three — knowing them is knowing the stack:

1. **GPU device-mapping.** Any module variable read inside a device kernel must be in
   `GPU_DECLARE` and pushed with `GPU_UPDATE(device=...)`. Miss it → nvfortran nvlink
   "Undefined reference", *or worse* the device silently reads the stale default (e.g.
   `sigma_model=0`, disabling the feature on GPU only). Both conduction (`6bc31ade`) and σ(T)
   (`d79874fc`) have a commit that is exactly this fix.
2. **Bitwise-identical fallback.** Every feature's "off" path executes the original code
   verbatim — the design invariant the golden tests protect, and the first thing to assert in
   review.
3. **Precision.** All new arithmetic is `wp` with generic intrinsics (`exp`, not `dexp`); the
   transient work-fields (`T_tc`, `c_sigma`) are `wp`, not `stp`, because they are computed,
   not stored/I-O.
4. **Scope guards over silent approximation.** Where the physics isn't wired (CFL `dt` doesn't
   see μ(T); conduction can't coexist with chemistry; isothermal-on-periodic is
   rank-count-dependent), the features *reject the case* rather than run it subtly wrong.

## How the commits read

Same rhythm in each feature: **land the feature → add the validation example → fix the bug the
example exposed → harden the guards → pin it with a golden test.**

- **Conduction (15):** foundation (`7e6123a2`) → examples → two GPU-correctness fixes
  (`6bc31ade` declare, `8e99d866` init/finalize) → diff-discipline revert (`6e272a4b`) →
  stability-number correction (`490315e5`) → periodic-BC prohibition (`3b40644b`) → docs + test.
- **μ(T) (7):** feature (`6743f72e`) → example → bugfix (`a20abe8f` double-free) → hardening
  (`9c6e192b` clamp + constraints) → doc fix → golden test → regenerate results.
- **σ(T) (6):** feature (`8e06deac`) → example → validator hardening (`966dc60d`) → GPU + macro
  fix (`d79874fc`) → README alignment → golden test.
