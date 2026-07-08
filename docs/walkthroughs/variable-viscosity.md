# Temperature-dependent (Arrhenius) viscosity — μ(T) = exp(C + D/T)

Branch `pr/variable-viscosity`, 7 commits on top of the thermal-conduction base
(`origin/pr/thermal-conduction..origin/pr/variable-viscosity`, `6743f72e`..`3b22de6a`). Line
numbers cite that branch tip.

The feature adds a per-fluid Arrhenius viscosity model. It adds **no new stored field, no new
equation, and no new flux path**. It changes exactly one line of the existing viscous code: the
place where the HLLC solver reads a fluid's constant viscosity, it instead computes μ(T) from
the local temperature. Everything downstream — the stress tensor, the momentum and energy
fluxes — is untouched.

## 1. The physics

The model is a single exponential (the liquid / activation-energy law):

```
μ(T) = exp(C + D/T)
```

- **D** (`visc_d`) is the activation-energy term, `D = Eₐ/R`. With `D > 0`, viscosity **falls
  as temperature rises** — hot liquid is thinner.
- **C** (`visc_c`) is the log-prefactor: `exp(C)` is the viscosity in the `T → ∞` limit.

This is deliberately *not* Sutherland's law (a gas law, which increases with T); that is flagged
as separate future work in the example README.

The example inverts C and D from physical targets rather than setting them raw
(`examples/2D_Couette_Variable_Viscosity/couette_config.py`):

```python
D = np.log(mu_contrast) / (1.0/T0 - 1.0/T1)   # fix μ(T0)/μ(T1) across the gap
C = np.log(mu_ref) - D / T_ref                 # fix μ at the mid-gap reference temperature
```

## 2. The injection point — understand "the constant path" first

This is the crux of the feature, and the part that reads as confusing until you know what MFC
stores. **MFC never carries μ directly.** It carries, per fluid, a **Reynolds number** `Re` from
the case file. In the code's nondimensionalization the shear viscosity is `μ = 1/Re`.

At start-up that per-fluid `Re` is copied into a flat array `Res_gs`
(`src/simulation/m_riemann_solvers.fpp:74`):

```fortran
Res_gs(i, j) = fluid_pp(Re_idx(i, j))%Re(i)     ! Res_gs = Re = 1/μ
```

### The constant path (stock MFC, unchanged)

For every face, the viscous block loops over fluids and accumulates a mixture inverse-viscosity,
then reciprocates it (`m_riemann_solver_hllc.fpp:263`):

```fortran
Re_L(i) = qL_prim_rsx_vf(..., eqn_idx%E + Re_idx(i,q)) / Res_gs(i,q) + Re_L(i)
!         └── α_q (reconstructed volume fraction) ──┘   └ 1/μ_q ┘
!         term = α_q / (1/μ_q) = α_q · μ_q
...
Re_L(i) = 1._wp / max(Re_L(i), sgm_eps)     ! Σ α_q·μ_q  →  reciprocal = mixture Reynolds no.
```

So the accumulator builds `Σ α_q · μ_q` — the volume-fraction-weighted **arithmetic** mean of
the per-fluid viscosities (equivalently, the harmonic mean of the Reynolds numbers) — and the
final reciprocal turns it back into a mixture Reynolds number stored at
`Re_avg_rsx_vf(...,i) = 2/(1/Re_L + 1/Re_R)` (`m_riemann_solver_hllc.fpp:350`). **This loop —
dividing the reconstructed volume fraction by the baked-in `Res_gs` — is "the constant path".**
It has always been in MFC and it runs for every viscous case.

### The override (this feature)

For a fluid flagged `visc_model == 1`, the feature replaces the single division `/Res_gs(i,q)`
(divide by `1/μ` = multiply by the *constant* μ) with a multiply by the *temperature-dependent*
μ(T) (`m_riemann_solver_hllc.fpp:255`):

```fortran
! fluids with visc_model == 1; Res_gs stores the constant 1/mu for all other cases.
if (i == 1 .and. viscous_T_dependent .and. visc_models(Re_idx(i,q)) == 1) then
    Re_L(i) = qL_prim_rsx_vf(..., eqn_idx%E + Re_idx(i,q)) &
              * exp(min(visc_cs(Re_idx(i,q)) + visc_ds(Re_idx(i,q))/max(T_L, sgm_eps), &
                        visc_exp_arg_cap)) + Re_L(i)
    ! ... same for Re_R with T_R
else
    Re_L(i) = qL_prim_rsx_vf(..., eqn_idx%E + Re_idx(i,q)) / Res_gs(i,q) + Re_L(i)  ! constant path
end if
```

The term is `α_q · μ(T)_q` — the exact temperature-dependent analogue of the constant path's
`α_q · μ_q`. **Same accumulator, same reciprocation, same downstream stress tensor and fluxes.**
The only thing that changed is the value of one fluid's viscosity at one face.

Two guards on the branch:
- `i == 1` restricts the override to **shear** viscosity; bulk viscosity (`i == 2`) always uses
  the constant path.
- When no fluid uses model 1, `viscous_T_dependent` is `.false.`, the `if` is never taken, and
  the code runs the original constant path **bit-for-bit**. That is why `visc_model = 0` is
  bitwise-identical to upstream.

### Where T comes from

`T_L` / `T_R` are computed per face from the reconstructed left/right states, inlining the same
mixture stiffened-gas EOS formula as `f_compute_mixture_temperature`
(`m_riemann_solver_hllc.fpp:242`):

```fortran
mCP_L = Σ_i  qL_prim_rsx_vf(..., i) * cvs(i) * gs_min(i)
T_L   = ((gamma_L + 1._wp)*pres_L + pi_inf_L) / max(mCP_L, sgm_eps)
```

It cannot call the scalar-field helper because the Riemann solver works on the reconstructed
face states, not on `q_prim_vf` — so it replicates the formula inline. The requirement `cv > 0`
(checker, below) exists precisely because this divides by `mCP`. **The temperature machinery is
inherited wholesale from the thermal-conduction base.**

## 3. How it threads the architecture (the 7 layers)

1. **Params** (`toolchain/mfc/params/definitions.py:868`): per-fluid `visc_model` (INT),
   `visc_c`, `visc_d` (REAL). Descriptions in `descriptions.py:329`. Fortran declarations and
   namelist bindings auto-generate.
2. **Derived type** (`src/common/m_derived_types.fpp:387`): `visc_model`, `visc_c`, `visc_d`
   added to `physical_parameters`, right after `cv`/`k_therm` (the conduction pattern).
3. **Flat device arrays** (`src/common/m_variables_conversion.fpp`): `visc_models`, `visc_cs`,
   `visc_ds`, plus a single master flag `viscous_T_dependent` (`.true.` if *any* fluid uses model
   1). Allocated and populated in the species-to-mixture setup; `GPU_DECLARE`d and pushed with
   `GPU_UPDATE(device=...)`. The hot HLLC loop checks the one boolean before touching the arrays.
4. **Injection point** — the HLLC override above.
5. **GPU**: arrays device-declared; `mCP_L`/`mCP_R` added to the loop's `private=[...]` list;
   inner loops use `GPU_LOOP(parallelism='[seq]')` inside the existing `GPU_PARALLEL_LOOP`.
6. **Validator + checker** (`case_validator.py` `check_visc_model` + `m_checker.fpp`
   `s_check_inputs_visc_model`, kept equivalent). Model 1 requires: `viscous=T`; `cv > 0` for
   **every** fluid; `riemann_solver=2` (HLLC) and `model_eqns=3` (the only branch the override
   lives in); `Re(1) > 0` and set (else the fluid never enters the viscous loop and the override
   silently does nothing); both `visc_c` and `visc_d` set. **Prohibited:** chemistry, cyl_coord,
   non-Newtonian, and CFL-based time stepping (`cfl_adap_dt`/`cfl_const_dt` — "the viscous dt
   limit does not yet see μ(T)", an honest scope guard).
7. **Example + test** — below.

## 4. Gotchas and subtleties

- **The exponent clamp** `min(..., 50)` (`visc_exp_arg_cap = 50._wp`, commit `9c6e192b`). `exp()`
  overflows to `Inf` above arg ≈ 88 in single precision (`log(huge)`). A transient bad `T` early
  in a run could push `C + D/T` past that and produce `Inf`/`NaN` stresses that silently poison
  the field. Capping the *argument* at 50 keeps μ finite in every precision build. It is a
  numerical guard, not physics — 50 is far above any meaningful exponent.
- **The double-free bug** (commit `a20abe8f`). The first cut duplicated the pre-existing
  fluid-property `@:DEALLOCATE` line in both preprocessor arms, double-freeing `gammas`,
  `gs_min`, etc. — aborting **every** run at finalization. Fixed by removing the stale line; the
  golden results were regenerated afterward (`3b22de6a`).
- **Precision**: `visc_cs`/`visc_ds`, `T_L`, `Re_L` are all `wp`. No `wp`/`stp` mixing; the whole
  computation lives in the `wp` Riemann-solver working set. `exp`/`min`/`max` are generic.
- **Interaction with the Re path**: the override is *inside* the existing `Re_idx` accumulation,
  only for `i==1` and only for model-1 fluids. A model-1 fluid with unset `Re(1)` would never
  enter the loop, so the checker makes `Re(1) > 0` a hard error. Model-1 and constant-viscosity
  fluids can be freely mixed in one case.

## 5. Validation

**Example — `2D_Couette_Variable_Viscosity`.** Steady planar Couette flow between two no-slip
walls, bottom at `T0=300 K`, top sliding at `U` and held at `T1=400 K`, periodic in x. The exact
solution solves the coupled two-point BVP (`reference.py`, scipy `solve_bvp` to ~1e-10):

```
momentum:  d/dy[ μ(T) du/dy ] = 0                  (uniform shear stress)
energy:    d/dy[ k dT/dy ] + μ(T)(du/dy)² = 0      (conduction + viscous heating)
```

Because μ(T) falls with T under uniform shear stress, `du/dy` varies inversely with μ and the
velocity profile **curves**: mid-gap `u = 0.352 U` versus the constant-μ `0.5 U` — that ~30%
deficit *is* the μ(T) signal. MFC matches at observed order **2.003** (`summary.json`), T exact to
~1e-6. It runs with `thermal_conduction=T` + isothermal walls to sustain the wall temperature
gradient — which is why it sits on the conduction base. The target `T(y)` is imposed via a
`y`-varying density at uniform pressure (temperature is not a state variable).

**Golden test — `6D29F444`** (`toolchain/mfc/test/cases.py`, 1D only): `visc_model=1`,
`visc_c=5.0`, `visc_d=10.0`, `cv=1.0` layered on the standard 1D viscous case (model_eqns=3, HLLC,
50 steps). The initial temperature contrast moves the shear Re ~2.7×, so the golden file genuinely
exercises the temperature branch.

## 6. Commit-by-commit

1. **`6743f72e`** Add temperature-dependent (Arrhenius) viscosity — params, derived-type members,
   device arrays + `viscous_T_dependent`, the inline mixture-T and the μ(T) override, first-cut
   validator/checker. (Introduced the double-free fixed in #3.)
2. **`6c5ca279`** Add `2D_Couette_Variable_Viscosity` — the exact-BVP validation case.
3. **`a20abe8f`** Remove duplicated deallocation — the double-free aborting every run.
4. **`9c6e192b`** Clamp the Arrhenius exponent; tighten checker/validator — cap at 50; `cv>0` for
   all fluids; `Re(1)>0`; both coefficients set; prohibit cyl_coord, non-Newtonian, CFL dt.
5. **`948c4830`** Fix example run instructions and step counts.
6. **`5f333eb5`** Add 1D `visc_model=1` golden case (`6D29F444`).
7. **`3b22de6a`** Regenerate example results at the current tip — order 2.003, errors unchanged.

## Key files

- `src/simulation/m_riemann_solver_hllc.fpp` — the override (lines 232-269), exponent cap (69-71)
- `src/common/m_variables_conversion.fpp` — device arrays, setup, deallocation
- `src/common/m_derived_types.fpp:387` — `visc_model`/`visc_c`/`visc_d`
- `src/simulation/m_checker.fpp` — `s_check_inputs_visc_model`
- `toolchain/mfc/case_validator.py` — `check_visc_model` + `PHYSICS_DOCS`
- `toolchain/mfc/params/definitions.py:868` — parameter registrations
- `examples/2D_Couette_Variable_Viscosity/` — `couette_config.py`, `reference.py`, `summary.json`
- `toolchain/mfc/test/cases.py`; golden files under `tests/6D29F444/`
