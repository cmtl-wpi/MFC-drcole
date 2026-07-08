# Chemistry-free bulk thermal conduction

Branch `pr/thermal-conduction`, the 15-commit foundation of the stack
(`f4ab1aab..origin/pr/thermal-conduction`, tip `6d36ba4d`). Line numbers cite that branch tip.

Adds an explicit Fourier heat-conduction term to the energy equation for non-reacting
multi-fluid flows. It also introduces the temperature-recovery helper and the flat-device-array
plumbing pattern that the μ(T) and σ(T) features reuse.

## 1. The physics

The added term is the divergence of the Fourier heat flux in the total-energy equation:

```
∂E/∂t + … = ∇·(k ∇T)
```

a **parabolic diffusion term** grafted onto the hyperbolic multi-fluid system. It is handled as a
**source flux** — the same machinery MFC already uses for viscous stress, surface tension, and
chemistry diffusion — not as part of the Riemann flux.

### Temperature recovery

Temperature is not a state variable; it is recovered algebraically from the mixture stiffened-gas
EOS in the new shared helper `f_compute_mixture_temperature`
(`src/simulation/m_sim_helpers.fpp:48`):

```fortran
gamma_mix  = Σ_i alpha_i * gammas(i)
pi_inf_mix = Σ_i alpha_i * pi_infs(i)
mCP        = Σ_i alpha_rho_i * cvs(i)*gs_min(i)          ! gs_min(i) = 1/gammas(i)+1 = γ_phys
T_cell     = ((gamma_mix + 1)*p + pi_inf_mix) / max(mCP, sgm_eps)
```

Because `gammas(i) = 1/(γ_phys − 1)`, `cvs(i)*gs_min(i) = c_v·γ = c_p`, so `mCP = ρ·c_p`. **This
`c_p` denominator (correct for recovering T) versus the `c_v` diffusivity that governs the update
is the crux of the stability-number fix (`490315e5`).** This helper is deliberately public and
shared: μ(T) and σ(T) both evaluate temperature through the same formula so the closures cannot
drift apart.

### Mixture conductivity and discretization

Face conductivity uses the **harmonic (series-resistance) mixture closure** (Samareh et al. 2014,
Eq. 8 — appropriate for heat crossing a diffuse interface), `m_thermal_conduction.fpp:189`:

```fortran
do i = 1, num_fluids
    k_L = k_L + clamp(alpha_i_L, 0, 1) / max(kappas(i), sgm_eps)     ! accumulate α/k
end do
k_L = 1 / max(k_L, sgm_eps)                                          ! reciprocal → harmonic mean
k_face = 0.5*(k_L + k_R)
```

Volume fractions clamped to `[0,1]`, reciprocals guarded with `sgm_eps` (α can overshoot near a
diffuse interface). `kappas(i)` is per-fluid `fluid_pp(i)%k_therm`.

The conductive flux at a face is a 2-point central difference divided by the **cell-center**
spacing (correct on stretched grids), `m_thermal_conduction.fpp:202`:

```fortran
dT_dxi = (T_tc(x+off) - T_tc(x)) / grid_spacing
flux_src_vf(eqn_idx%E)%sf(x,y,z) = flux_src_vf(eqn_idx%E)%sf(x,y,z) - k_face*dT_dxi
```

The RHS takes the divergence with a second first-difference (`m_rhs.fpp:1422`). Composing the two
yields the standard **FTCS** 3-point Laplacian `(T_{j-1} − 2T_j + T_{j+1})/dx²` — explicit in
time, second-order central in space.

### The conduction stability number (`490315e5`)

FTCS diffusion is conditionally stable: the diffusion number `d = α·Δt/Δx²` must satisfy
`2·(#dims)·d ≤ 1`. The commit makes the reported `tcfl` equal that dimensionless group, so `1.0`
is the true stability edge. It fixes **two** bugs at once:

1. **Wrong diffusivity.** The per-step diffusivity of `dE = Δt·∇·(k∇T)` is `α = k/(ρ·c_v)`, not
   `k/(ρ·c_p)`. The denominator changed from `mCP` (= ρc_p) to `mCV = Σ ρ_i c_v,i` (= ρc_v) — a
   factor γ larger.
2. **Missing FTCS normalization.** The reported number and the adaptive-`dt` limit gained the
   `2·#dims` factor (1D → `2`, 2D → `4`, 3D → `6`). The adaptive `dt` is the exact inverse of
   `tcfl`, so requesting `cfl_target` produces exactly that diffusion number.

The conduction diffusion number rides in the **existing VCFL field** (not a new column), so one
stability report covers both diffusive limits. Isothermal-wall rows are ~1.5× stiffer than the
interior (ghost reflection sharpens the near-wall stencil), so keep margin below 1.

## 2. How it threads the architecture

- **Params** (`toolchain/mfc/params/definitions.py`): global flag `thermal_conduction` (LOG);
  per-fluid `k_therm` (REAL) in the `fluid_pp(f)%` loop. Declarations/namelist auto-generate.
- **Derived type** (`src/common/m_derived_types.fpp:386`): `k_therm` on `physical_parameters`,
  next to `cv`. Defaults `k_therm = 0` (all three targets) and `thermal_conduction = .false.`
  (shared `s_assign_default_values_to_user_inputs`).
- **`kappas` array** (`src/common/m_variables_conversion.fpp`): mirrors `cvs` exactly —
  declared/`GPU_DECLARE`d, `@:ALLOCATE(kappas(1:num_fluids))`, populated
  `kappas(i)=fluid_pp(i)%k_therm`, pushed via `GPU_UPDATE`, `@:DEALLOCATE`d in both
  `#ifdef MFC_SIMULATION` arms. **This is the template μ(T)/σ(T) copy.**
- **The computation** — new module `src/simulation/m_thermal_conduction.fpp` (224 lines):
  `s_get_thermal_conduction` fills the ghost-extent work array `T_tc` once per RHS eval (ghost
  cells already populated, so no extra halo exchange), applies isothermal BCs;
  `s_compute_conductive_flux(idir,…)` accumulates the per-direction face flux; init/finalize
  `@:ALLOCATE`/`@:DEALLOCATE` `T_tc`.
- **RHS wiring** (`m_rhs.fpp`): `use m_thermal_conduction`; once-per-step T refresh (583); the
  per-direction flux inside the sweep (751-770); the additional-physics gate extended everywhere
  from `viscous .or. surface_tension .or. chem…` to also include `.or. thermal_conduction`; a
  dedicated conduction-only E-slot divergence block for when viscous/ST are both off (1422/1522/
  1614), skipped otherwise to avoid double-counting.
- **CFL/`dt` plumbing**: `s_compute_stability_from_dt`/`s_compute_dt_from_cfl` gain an optional
  `alpha_T`; `m_data_output.fpp`/`m_time_steppers.fpp` compute `alpha_T_cell` (harmonic k, `mCV`);
  VCFL guards extended to `viscous .or. thermal_conduction`; the `Rc` column stays viscous-only;
  `m_mpi_common.fpp` VCFL reduce likewise extended.
- **GPU** (`6bc31ade`): `thermal_conduction` is read *inside* the CFL acc routines, so it must be
  in an `acc declare`, and because it is a **runtime** flag it must be declared
  **unconditionally** (`m_global_parameters_common.fpp:101`):
  ```fortran
  $:GPU_DECLARE(create='[thermal_conduction]')
  ```
  Without it, any `--gpu acc` build linking those routines fails at nvlink "Undefined reference".
- **Init/finalize** (`8e99d866`, `m_start_up.fpp:839/1107`): gated
  `s_initialize/finalize_thermal_conduction_module`. Without init, `T_tc` is never
  device-allocated and the kernel writes a null device pointer — a step-1 CUDA illegal-address
  crash.
- **Boundary conditions**:
  - Isothermal (Dirichlet) walls: ghost reflection `T_ghost = 2·Twall − T_interior` so the face
    temperature equals `Twall`. Applied **only** when the per-rank BC code is negative
    (`bc < 0`, i.e. the rank owns the physical face); on interior ranks the "ghost" cells are
    valid cross-rank halo data and must not be overwritten.
  - `54e0d442` — allow isothermal BCs with `thermal_conduction`: bulk conduction is now a valid
    heat-conduction provider (previously only chemistry diffusion was), at any boundary type.
  - `3b40644b` — prohibit isothermal on periodic boundaries: under MPI a periodic axis is
    rewritten to a neighbor-rank id (`bc ≥ 0`), so the `bc < 0` guard would fire on 1 rank and
    not on many — a rank-count-dependent result. Rejected in both the Python validator and the
    Fortran checker (which runs before decomposition rewrites the codes).
- **Validator/checker**: `case_validator.py check_thermal_conduction` — if any `k_therm` set,
  require `thermal_conduction`; require `model_eqns ∈ {2,3}`, per-fluid `k_therm>0` and `cv>0`;
  prohibit chemistry, igr, cyl_coord, bubbles, hypo/hyperelasticity, mhd, relax, ib. Mirrored in
  `src/simulation/m_checker.fpp`.

## 3. Design decisions

- **Why "chemistry-free".** MFC already had heat conduction, but welded to the chemistry
  species-diffusion path (`chem_params%diffusion`, a single reacting gas with a Cantera-backed
  temperature). This provides bulk Fourier conduction for the plain multi-fluid stiffened-gas
  system, independent of chemistry — hence the hard `@:PROHIBIT(chemistry)`; the two are mutually
  exclusive providers.
- **`6e272a4b` "restore upstream q_prim_qp allocation block"** — a diff-discipline repair. The
  foundation commit had rewritten the `q_prim_qp` allocation loop in a way that carried an
  unrelated MHD (`hyper_cleaning`) fix. This reverts the block to the exact upstream form so the
  PR diff stays conduction-scoped; the MHD fix belongs in its own PR.
- **Precision**: all `wp`, generic intrinsics; `T_tc` and `kappas` are `real(wp)`; reciprocals
  guarded by `sgm_eps`; `k_therm>0`, `cv>0` enforced so no division by zero.

## 4. Commit-by-commit

1. `809696a5` test: refresh coverage map — no logic.
2. `7e6123a2` Add chemistry-free bulk thermal conduction — the foundation (module, helper, RHS
   wiring, `kappas`, params, derived-type member, validator + checker).
3. `046081d9` Add validation example cases — five 1D/2D/3D cases vs analytic heat-equation.
4. `e4402437`, `a59f49e2`, `c92c4843` — example diagram/geometry polish.
5. `54e0d442` allow isothermal BCs with thermal_conduction.
6. `6bc31ade` GPU_DECLARE thermal_conduction so the CFL acc routines device-link.
7. `8e99d866` wire init/finalize — fixes a step-1 CUDA illegal-address crash.
8. `6e272a4b` restore upstream q_prim_qp allocation block — back out an out-of-scope MHD change.
9. `490315e5` conduction stability number: cv-based diffusivity, FTCS-normalized.
10. `3b40644b` prohibit isothermal BCs on periodic boundaries.
11. `a53e54d5` docs: `fluid_pp%k_therm`, isothermal-BC section.
12. `c2b9818d` add 1D thermal conduction golden case (`B38D8D17`).
13. `6d36ba4d` examples: suite fixes and refreshed results.

## Key files

- `src/simulation/m_thermal_conduction.fpp` — the module
- `src/simulation/m_sim_helpers.fpp:48` — `f_compute_mixture_temperature` (the shared T helper)
- `src/simulation/m_rhs.fpp` — RHS wiring (583, 751-770, 1422/1522/1614)
- `src/simulation/m_data_output.fpp`, `m_time_steppers.fpp` — CFL report / adaptive dt
- `src/simulation/m_checker.fpp` — `s_check_inputs_thermal_conduction`
- `src/simulation/m_start_up.fpp` — init/finalize
- `src/common/m_variables_conversion.fpp` — `kappas` array
- `src/common/m_global_parameters_common.fpp:101` — unconditional GPU declare
- `src/common/m_derived_types.fpp:386` — `k_therm`
- `toolchain/mfc/case_validator.py` — `check_thermal_conduction`
- `examples/{1D,2D,3D}_thermal_conduction*/`; golden test `tests/B38D8D17/`
