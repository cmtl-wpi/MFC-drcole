I have verified all the doc's hook points against the source. Here is the review.

---

# Adversarial feasibility review: droplet-bouncing sub-grid film model in MFC

## Line-reference audit — essentially all correct

I read every cited location. The file:line references are unusually accurate; reviewers can trust them:

- `m_thinc.fpp`: storage `mthinc_nhat`/`mthinc_d` at 36–38 ✓; alloc 199–201 ✓; `s_compute_mthinc_normals` 207–291 ✓; normal built from `eqn_idx%adv%beg` at 242, 245–258 ✓; `s_thinc_compression` 295–406 ✓; the `%end`/`1-aC` closure writes at 335, 343, 345, 353, 355 (MTHINC branch) and 377, 384, 386, 393, 395 (THINC branch) ✓.
- `m_rhs.fpp` source block 828–862 ✓ (hosts `acoustic_source`, `bubbles_euler`, `bubbles_lagrange`, `chemistry`, `cont_damage`).
- `m_body_forces.fpp:79–137` ✓ (canonical mom+energy source).
- `m_surface_tension.fpp`: flux kernel from 74, `c_divs(num_dims+1)` used at 131, color grads 255–279 ✓.
- `case_validator.py:427` (`int_comp != 0 and num_fluids != 2`) ✓; `:772` already relaxed to `num_fluids not in (2,3)` with a comment ✓.
- `m_grid.f90:51` serial / `178` parallel `stretch_x` ✓.
- `m_global_parameters.fpp:887` `GPU_UPDATE(device='[int_comp, ic_eps, ic_beta]')` ✓; all `parallel_macros.fpp` macro names exist ✓.

Physics sizing (acoustic `dt`, step counts) checks out arithmetically.

---

## Prioritized problems

### P1 — The "markers stay distinct" premise is weaker than the doc admits, and the single color field is *actively attractive* in the film

This is the load-bearing physics question and the doc's framing is too optimistic.

In `model_eqns=2` (5-equation) MFC is **single-pressure, single-velocity**. Two identical-EOS liquid markers have *nothing* thermodynamically or kinematically distinguishing them: same EOS, same shared cell pressure, same shared cell velocity (zero relative slip). They remain distinct **only** as separate advected `α` fields *while they occupy different cells*. Once the resolved gas layer (`α₃`) thins below the numerical diffusion width (~2–3 cells ≈ D/100–D/70 at your resolution), markers 1 and 2 co-locate, and **nothing in the scheme un-mixes them** — mTHINC re-sharpens each marker toward its own `α=1` independently, their sum exceeds 1, and `mpp_lim` clipping (`m_riemann_solver_hllc.fpp:839`, the `num_fluids > 2` branch) silently redistributes/fuses them while conserving mass. So distinctness is not a property of the multi-marker representation; it is entirely contingent on the injected force holding a ≥2-cell resolved gas gap open at every instant. State this as the central risk, not S1's "confirm markers stay distinct."

Worse: I traced the color function. `eqn_idx%c` is a **single user-painted advected scalar** (`m_global_parameters_common.fpp:296`; filled from `patch_icpp%cf_val` at `m_assign_variables.fpp:545`), and surface tension acts on `∇c` (`m_surface_tension.fpp:255–279`). With the natural painting (liquid `c=1`, gas `c=0`) the film is a `1→0→1` trough in `c`. MFC's CSF capillary force acts to **minimize color-interface area — i.e. it pulls the two drop surfaces together to erase the gas trough.** So the resolved surface-tension model is not merely "unable to express opposing curvature" (the doc's neutral wording in §4.1b/R6c); **in the film it is a coalescence-promoting force that the injected lubrication source must overcome on top of inertia.** This materially raises the required injected force and makes the per-marker-σ refactor (R6c) look non-optional rather than an escalation path. Recommend promoting per-marker smoothed curvature into the v1 critical path, or at minimum redesign the color painting (e.g. `c` tracing only one drop) and prove the film capillary force is non-attractive before Phase 2.

**Fix to Phase-1 validation:** the current gate ("initialize two markers 1–2 cells apart and advect; pass if they stay distinct over 10⁴ steps") is too weak — with no relative motion and no compression driving them together it passes trivially and never exercises the failure mode. Replace it with: two markers **driven together at `U_r/2`** by initial velocity, force OFF, and measure whether the gas gap collapses and the markers auto-merge. That is the actual R9 decision gate and also the R1 experiment.

### P2 — Reynolds ODE integrated in the RHS will be triple-counted per RK step

The doc puts film-thickness (`film_h`) time-integration in "Phase 2b … each step" but hooks the whole film machinery into `m_rhs.fpp:828–862`. `s_compute_rhs` is called **once per Runge–Kutta substage** (3× for RK3), and `s_compute_mthinc_normals` is already called per-substage at `m_rhs.fpp:579`. A history-dependent ODE (`∂(ρh)/∂t = …`) integrated inside the RHS advances drainage 3× per physical step with the substage `dt`, corrupting `h(t)` and hence the `h_crit` gate and `F_lub∝1/h³`. Momentum/energy *sources* belong in the RHS (they get RK-weighted correctly); the **thickness state update and the rupture min-reduction must live once-per-step in the time-stepper** (`m_time_steppers.fpp`), not in `s_compute_rhs`. The doc does not specify this split and as written it is a correctness bug. Make the per-substage vs per-step placement explicit.

### P3 — Missing timestep constraint: the `1/h³` injected source is stiff

Section 2's entire `dt` argument is acoustic-CFL only. But an explicit body force `F_lub = 3πμ_g V_n a⁴/(2h³)` **diverges as `h→0`** and is added explicitly. Near the dimple/rim (`h ~ 30–300 nm`) this imposes its own explicit-stability limit that can be *tighter* than the acoustic CFL, and it will manifest as ringing/blow-up in the momentum RHS rather than a clean CFL violation. The energy work term `u·f_film` inherits the same `1/h³` stiffness. This directly answers "does the acoustic timestep bite once you add near-contact stretch": the binding constraint may not be acoustic at all near contact — it may be source stiffness. The design needs either an implicit/analytic sub-step for the drainage+force coupling (Stefan–Reynolds is analytic — exploit that) or a documented source-based `dt` limiter. Currently unaddressed.

### P4 — The N-fluid closure can drive the gas volume fraction negative — the real robustness crux, presented as clean

The doc's "running-sum closure" (`α_N = 1 − Σ_f α_f`, `cont%end = ρ_N·α_N`) is exactly where this breaks. In the film, both liquid markers are transitioning in the *same* cells and mTHINC sharpens each **independently** (each solves its own `d` in `f_mthinc_solve_d`). Two independently over-steepened tanh profiles can give reconstructed face values with `α₁+α₂ > 1`, so the closure fluid (gas) gets `α₃_face < 0` and `ρ₃α₃ < 0` — negative partial density fed to the Riemann solver right where you most need the gas film. `mpp_lim` will clip it, but clipping in the film = the numerical merge you are trying to prevent. This is not a mechanical index edit; it is the crux robustness problem and deserves top billing in §4.1a, with an explicit consistency strategy (e.g. renormalize the liquid markers so `Σα ≤ 1−α_min` before closure, and prove `α₃_face ≥ α_min`).

### P5 — MTHINC branch specifics the doc misattributes

Two concrete corrections for §4.1a so the implementer edits the right code:

- The doc says "the monotonicity test (358) and A/B/C THINC algebra (372–395) apply per-marker unchanged." Those lines are the **plain-THINC (`int_comp=1`) branch** (`else` at 357). The path you must use for 3-D oriented interfaces is **MTHINC (`int_comp=2`), lines 324–356**, which has **no `moncon` monotonicity test** and no A/B/C algebra — it uses `f_mthinc_face_average`. So "358/372–395 unchanged" is the wrong branch; the per-marker loop wraps 327–355.
- **Reshaped-vs-physical indexing (inherited hazard):** `mthinc_nhat`/`mthinc_d` are stored in physical `(j,k,l)` (computed in `s_compute_mthinc_normals` over `q_prim_qp`), but read inside `s_thinc_compression` using the **reshaped** loop indices, with only a bare comment "Map reshaped (j,k,l) to physical" at line 325 and no actual remap. For `recon_dir=1` reshaped==physical so it works; the per-marker extension must be validated in **all three reconstruction directions** (`m_weno.fpp:1461`, `m_muscl.fpp:226`) because the normal components are in physical x/y/z while the reconstruction axis is permuted. The doc's "applies unchanged" glosses over this; treat directional correctness as a test requirement.

### P6 — Body-force template copy trap (energy sign / RHS zeroing)

The doc says implement the source by "copying the `m_body_forces.fpp:79–137` template." That template **zeroes `rhs_vf(mom:E)` first** (lines 89–99) because body forces are assembled standalone. At the `828–862` hook the RHS **already holds the flux divergence**; a literal copy would wipe out all transport. The doc's snippet correctly uses `+=`, but call out explicitly: do **not** copy the zero-init loop. Also note the template reads `q_prim_vf`/`q_cons_vf`, whereas the hook neighborhood uses `q_prim_qp%vf`/`q_cons_qp%vf` (pressure lives at index `E` in **q_prim**, total energy at `E` in q_cons) — the doc's call signature gets this right but the "copy the template" instruction invites the wrong array.

### P7 — Under-specified / smaller items

- **`h_Est` ray-integral vs the "`c_divs` shortcut":** §4.2a proposes masking `β` from `c_divs(num_dims+1)` "both markers present and anti-parallel." But with a single color field (P1) there is only **one** `∇c`, not `∇α₁` and `∇α₂`. The "both markers' gradients anti-parallel" test therefore requires the **per-marker** `∇α_θ` (available from `q_prim_qp%vf(adv%beg+f-1)`), not `c_divs`. Fix the field the shortcut reads.
- **R5 band dropout is likely, not a tail risk.** At D/200 the max-deformation film (~0.3 µm) is <1 cell (1.5 µm), so `β` is 0–1 cells thick and the tangential-gradient stencils in the Reynolds solve (which need ≥2 in-plane cells) have no room. The "fall back to the ray integral" mitigation doesn't create resolution that isn't there. This tension between "stop at D/200 for `dt`" (R4) and "need in-plane stencil in the film" (R5) is not reconciled; it may force D/500 (4× the steps) purely to make the film solver well-posed, changing the cost story.
- **Grid stretching symmetry:** §5 requires `x_domain%beg = −x_domain%end` and tanh clustering. Note `dx = minval(...)` at `m_grid.f90:~68` sets the global `dt` from the finest (plane) cell, and `s_mpi_reduce_min` reduces it across ranks — consistent with R4, worth stating that clustering buys locality but not a larger `dt`.
- **Per-marker Helmholtz smoothing cost (R8):** a linear solve per marker per substage on GPU is not sketched at all (no solver named; no `GPU` strategy). For an explicit code that otherwise does only local stencils, introducing a global implicit solve per step is the single biggest performance and complexity unknown — larger than the doc's "30–50% overhead." Either specify a matrix-free Jacobi/multigrid smoother with iteration count, or drop to an explicit diffusion smoothing of a few iterations.

---

## Bottom line on the four questions asked

1. **mTHINC N-fluid scoping:** Correctly localized to `m_thinc.fpp` for the *mechanics*, and the callers (`m_weno`/`m_muscl`) pass full reshaped state so the fields are available. But it is **not** as clean as "indices only": the independent per-marker sharpening + running-sum closure is a genuine `α₃<0` robustness problem (P4), the doc points at the wrong branch/lines for MTHINC (P5), and directional (reshaped-index) correctness must be tested (P5).
2. **RHS hook:** Hook point is correct and all fields are in scope/device-resident there. Two real defects: history-dependent drainage must not live in the per-substage RHS (P2), and the body-force template's zero-init must not be copied (P6).
3. **Two identical-EOS fluids staying distinct:** **No, not intrinsically.** Single-pressure/single-velocity gives zero relative slip; distinctness survives only while a ≥2–3-cell resolved gas gap persists, and the single color-function capillary force actively *pulls the film shut* (P1). The whole design rests on the injected force defeating both inertia and the spurious capillary attraction — validate that first with a proper driven-together test, not passive advection.
4. **GPU offload:** All macros (`GPU_PARALLEL_LOOP`, `GPU_DECLARE`, `GPU_UPDATE`, `GPU_ROUTINE`, `GPU_LOOP`, `END_GPU_PARALLEL_LOOP`) exist as cited; per-marker inner loops need `GPU_LOOP(parallelism='[seq]')` as the doc says. The unaddressed GPU risk is the per-marker Helmholtz **implicit solve** (P7), not the source kernels. And the acoustic-CFL timestep story omits the stiff `1/h³` source limit (P3), which may bite before acoustics do near contact.

Most consequential fixes, in order: **P1** (redo the marker-distinctness premise and the color-force sign; fix the Phase-1 test), **P2** (move drainage integration out of the RHS), **P4** (closure positivity), **P3** (source-stiffness timestep).