# Extending MFC to Reproduce Droplet Bouncing via a Sub-Grid Gas-Film-Drainage Model

**Status:** Design proposal (v2 — revised to incorporate physics and MFC-implementation review)
**Author:** Lab assistant (MFC extension)
**Date:** 2026-06-30
**Audience:** MFC domain experts committed to the MFC codebase (compressible, structured-Cartesian, no AMR, diffuse-interface)

**Change log v1→v2 (load-bearing revisions):**
1. **Collapsed the two lubrication descriptions into one.** The film pressure that drains `h`, the pressure injected as a force, and the diagnostic pressure are now a **single dimpled field** from one coupled Reynolds solve on the band (was: an inconsistent flat-film Stefan–Reynolds force bolted onto a spatially-resolved Reynolds thickness equation). Stefan–Reynolds is demoted to a flat-plate sanity check. (Physics P1, P3.)
2. **Reframed what the injected force does.** It supplies missing sub-grid **load transmission / deceleration and prevents premature numerical merge** — it is *not* the rebound force. **Rebound is capillary** (resolved drop deformation) and the lubrication force actually *vanishes* at maximum compression and *opposes* separation. R1's decision metric changed accordingly. (Physics P2.)
3. **Rupture gate now tracks the rim minimum** `min_rim h(t) < h_crit` on the dimpled field, and the high-We re-merge (region II→III) is restated as a **rim-thinning** mechanism, not `t_drain ∝ a⁴`. (Physics P4.)
4. **Added stiffness handling and a magnitude guard.** The `1/h³` source is stiff and, uncorrected, over-predicts pressure by orders of magnitude; v1 self-limits via two-way feedback, a `ρU_r²` cap, and sub-cycling / semi-implicit treatment. The binding timestep near contact may be **source stiffness, not acoustics**. (Physics P5, MFC P3.)
5. **Added explicit sub-grid dissipation.** The bounce is lossy; the Poiseuille dissipation is deposited as internal energy, not returned as reversible work. Corrects the v1 "add `u·f`" instruction, which produces an artificially elastic bounce. (Physics P6.)
6. **Marker-distinctness premise rewritten as the central numerical risk**, including that MFC's single color-function capillary force is **coalescence-promoting** in the film. Added a symmetry-plane path that removes the marker machinery entirely for head-on rungs. (MFC P1; Physics "opportunity".)
7. **Reynolds thickness integration moved out of the per-substage RHS** into the once-per-step time-stepper (else it is triple-counted under RK3). (MFC P2.)
8. **N-fluid closure positivity** (`α₃ ≥ 0`) promoted to the crux robustness problem, with a renormalization strategy. Corrected the mTHINC branch/line attribution and flagged reshaped-vs-physical index correctness. (MFC P4, P5.)
9. Corrected `λ` (66 nm, not 100 nm), softened the "exact" velocity claim, re-sourced the S2 transition-`We` targets, fixed the body-force-template copy trap and the `c_divs`-vs-`∇α_θ` detection field, and specified an explicit (not implicit) smoother. (Physics P7, P8; MFC P5–P7.)

---

## 1. Goal and scope

### 1.1 Objective

Reproduce **droplet bouncing** — two liquid drops that approach, deform, and rebound **without ever touching**, cushioned by a thin gas film — inside MFC, a compressible structured-Cartesian diffuse-interface solver with **no AMR**. Chen & Yang (2020) reproduced this in Gerris by *resolving* the film to ~15 nm with octree AMR and balanced-force surface tension. That route is infeasible in MFC (Section 2). MFC's route is to **model the film sub-grid** — track its thickness and drainage below the smallest cell, gate coalescence on a physical rupture criterion, and keep the two approaching liquid bodies mechanically separated so the *resolved* capillary spring can rebound them.

One correction of emphasis up front, because it reorders everything downstream: **the resolved physics already contains the restoring force.** The rebound is stored as surface-tension energy in the deformed drop (deformation scale `~D`, well resolved at `D/200`) and released elastically. The sub-grid film's jobs are (i) to **transmit the collision load** so the drops decelerate and deform correctly and (ii) to **gate merging** via `h(t)` against a rupture thickness `h_crit`. The film force is dissipative and *zero at the turning point*; it does not, by itself, push the drops apart (Section 3.2). The whole design is organized around this.

### 1.2 Fluid system and dimensionless targets

Throughout, the physical system is **tetradecane droplets in nitrogen at 1 atm** (Chen & Yang's system): density ratio `ρ_l/ρ_g = 666`, viscosity ratio `μ_l/μ_g = 119`. Properties used for sizing: `ρ_l ≈ 763 kg/m³`, `σ ≈ 0.0265 N/m`, `μ_l ≈ 2.14×10⁻³ Pa·s`, `μ_g ≈ 1.8×10⁻⁵ Pa·s`, liquid sound speed `c_l ≈ 1300 m/s`, gas sound speed `c_g ≈ 349 m/s`, **N₂ mean free path `λ ≈ 66 nm` at 1 atm, 293 K** (the v1 value of 100 nm was ~1.5× too high and inflated every `Kn` and therefore every rarefaction correction and merge threshold — see Section 3.4, Physics P7).

Governing groups (Chen & Yang Eqs. 1–4), with `U_r` the relative (collision) velocity and `D` the initial drop diameter:

```
B  = X/D                     (impact parameter, X = centroid offset ⟂ motion)
We = ρ_l U_r² D / σ          (inertia vs surface tension)
Re = ρ_l U_r D / μ_l
Oh = μ_l/√(ρ_l σ D) = √We/Re  (liquid viscous dissipation)
τ_c = D/U_r                  (convective time)
```

> **Note on the velocity label (do not advertise as exact).** Chen & Yang's Fig. 6 lists `U = 0.97 m/s, D = 306 µm, We = 8.6, Re = 105.9`. Back-solving `U_r = √(We·σ/(ρ_l D))` with tetradecane properties gives `U_r ≈ 0.99 m/s`, i.e. the tabulated `0.97 m/s` is the **relative/collision** velocity `U_r`, with each drop at `U_r/2 ≈ 0.49 m/s`. We adopt `U_r` as the characteristic velocity in `We`, `Re`, `τ_c`. Caveat verified in review: `We=8.6` back-solves to `U≈0.99` but `Re=105.9` back-solves to `U≈0.97` — a **~2 % internal inconsistency**, well within property uncertainty but *not* an exact identity; the Eq. 10 / Eq. 11 self-consistency check in Section 5 is at the ~few-% level, not exact.

### 1.3 Staged scope (each stage gates the next)

| Stage | Target (Chen & Yang) | We | Re | B | Purpose |
|---|---|---|---|---|---|
| **S1** | **Fig. 6 head-on bounce** | 8.6 | 105.9 | 0 | First light: one clean head-on bounce; confirm the drops decelerate, deform, and rebound with the film gating merge. Runnable on the **symmetry-plane / single-drop path** (Section 3.9) — no two-marker machinery. |
| **S2** | **Nonmonotonic head-on transition** | sweep (see below) | — | 0 | Reproduce merge→bounce (low-We boundary) and bounce→merge (high-We boundary). Tests the drainage-time / rim-`h_min` / `h_crit` gate. Symmetry-plane path. |
| **S3** | **Oblique bounce** | ~30–50 | ~200–260 | **0.6** | 3-D, no symmetry plane; rotating film. This is the only stage that *requires* the general (Musehane) multi-marker film model over the axisymmetric (Mason) one. Chen & Yang's rendered off-center case is Fig. 10 (B=0.9, We=48.8) for qualitative comparison. |

> **S2 transition targets — re-sourced (Physics P8).** The commonly-quoted `We≈1.88` (merge→bounce) and `We≈10.3` (bounce→merge) are **Qian & Law's values for their specific gas/pressure condition**, and their central finding is that the bouncing window is a *strong function of ambient gas pressure*. Do **not** anchor S2 pass/fail to 1.88/10.3 unless the condition matches. **Read the transition `We` from Chen & Yang's own regime map at the tetradecane/N₂-1-atm condition** (their Fig. 2 / regime diagram) and cite that condition explicitly; use Qian & Law only where the gas/pressure matches. The *shape* of the requirement (a merge→bounce lower boundary and a bounce→merge upper boundary bracketing a bouncing window) is robust; the numeric `We` at each boundary is condition-specific and must be pinned to the target system before it becomes a pass criterion.

**Diagnostic quartet to reproduce** (Chen & Yang Figs. 13, 15): `h_center` (axis, max film thickness), **`h_min`** (rim of dimple, min film thickness — now the primary gate variable), `r_h` (dimple radius), and `p(0,0)` (gauge pressure at collision center), all vs `t/(D/U)`. The **dimpled** film (thicker at center, thinner at the rim), the **nonmonotonic `h_min`**, and the **rim pressure peak / negative rim gauge pressure during bounce-off** are the load-bearing physics — and they are exactly what a flat-film closure cannot represent, which is why v1 uses a coupled Reynolds solve (Section 3.3).

**Quantitative pass criteria:** bouncing-time correlation **Eq. 10, `τ_b/τ_c = 0.82√We` (valid 0.2 < We < 8.6)**; the **Fig. 16 pressure collapse** (film-center pressure ∝ `ρU_r²`); and the **restitution correlation Eq. 13** (`τ_b/τ_o = 0.29 We^(−0.35) + 0.88`), which is the direct test that the sub-grid **dissipation** is right (Section 3.5, Physics P6).

### 1.4 Out of scope (v1)

Shattering; unequal drop sizes; temperature-dependent properties; explicit van der Waals *body* force (we use a fixed `h_crit` cutoff instead, Section 3.6). Reflexive/stretching separation (regions IV/V) are natural follow-ons once merging works but are not v1 targets.

---

## 2. Why model the film instead of resolving it

### 2.1 The two levers Chen & Yang used, and that MFC lacks

Chen & Yang reached ~15 nm (`≈10⁻⁵ D`, level 15) using **two levers simultaneously**:

1. **Octree AMR** localizes the 15 nm cells to the film only — a typical 3-D run is ~6×10⁵ AMR cells, equivalent to ~5×10⁷ uniform cells. MFC is **structured-Cartesian with no AMR**, so it cannot localize fine cells: uniform 15 nm over the 9D×9D×9D domain is `(9·306µm / 15nm)³ ≈ (1.84×10⁵)³ ≈ 6×10¹⁵ cells` — impossible. Grid stretching (Section 4, `m_grid.f90`) can cluster cells at the collision plane, but every cell that reaches 15 nm still forces the global timestep down (lever 2), and `dx = minval(...)` with `s_mpi_reduce_min` sets the global `dt` from the finest (plane) cell regardless of where it sits (`m_grid.f90:~68`).

2. **Incompressible timestep.** Gerris advances on a capillary/viscous/advective CFL. MFC is **compressible and explicit**: `dt = CFL · dx / (|u| + c)`, limited by the **liquid acoustic speed** `c_l ≈ 1300 m/s` (larger than the gas speed, so liquid cells set `dt`).

### 2.2 The concrete step-count argument

At `dx = 15 nm`, CFL = 0.8:

```
dt_acoustic ≈ 0.8 · 15×10⁻⁹ / 1300 ≈ 9.2×10⁻¹² s   (9.2 ps)   [verified in review]
```

The full collision event is approach (close a `2D` surface gap at `U_r`) + bounce + separation:

```
t_approach ≈ 2D/U_r = 612µm / 0.99 ≈ 6.2×10⁻⁴ s
τ_b (We=8.6) ≈ 0.82√We · D/U_r ≈ 7.45×10⁻⁴ s
t_event ≈ 1.5×10⁻³ s  (approach + bounce + partial separation)
N_steps ≈ 1.5×10⁻³ / 9.2×10⁻¹² ≈ 1.6×10⁸ steps
```

Resolving to 15 nm costs **~10⁸ explicit acoustic steps** — *per the timestep alone*, before multiplying by the ~10¹⁵ uniform cells MFC would need without AMR. **The acoustic `dt` kills it even if the cell count were free**, and MFC has no AMR lever to reduce the cell count either. Both factors are structural to MFC; neither is a tuning knob.

There is a deeper mismatch: a compressible code pays the acoustic `dt` to resolve *fast liquid acoustics* that are **irrelevant** to a quasi-static lubrication film draining at `O(1 m/s)`. The film physics is essentially incompressible-lubrication; the solver is compressible-acoustic. Resolving the film forces the solver to spend its entire budget on the wrong physics. (Corollary flagged in Section 3.3: because the physical film pressure jump is only `ΔP ~ 10² Pa`, **gas compressibility in the film is negligible** — the compressible-conservative `ρh` form we adopt is harmless over-engineering, kept only for EOS consistency with the resolved gas.)

### 2.3 What the sub-grid approach buys

Stop refinement at the collision plane at **`D/200`–`D/500` (≈0.6–1.5 µm)** — enough to capture drop deformation and the resolved gas approach — and **model everything below that**. At `dx = D/200 = 1.53 µm`:

```
dt ≈ 0.8·1.53×10⁻⁶/1300 ≈ 9.4×10⁻¹⁰ s
N_steps ≈ 1.5×10⁻³ / 9.4×10⁻¹⁰ ≈ 1.6×10⁶ steps   (~10⁶, verified)
```

**Two decades fewer steps than resolving** (10⁶ vs 10⁸). **Caveat (MFC P3, added in v2):** this counts *acoustic* steps only. The injected `1/h³` source is stiff; near contact its explicit-stability limit can be **tighter than the acoustic CFL**, so the binding constraint on `dt` may not be acoustic at all in the last decade of `h`. Section 3.8 addresses this with sub-cycling / semi-implicit treatment and a source-based `dt` limiter; the `~10⁶` figure holds only if the film sub-system is integrated on its own (sub-cycled) clock rather than the global one.

The film **must** be sub-grid regardless of how far we refine: at maximum deformation Chen & Yang's film is `z/D ~ 10⁻³` (~0.3 µm at D=306 µm), already **below** `D/500`; at the rim just before rebound `h_min ~ 10⁻⁴ D ≈ 30 nm`. So across the entire load-bearing phase the film is thinner than the smallest affordable cell. Refinement to `D/200` locates the band and captures deformation; the dimple, rim, and drainage below ~0.6 µm are modeled. **This creates a genuine tension (Section 6, R5, now promoted): at `D/200` the film is 0–1 cells thick, but the in-plane Reynolds stencil needs ≥2 cells — which may force `D/500` purely to make the film solver well-posed, at 4× the step count.** This is reconciled in Section 5 (Rung 2), not hand-waved.

---

## 3. Physical model

### 3.1 Model selection: Musehane-general (multiple-marker) as the backbone — but only where needed

Two documented sub-grid architectures exist:

| | **Musehane / multiple-marker** (thesis 2016; JCP 2018) | **Mason / quasi-1D VOF** (CSIRO 2012) |
|---|---|---|
| Interfaces | N separate markers `α_θ`, each advected + smoothed; per-marker surface tension mechanically forbids merging | Single liquid marker; drops *can* merge numerically |
| Geometry | **General** — not restricted to head-on/equal drops; film detected by marker overlap | **Axisymmetric**, film pinned to an imposed `z=0` symmetry plane |
| Film pressure | Reynolds equation inverted on the band | Solved from reduced NS on the plane |
| Force coupling | One-way / diagnostic (his solver resolves the cushion) | Two-way body force `p_film n̂` |
| Rupture | `h < h_crit` fuses markers | `h < h_crit` (40 nm) switches a collision-plane ghost BC |
| Validation | Min film thickness for bouncing **matched thickness-based-AMR DNS** | Reproduced Pan et al. 2008 head-on tetradecane |

**Decision:** build the **general 3-D capability on Musehane's multiple-marker mechanics** (per-marker `α_θ`, geometric film detection, `h_crit` merge) **for S3**, because Mason's symmetry-plane assumption cannot represent the *rotating* off-center film (Chen & Yang Fig. 11) that defines region-II bouncing. Committing to Mason would forfeit S3.

**But note the review's key structural insight (Physics "opportunity"; MFC P1):** for **head-on S1/S2 you do not need the multi-marker machinery at all.** A reflective symmetry plane at the collision midplane physically holds the drops apart, enforces `V_n` symmetry for free, halves the domain, and — critically — **side-steps the entire "does the single color field hold two markers apart" risk** (R9). It also cleanly isolates whether the bounce is real physics (resolved capillary spring + film load transmission) versus an artifact of the marker trick. **v2 therefore stages the work so that the film physics is validated first on the symmetry-plane single-drop path (Section 3.9 / Phase 0), and the multi-marker machinery is only built for S3.** This is a major de-risking over v1, which coupled the two.

The multi-marker port still maps onto MFC's machinery — mTHINC is a per-marker extension away from N-fluid, and two liquid markers of the *same* material are exactly what MFC's `cont`/`adv` index ranges support — but with two hard caveats that v1 understated and Section 4.1 now treats as the crux (MFC P1, P4).

### 3.2 What the injected film force actually does (reframed — Physics P1, P2)

v1 framed the injected lubrication force as "the rebound force MFC must add to supply the repulsion." **This is conceptually wrong and is corrected here.**

- The lubrication pressure scales `p_film ∝ V_n` (normal approach velocity). At **maximum compression `V_n = 0`, so `p_film = 0`** — the film contributes *zero* force at the turning point. It cannot be the spring.
- During separation `V_n < 0`, so `p_film < 0`: the film **sucks the drops back together** (this *is* the physical Fig. 15b negative rim pressure). The lubrication force *opposes* rebound.
- The **actual restoring force is capillary**: surface-tension energy stored in the deformed drop, at deformation scale `~D`, which MFC **resolves** at `D/200`. Bouncing can therefore emerge in MFC — but because the capillary spring is resolved, *not* because we inject a repulsion.

So the injected force's real roles are: **(a) transmit the collision load** so the resolved drops decelerate and deform on the correct timescale (without it, an under-resolved `D/200` gas cushion passes too little load and the drops interpenetrate/merge numerically before deforming), and **(b) keep the resolved gas gap open** so the markers/interfaces do not numerically fuse before the capillary spring acts. It is a **load-transmission and anti-premature-merge** term, not a rebound term.

**Consequence for R1's decision gate:** the Phase-0 "force off" experiment is now judged by *"do the drops merge prematurely / is the deceleration and `h(t)` history correct?"* — **not** by *"is `τ_b` under-predicted / is rebound adequate?"* The rebound test is a *capillary* test and will (if the drop is well resolved) succeed with or without the force; the force's necessity shows up as **premature merge and wrong deceleration**, which is what to measure.

### 3.3 One coupled Reynolds solve — drainage, injected pressure, and diagnostic are the *same* field (Physics P1, P3)

v1's fatal inconsistency (Physics P1): it evolved `h` with a *spatially-resolved* Reynolds PDE (§3.3) but injected force from a *flat-film Stefan–Reynolds* formula (§3.4) whose pressure is **maximum at the center and zero at the rim** — the qualitative *opposite* of the dimpled film (rim minimum, rim pressure peak) that the validation targets demand. Two mutually inconsistent lubrication theories, neither able to produce the diagnostics. **v2 uses one description.**

On the interface band we solve a **single elliptic Reynolds equation for the film pressure `p_film`**, given the current sub-grid thickness field `h(x)` and the resolved kinematics (mean interface velocity `Q` and normal approach `V_n`):

```
∇_t · [ (ρ h³ / 12 μ_eff) ∇_t p_film ]  =  ∂(ρh)/∂t|_kinematic  +  Q · ∇_t(ρh)
    with   ∂(ρh)/∂t|_kinematic = −ρ V_n(x)            (rigid approach imposed by resolved interfaces)
    BC:    p_film → p_resolved,gas   at the band edge   (σ-jump-subtracted, smoothed — see below)
    ∇_t(·) ≡ ∇(·) − (∇(·)·n) n        (in-plane gradient),   n = unit interface normal
    Q = (u₁ + u₂)/2                    (mean interface velocity)
    μ_eff = μ_g / f(Kn)                (rarefaction correction, Section 3.4)
```

This is Musehane's coordinate-free Reynolds equation (his Eq. 3.44) written in **`ρh` conservative form** (compressible-consistent; gas compressibility itself is negligible here per Section 2.2, but the form keeps the EOS consistent). The **dimple emerges naturally**: gas squeezed from the center must travel farther to escape than gas at the rim, so pressure builds toward the center, the center resists approach more and thins slower, and the film ends up thicker at the center with `h_min` at the rim — exactly the observed structure.

The **same `p_film`** then does three jobs, eliminating v1's underdetermination (Physics P3):

1. **Drains the thickness field** (integrated once per step, in the time-stepper — see MFC P2 / Section 4):
   ```
   ∂(ρh)/∂t  =  −ρ V_n  −  Q·∇_t(ρh)  +  ∇_t·[ (ρh³/12 μ_eff) ∇_t p_film ]
   ```
2. **Is injected as the momentum source** `f_film = p_film n̂_θ` (Section 3.7 fixes the sign), producing the *dimpled* pressure profile at the resolved interface.
3. **Is the diagnostic pressure** compared to Chen & Yang's `p(0,0)` and Fig. 15/16.

Because all three are one object, the drainage driver is no longer ambiguous. The **resolved pressure enters only as the outer boundary condition** for the elliptic solve (subtract the surface-tension jump and smooth it before use, Musehane Eq. 3.65, to avoid feeding resolved acoustic noise into the film — this is the surviving, correctly-scoped remnant of v1's R6b).

**Head-on axisymmetric reduction (S1/S2 sanity check).** With the symmetry-plane path (Section 3.9) the solve collapses to the classical 1-D radial form, which is a cheap, well-posed check of the full solver:
```
(1/r) ∂/∂r [ r · (h³/12 μ_eff) · ∂p_film/∂r ]  =  −V_n
∂h/∂t = −V_n + (1/r) ∂/∂r [ r · (h³/12 μ_eff) · ∂p_film/∂r ]
```

**Flat-plate Stefan–Reynolds — demoted to a sanity check only.** For a *flat* parallel film of lamella radius `a` and gap `h` squeezed at `V_n` with immobile interfaces, the analytic result
```
p_flat(r) = 3 μ_eff V_n (a² − r²)/h³ ,   F_lub = 3π μ_eff V_n a⁴/(2 h³) ,   t_drain ∝ μ_eff a⁴/(F h_crit²)
```
is retained **only** to (i) unit-test the elliptic solver in the flat limit and (ii) reason about scalings (`F ∝ 1/h³`, `t_drain ∝ a⁴`). It is **never injected as the force** — it cannot produce a dimple, a rim minimum, or a rim pressure peak, and would fail Rungs 4–5 by construction.

Van der Waals disjoining pressure is **omitted from the drainage driver** (following Musehane/Malavé et al.: London–vdW forces do not aid dimple formation); vdW enters only through the fixed `h_crit` rupture cutoff (Section 3.6). This is a standard, defensible modeling choice (Mason, Kwakkel); the caveat is that the *merge outcome* is sensitive to `h_crit` below ~50 nm (R3).

### 3.4 Rarefaction correction (μ_eff) — not a harmless placeholder (Physics P7)

At `h → tens of nm`, `Kn = λ/h → O(1)` and continuum lubrication over-predicts resistance; the correction enters through `μ_eff = μ_g/f(Kn)` in the elliptic solve, the drainage flux, and the injected pressure. With the corrected `λ ≈ 66 nm`:
- at the dimple, `h ~ 0.3 µm`: `Kn ≈ 0.22`, first-order slip `f = 1+6Kn ≈ 2.3` (already marginal);
- at rupture, `h ≈ h_crit = 30 nm`: `Kn ≈ 2.2`, `f ≈ 14` — deep in transition/free-molecular, where first-order slip is **meaningless**.

Because `f(Kn)` rescales `μ_g` in *both* the force and the drainage flux, **it directly sets how fast `h` reaches `h_crit`, i.e. it directly sets the S2 transition `We`.** So `f(Kn)=1+6Kn` is **not** a benign placeholder — S2 outcomes are a strong function of this unknown function. **v2 requirements:**
1. Obtain the actual `Q(Kn)`/`f(Kn)` (Zhang & Law 2011, BGK-derived; or the 2025 rarefied dual-VOF / augmented-vdW papers) **before running S2**.
2. **Report S2 sensitivity to `f(Kn)`** as a first-class result, not an afterthought.
3. Use it only for gating; the *load-bearing* pressure is built at the dimple (`Kn ≈ 0.22`), where continuum-plus-slip is defensible, so the force is on firmer ground than the final-rupture instant (this is R2's mitigation, still valid).

### 3.5 Sub-grid dissipation — the bounce is lossy (Physics P6)

The physical lubrication resistance **dissipates energy inside the sub-grid film**; that dissipation is why bouncing is lossy and why `τ_b/τ_o` deviates from unity (Eq. 13). v1's instruction — "add the reversible work `u·f_film` to the energy RHS" — makes the force **reversible/elastic**: it conserves total energy but generates *no heat*, yielding an artificially elastic bounce (restitution too high, `τ_b` too long). That is exactly wrong for a dissipative film.

**Correct energy treatment (v2).** Split the coupling:
- **Reversible load transmission:** the momentum source `f_film` with its matching work term transmits the collision load between the approaching bodies (or, on the symmetry-plane path, reflects it) — this is real and must be present, or the drops don't decelerate correctly.
- **Irreversible Poiseuille dissipation:** the sub-grid viscous loss, per unit film area,
  ```
  Φ = (ρ h³ / 12 μ_eff) |∇_t p_film|²          (Poiseuille dissipation)
  ```
  must be **removed from the resolved liquid kinetic energy and deposited as internal energy in the band gas cells** — it does not return to the drop. Adding only the reversible term over-predicts restitution; omitting the reversible term breaks load transmission. v1 implements both and **validates net restitution against Eq. 13** (`τ_b/τ_o = 0.29 We^(−0.35) + 0.88`) — restitution is the *direct* test that the dissipation partition is right.

**Honesty flag (R6/R-diss):** the exact *conservative* discretization of "remove KE from resolved liquid, deposit as gas internal energy" in MFC's single-velocity mixture (where liquid and gas share a cell velocity, so the partition is KE↔mixture-internal-energy) is genuinely subtle and is an explicit design risk. v1 target: enforce global energy conservation to machine precision and let Eq. 13 arbitrate the physical partition. Getting this term wrong shows up as the wrong restitution, so it is observable.

### 3.6 Rupture criterion and the coalescence gate — track the *rim* (Physics P4)

**Critical thickness `h_crit`** (calibratable, default **30 nm**): Chen & Yang's vdW effective distance is 30 nm; Mason uses 40 nm; Deblais et al. (2024) interferometry measures jump-to-contact at ≈30 nm; Musehane reaches `O(10⁻⁷ m)`.

**Gate on the rim minimum, on the dimpled field:** merge when `min_rim h(t) < h_crit` (the minimum over the film band, which physically sits at the *rim* of the dimple), evaluated on the coupled-solve dimpled `h` — **not** on a single scalar `h`. When the gate trips, fuse the markers (Section 4, Phase 3); otherwise they stay distinct and the drops **bounce**.

This is the threshold form of **Kwakkel's `t_c ≥ t_drain` gate** (if `h` reaches `h_crit` while the drops are still pressed, they merge; if they separate before, they bounce); we implement the threshold form because `h(x,t)` is tracked directly.

**The We-nonmonotonicity mechanism — corrected (Physics P4).** v1 claimed *both* S2 transitions "emerge from `t_drain ∝ a⁴`." That is wrong for the high-We branch:
- **Low-We merge→bounce:** correct as stated — higher We → larger lamella `a` → `t_drain ∝ a⁴` grows → the film survives → bounce.
- **High-We bounce→merge:** does **not** follow from `t_drain ∝ a⁴` (that scaling predicts drops stay *bounced* as We rises — the opposite). The real mechanism is that increased inertia drives the **rim `h_min` below `h_crit` while the drops are still pressed**: the rim thins faster than the wider lamella can compensate. This is a competition of **`h_min(t)` (at the rim) against `h_crit`**, i.e. precisely the dimple/rim structure the demoted flat-film closure could not represent. State this explicitly in reporting S2, and **confirm the coupled solver actually produces a rim minimum before claiming S2** (a solver that only yields a center-max pressure would fail the high-We branch by construction).

### 3.7 Sign convention of the injected force (Physics P9)

`f_film = p_film n̂` is **repulsive only if `n̂` is the liquid-inward normal.** Define per marker:
```
n̂_θ = ∇α_θ / |∇α_θ|          (points into liquid θ)
```
so `p_film > 0` pushes the interface *into* its own drop, away from the film. With the gas-film-outward normal the sign flips and the force would be *attractive* and collapse the drops. Because the two markers have **opposite normals across the gap**, the two injected forces must be verified **anti-parallel** (drops pushed apart) in a unit test before Phase 2. On the symmetry-plane path there is one interface and one normal (into the single drop), which removes the ambiguity — another reason to validate there first.

### 3.8 Magnitude and stiffness — self-limiting and a hard guard (Physics P5; MFC P3)

A naive injection over-predicts by orders of magnitude. Plugging nominal numbers into the flat estimate at the dimple (`μ_eff≈2×μ_g`, `V_n=0.5`, `a≈90 µm`, `h=0.3 µm`) gives **~8 MPa** — four orders above the physical film pressure, which the Fig. 16 target says is `~ρU_r² ≈ 748 Pa = O(10²–10³ Pa)`. The mismatch is diagnostic: the injected closure scales **viscously** (`μ_eff/h³`) while the physical pressure is **inertially capped at `~ρU²`** by drop deformation. The resolution is that by the time `h ~ 0.3 µm` **the drop has already decelerated, so the instantaneous resolved `V_n ≪ U_r/2`** — the pressure never reaches 8 MPa physically. **Mitigations (all v1):**

- **(a) Two-way self-limiting.** Use the **resolved instantaneous `V_n`** (from the resolved interface velocities), so that injecting `f_film` decelerates the drop → `V_n → 0` → `p_film → 0`. The feedback caps the pressure automatically at the inertial scale — *provided the coupling is tight enough within a step*.
- **(b) Runtime magnitude guard.** Enforce `peak p_film ≲ C·ρU_r²` (`C = O(1–10)`) as a **stability guard**, not merely a validation output. This *is* Rung 5 promoted to a runtime assertion: if the injected pressure does not collapse on `ρU_r²`, the run is unphysical and should trip, not silently blow up.
- **(c) Stiffness control.** `p_film ∝ 1/h³` spanning two decades of `h` (1.5 µm → 30 nm) is a genuinely **stiff explicit source**; dumping ~MPa into a cell in one `~9×10⁻¹⁰ s` step blows up the momentum/energy RHS. The film sub-system (elliptic solve + drainage + force) is therefore **sub-cycled on its own stability-limited `dt_film`** and/or the `1/h³` term is treated **semi-implicitly** (Stefan–Reynolds is analytic — exploit it for the sub-step). A **source-based `dt` limiter** is added alongside the acoustic CFL. **Near contact the binding constraint may be source stiffness, not acoustics** — the Section 2 step-count budget assumes the film is sub-cycled off the global clock; if it is not, the global `dt` collapses and the `~10⁶`-step story fails.

### 3.9 The symmetry-plane path for head-on (S1/S2) — recommended, and it de-risks the marker machinery

For head-on collisions the collision midplane is a plane of symmetry: zero normal velocity, free tangential velocity — exactly a **reflective (free-slip) symmetry BC**. Simulate **one drop in the half-domain**; the "film" is the symmetric squeeze film between the drop and its mirror image (thickness `2×` the modeled `h`), and the free-slip plane correctly reproduces the *mobile-interface* film of two identical drops. Benefits:
- **Halves the cost** and enforces `V_n` symmetry for free.
- **Removes the two-marker distinctness problem entirely** for S1/S2 — there is one drop, one marker, one color normal, so R9 and the coalescence-promoting color force (Section 4.1, MFC P1) do not arise on this path.
- **Isolates the physics from the numerics:** it directly tests whether the bounce comes from the resolved capillary spring + film load transmission, independent of the marker trick.

**v2 staging consequence:** Phase 0 (Section 4) validates the *film physics* on this single-drop path; the *multi-marker machinery* (Phases 1–3) is built and validated only for S3, where there is no symmetry plane. This is a strictly better ordering than v1.

---

## 4. MFC implementation plan (phased, with concrete hooks)

All file:line references are to the current MFC tree and were **verified against source in review** (the audit found them accurate). Every new device array needs `GPU_DECLARE`; every new kernel loop needs `GPU_PARALLEL_LOOP`/`END_GPU_PARALLEL_LOOP`; every new case-file scalar needs `GPU_UPDATE(device=...)` (macros in `src/common/include/parallel_macros.fpp`). New references introduced by the review are cited inline.

### Phase 0 — Symmetry-plane head-on film physics (de-risk first; no markers)

Single drop in a half-domain with a reflective symmetry BC at the collision plane, `num_fluids = 2` (tetradecane drop + N₂), plus the film sub-system (Sections 3.3–3.8). This runs **S1 and S2** and validates the coupled Reynolds solve, the dimple/rim diagnostics, the dissipation/restitution, and the stiffness handling **without any of the multi-marker risk**. Only after Phase 0 passes do the multi-marker phases below (needed solely for S3) become worth their cost. Phase 0 also *is* the R1 experiment: run S1 with the film source **off** and confirm the drop merges into the plane prematurely / decelerates wrongly, then on and confirm correct `h(t)` and gating.

### Phase 1 — N-fluid interface sharpening (for S3): the two hard problems first

**Setup:** `num_fluids = 3` in a 5-equation model: `fluid 1 = tetradecane (drop-A marker)`, `fluid 2 = tetradecane (drop-B marker, identical EOS)`, `fluid 3 = nitrogen (continuous phase, residual)`.

**Reality check the crux — two identical-EOS markers are *not* intrinsically distinct (MFC P1).** In `model_eqns=2` (5-equation) MFC is **single-pressure, single-velocity**. Two identical-EOS liquid markers share the cell's pressure and velocity (zero relative slip) and are distinguished *only* by being separate advected `α` fields *while they occupy different cells*. Once the resolved gas layer `α₃` thins below the numerical diffusion width (~2–3 cells ≈ `D/100`–`D/70` at these resolutions), markers 1 and 2 co-locate and **nothing in the scheme un-mixes them**: mTHINC re-sharpens each marker toward its own `α=1` independently, their sum exceeds 1, and `mpp_lim` clipping in the `num_fluids > 2` branch (`src/simulation/m_riemann_solver_hllc.fpp:839`) silently redistributes/fuses them while conserving mass. **Distinctness is therefore contingent on the injected force holding a ≥2-cell resolved gas gap open at every instant** — it is not a property of the representation. This is the central S3 risk; state it as such.

**Worse — the single color-function force is coalescence-promoting in the film (MFC P1).** `eqn_idx%c` is a **single user-painted advected scalar** (`src/common/m_global_parameters_common.fpp:296`, filled from `patch_icpp%cf_val` at `src/pre_process/m_assign_variables.fpp:545`), and surface tension acts on `∇c` (`src/simulation/m_surface_tension.fpp:255–279`). With the natural painting (liquid `c=1`, gas `c=0`) the film is a `1→0→1` trough in `c`, and MFC's CSF capillary force acts to **minimize color-interface area — i.e. it pulls the two drop surfaces together to erase the gas trough.** So the resolved surface-tension model is not merely "unable to express opposing curvature" (v1's neutral wording); **in the film it is an *attractive*, coalescence-promoting force the injected lubrication source must overcome on top of inertia.** This materially raises the required injected force and makes the per-marker-σ refactor **non-optional for S3**, not an escalation path.

**Decision for S3 (revised):** promote **per-marker smoothed curvature** `κ_θˢ = −∇·(∇α_θˢ/|∇α_θˢ|)` (Musehane Eq. 3.7) into the S3 critical path, **or** redesign the color painting so `c` traces only one drop and **prove the film capillary force is non-attractive** before relying on it. Do not assume the single color field holds the markers apart.

**1a. Extend mTHINC to N markers — `src/simulation/m_thinc.fpp` — and guarantee closure positivity (MFC P4, P5).**

The module is hardcoded 2-fluid: one normal `mthinc_nhat(1:3,i,j,k)` and one interface position `mthinc_d(i,j,k)` (36–38), normal from the *single* first marker `eqn_idx%adv%beg` (242, 245–258), and four face fields assuming the complement of fluid-1 *is* fluid-2.

- **Storage:** promote `mthinc_nhat → mthinc_nhat(1:3, 1:num_fluids-1, i,j,k)`, `mthinc_d → mthinc_d(1:num_fluids-1, i,j,k)` (one normal + `d` per liquid marker; last fluid is the closure complement). Update allocations (199–204) and the `GPU_DECLARE` at 38.
- **`s_compute_mthinc_normals` (207–291):** wrap in `do f = 1, num_fluids-1`, `eqn_idx%adv%beg → eqn_idx%adv%beg + f - 1`, store into `mthinc_nhat(:,f,…)`, `mthinc_d(f,…)`.
- **`s_thinc_compression` — use the MTHINC branch, not plain-THINC (MFC P5, correction).** For 3-D oriented interfaces the required path is the **MTHINC (`int_comp=2`) branch, lines 324–356**, which uses `f_mthinc_face_average` and has **no** `moncon` monotonicity test and **no** A/B/C THINC algebra. v1's claim that "the monotonicity test (358) and A/B/C algebra (372–395) apply per-marker unchanged" points at the **plain-THINC (`int_comp=1`) branch** (`else` at 357) — the wrong branch. The **per-marker loop wraps 327–355**. For each marker compute its own `aTHINC` from its own normal/`d` and write `adv%beg+f-1`, `cont%beg+f-1` (density `rho_f = v_rs_ws(…,cont%beg+f-1)/alpha_f`). Decorate the per-marker inner loop with `GPU_LOOP(parallelism='[seq]')`.
- **Closure with an enforced positivity floor (MFC P4 — the real robustness crux, not a mechanical edit).** The naive running-sum closure `α_N = 1 − Σ_f α_f` is exactly where the scheme breaks: both liquid markers transition in the *same* film cells, mTHINC sharpens each independently (each solving its own `d` in `f_mthinc_solve_d`), and two over-steepened tanh profiles can reconstruct face values with `α₁+α₂ > 1` → `α₃_face < 0` → `ρ₃α₃ < 0` fed to the Riemann solver **right where the gas film matters most**; `mpp_lim` then clips it, and clipping in the film *is* the numerical merge we are trying to prevent. **Required strategy:** before forming the closure, **renormalize the liquid markers so `Σ_f α_f ≤ 1 − α_min`** (small floor), then set `α_N = 1 − Σ_f α_f ≥ α_min`, `cont%end = ρ_N·α_N`. **Prove `α₃_face ≥ α_min` at every reconstructed face** as a unit test. Give this top billing in the module changes.
- **Reshaped-vs-physical indexing must be tested in all three directions (MFC P5).** `mthinc_nhat`/`mthinc_d` are computed in **physical** `(j,k,l)` (over `q_prim_qp`) but read inside `s_thinc_compression` using the **reshaped** loop indices, with only a bare comment at line 325 and no actual remap. For `recon_dir=1` reshaped == physical so it works; the per-marker extension must be validated in **all three reconstruction directions** — callers `src/simulation/m_weno.fpp:1461` and `src/simulation/m_muscl.fpp:226` permute the reconstruction axis while the normal components stay in physical x/y/z. Treat directional correctness as an explicit test requirement, not an "applies unchanged."

**1b. Surface tension — per-marker curvature is now on the S3 critical path (see the crux above).** `m_surface_tension.fpp` remains `num_fluids`-agnostic for the *outer* interfaces (capillary force from the single color field, grads 255–279, flux 74–131), and that is fine away from the film. But because the single color field is *attractive* in the film (MFC P1), S3 requires either per-marker smoothed curvature (`s_get_capillary` 236–327 refactored to per-`α_θˢ`) or a redesigned color painting with a proven non-attractive film force. This is the larger surface-tension refactor; unlike v1 it is **not deferred** for S3 — it is a prerequisite. (S1/S2 avoid it entirely via the symmetry plane, Section 3.9.)

**1c. Relax validator guards — `toolchain/mfc/case_validator.py`.**
- Line 427: `int_comp != 0 and num_fluids != 2` → allow the N-fluid range (`num_fluids not in (2,3)` or a configured max). This is the guard that blocks N-fluid THINC.
- Line 772: surface-tension guard already relaxed to `num_fluids not in (2,3)` with a comment noting THINC was 2-fluid-only — Phase 1a removes that limitation, so update the comment.

**Phase 1 validation (revised — the passive test is too weak, MFC P1).** v1's "initialize two markers 1–2 cells apart and advect; pass if they stay distinct over 10⁴ steps" passes **trivially** with no relative motion and never exercises the failure mode. **Replace it with the real R9/R1 gate:** initialize two markers **driven together at `U_r/2`** by initial velocity, **film force OFF**, and measure whether the resolved gas gap collapses and the markers auto-merge (and whether `α₃` goes negative / gets clipped). Pass = the gap survives and markers stay distinct *only if* the force+per-marker-σ machinery is doing its job; a merge here quantifies exactly how much injected force / curvature is needed.

### Phase 2 — Film detection + thickness tracking + coupled solve + source injection

**2a. Film detection — use per-marker `∇α_θ`, not the single `c_divs` (MFC P7).** For each marker compute a smoothed field `α_θˢ` (Section 2b) and detect the thin-film band:
```
β = pos(ξ·δh − h_Est) · pos(α₁ˢ − ε) · pos(α₂ˢ − ε),   ξ ≈ 3–5,  ε ≈ 10⁻³
h_Est = h_Est1 + h_Est2  from the ray integrals  ∇h_Estθ·n_sθ = (1−α_θ)  (Musehane Eqs. 3.46–3.58)
```
v1's "MFC shortcut" proposed masking `β` from `c_divs(num_dims+1)` with a "both markers anti-parallel" test — but there is only **one** `∇c`, not `∇α₁` and `∇α₂`. The shortcut must instead read the **per-marker gradients** `∇α_θ` from `q_prim_qp%vf(adv%beg+f-1)` (both present within `ξ` cells and anti-parallel), falling back to the full ray integral if too coarse. On the symmetry-plane path (S1/S2) detection is trivial (band = cells within `ξ·δh` of the plane), which is a further reason to validate the solver there first.

**2b. Smoothing solver — explicit, not implicit (MFC P7 / R8).** The per-marker smoothing `α_θˢ − α_θ = ∇²(η α_θˢ)` (Musehane Eq. 3.3) is written as a Helmholtz solve, but a **global implicit solve per marker per substage on GPU is the single largest performance/complexity unknown** in the design — larger than v1's "30–50 % overhead" estimate, and out of character for an otherwise local-stencil explicit code. **v2 specifies an explicit smoother:** a **few (3–5) fixed Jacobi / explicit-diffusion iterations** with a documented iteration count, matrix-free, decorated with `GPU_PARALLEL_LOOP`. If accuracy demands more, escalate to matrix-free multigrid — but do not introduce a global implicit solve without accounting for it as the dominant cost.

**2c. Sub-grid thickness field and the coupled solve — integrate ONCE PER STEP in the time-stepper (MFC P2).** Add a device field `film_h(i,j,k)` (`GPU_DECLARE`). **Critical placement correction:** the thickness field is a **history-dependent state variable**; its ODE and the elliptic `p_film` solve **must be integrated once per physical step in `src/simulation/m_time_steppers.fpp`, NOT inside `s_compute_rhs`.** `s_compute_rhs` (and `s_compute_mthinc_normals` at `m_rhs.fpp:579`) is called **once per Runge–Kutta substage (3× for RK3)**; integrating a history-dependent ODE there would advance drainage **3× per step with the substage `dt`**, corrupting `h(t)`, the `h_crit` gate, and `p_film∝1/h³`. The split:
- **Once per step (time-stepper):** solve the elliptic `p_film` (Section 3.3) from step-start `film_h` and resolved kinematics; integrate `film_h` (sub-cycled on `dt_film` for stiffness, Section 3.8); compute the rupture min-reduction (Section 3.6); **store `p_film` as a device field.**
- **Per substage (RHS):** the momentum/energy **sources** read the stored `p_film` and are added with `+=` so they get RK-weighted correctly (below). Freezing `p_film` across the substages relaxes strict RK order for the source; this is acceptable for an operator-split source and is the pragmatic default, but flag the order caveat.

On band turn-on initialize `film_h = h_Est`; outside `β`, peg `film_h = h_Est`. Use a one-cell-eroded mask `β_reduced` for the `∇_t` terms (Musehane Eq. 3.63) to avoid a band-edge discontinuity. Local `a = r_h` and `V_n` come from the resolved interface velocities (`mom%beg..mom%end`).

**2d. Lubrication source — RHS hook, and do NOT copy the body-force zero-init (MFC P6).** Add a guarded call in the **post-flux global-source block of `src/simulation/m_rhs.fpp` (828–862)** — the block that already hosts `acoustic_source`, `bubbles_euler`, `bubbles_lagrange`, `chemistry`, `cont_damage`:
```fortran
if (film_drainage) call s_compute_film_drainage_source(q_prim_qp%vf, q_cons_qp%vf, rhs_vf)
```
Implement it **modeled on** `src/simulation/m_body_forces.fpp:79–137` (`s_compute_body_forces_rhs`), the canonical momentum+energy source pattern, **but with two copy-traps called out:**
- **Do not copy the template's zero-init loop (lines 89–99).** The template zeroes `rhs_vf(mom:E)` because body forces are assembled standalone; at the 828–862 hook **`rhs_vf` already holds the flux divergence**, so a literal copy would wipe out all transport. Use `+=` only.
- **Read `q_prim_qp%vf`/`q_cons_qp%vf`** (the hook neighborhood's arrays; pressure at index `E` in `q_prim`, total energy at `E` in `q_cons`), **not** the template's `q_prim_vf`/`q_cons_vf`.

The source, over `do l=0,p / k=0,n / j=0,m` under `GPU_PARALLEL_LOOP`, applies (Sections 3.5, 3.7):
```
rhs_vf(eqn_idx%mom%beg + d)%sf(j,k,l) += p_film · n̂_θ,d           ! load transmission, n̂_θ = ∇α_θ/|∇α_θ|
rhs_vf(eqn_idx%E)%sf(j,k,l)          += (reversible work term)     ! keeps load transmission energy-consistent
                                     +  Φ deposited as band-gas internal energy   ! irreversible Poiseuille loss (Sec 3.5)
```
All needed fields are in scope and device-resident: velocities (`mom`), pressure (`E` in `q_prim`), volume fractions (`adv`), partial densities (`cont`), color (`c`), geometry (`x_cc,dx,…`), and `film_h`/`p_film`. Enforce the magnitude guard `p_film ≲ C·ρU_r²` (Section 3.8) here.

**2e. New case parameters.** `film_drainage` (logical), `h_crit`, `A_hamaker` (reserved, unused in v1), `xi_film` (ξ), `lambda_mfp` (λ = 66 nm), `f_kn_model`, `p_film_cap` (the `C·ρU_r²` guard), `film_subcycle` (sub-cycle count / `dt_film` policy). Read in `m_global_parameters.fpp` and `GPU_UPDATE(device=[…])` near line 887 (alongside `int_comp, ic_eps, ic_beta`).

### Phase 3 — Rupture criterion + merge handoff

**3a. Rupture test.** Once per step (time-stepper, alongside 2c), reduce `min_rim(film_h)` over each film band (masked min-reduction) and set a merge flag when `< h_crit`.

**3b. Merge handoff (Musehane fuse — default).** In flagged cells fuse the markers: `α_new = α₁ + α₂` into marker 1, zero marker 2; recompute partial densities so **mass, momentum, and energy are conserved cell-by-cell** in the fused cells (the delicate step — a topology change in a conservative compressible scheme; see R7). Turn off the injected lubrication source and dissipation in the fused band (single interface now). mTHINC then sharpens the merged interface naturally.

**3c. Alternative (Mason mode).** Instead of fusing, **switch off the lubrication source** in the ruptured band and let mTHINC / the diffuse interface reconnect across the now-thin gap. Cheaper, avoids field surgery, but relies on the auto-merge behaving (and, per MFC P1, the color force *wanting* to merge helps here). Provide both; **default to the fuse** for control, keep Mason mode as a fallback if the fuse proves non-conservative or noisy.

---

## 5. Validation ladder

Run in order; each rung gates the next. Domain `9D×9D×9D` (half-domain for the symmetry-plane path), initial center-to-center `3D`, outflow BCs on the open sides, Courant 0.8. **Head-on rungs use the reflective symmetry plane at `x=0`** (Section 3.9): free-slip BC at the plane, single drop, no markers. **tanh grid clustering at the collision plane** (`src/pre_process/m_grid.f90:51–70` serial / `178–194` parallel): transition band `x_a=−δ, x_b=+δ` around 0, tune `a_x`/`loops_x` to reach `D/200`–`D/500` at the plane, coarsening to `~D/10` in the bulk; requires `x_domain%beg = −x_domain%end`. Note (MFC P7): `dx = minval(...)` (`m_grid.f90:~68`) with `s_mpi_reduce_min` sets the global `dt` from the finest cell — clustering buys *locality*, not a larger `dt` (consistent with R4); the only relief from the finest-cell `dt` is sub-cycling the film sub-system (Section 3.8).

**Rung 0 — coupled-solver unit tests (before any collision).** (i) Flat-plate limit: the elliptic `p_film` solver reproduces Stefan–Reynolds `p_flat(r)` and `F_lub` to tolerance. (ii) Dimple existence: a prescribed squeeze produces a **center-thick, rim-thin** `h` and a **rim pressure feature** — if it does not, S2/Rung 4–5 cannot pass, stop here. (iii) Sign test (Physics P9): on two markers, verify `f_film` is anti-parallel (drops pushed apart). (iv) Positivity test (MFC P4): `α₃_face ≥ α_min` across all reconstructed faces in all three directions.

**Rung 1 — S1: Fig. 6 head-on bounce (We=8.6), symmetry-plane path.** Confirm the drop **decelerates, deforms, and rebounds** with the film gating merge. Measure `τ_b` as the interval during which `p(0,0) > p_ambient` (Chen & Yang's definition). Self-consistent targets (recomputed in review):
```
τ_c = D/U_r ≈ 306µm / 0.99 ≈ 309 µs
Eq.10: τ_b/τ_c = 0.82√8.6 = 2.40   →  τ_b ≈ 745 µs
Eq.11 check: τ_b = 0.82√(ρ_l D³/σ) ≈ 745 µs   (agrees to ~% with Eq.10 — see §1.2 velocity caveat, not exact)
τ_o (Eq.12) = (π/4)√(ρ_l D³/σ) ≈ 713 µs
```
**Also run the R1 experiment here (Section 3.2):** film source OFF → expect premature merge / wrong deceleration; ON → correct `h(t)`. Judge by *merge outcome and deceleration*, not by whether it rebounds (the capillary spring rebounds it either way on this path).

**Rung 2 — near-contact resolution study, reconciled with the film-solver well-posedness tension (MFC P7 / R5).** Rerun S1 at plane resolutions `D/100, D/200, D/500`. **Pass = bounce outcome and `τ_b` converge and are *insensitive* to cell size** above the band-detection resolution — the load-bearing verification that the model *works*: unlike Chen & Yang (where outcome is a deliberate function of min grid size), a correct sub-grid model must be **grid-independent**. **But note the real tension:** at `D/200` the max-deformation film (~0.3 µm) is 0–1 cells thick, and the in-plane Reynolds stencil needs ≥2 cells — "fall back to the ray integral" does not create resolution that isn't there. So this rung has a second job: **determine whether `D/200` even makes the film solver well-posed, or whether `D/500` (4× the steps) is required** purely for stencil support. Decide `D/200` vs `D/500` here; this may change the cost story from `~10⁶` to `~4×10⁶` steps, and that is an honest possibility, not a failure.

**Rung 3 — S2: nonmonotonic head-on transition, symmetry-plane path.** Sweep `We` across the bouncing window **using transition targets read from Chen & Yang's regime map at tetradecane/N₂-1atm** (Physics P8 — do not assume Qian & Law's 1.88/10.3 unless the condition matches). Reproduce the merge→bounce lower boundary and the bounce→merge upper boundary. The lower boundary must emerge from the drainage gate (`t_drain ∝ a⁴`); the **upper boundary must emerge from rim thinning** (`min_rim h(t) < h_crit` while pressed, Section 3.6), **not** from `t_drain ∝ a⁴` and **not** from grid tuning. **Report S2 sensitivity to `f(Kn)` and to `λ`** (Physics P7) as first-class results. Cross-check the `r_h`↑ / `h_min`-nonmonotonic trends (Chen & Yang Fig. 17).

**Rung 4 — diagnostic quartet.** For S1 and the S2 cases extract `h_center`, `h_min`, `r_h`, `p(0,0)` vs `t/(D/U)` and compare *shapes* to Chen & Yang Figs. 13, 15: `h_min` decreasing faster than `h_center`, the **nonmonotonic `h_min`** (rapid rise then sudden drop at rebound), and the **rim pressure peak / negative rim gauge pressure during bounce-off** (Fig. 15b, t5–t7). These shapes are the physics and are the direct test that the coupled solve produces a real dimple; matching them is stronger evidence than matching `τ_b`. (A flat-film closure would fail this rung — which is why v1's closure was replaced.)

**Rung 5 — Fig. 16 pressure collapse (also the runtime stability guard).** Confirm the film-center pressure is roughly constant during a bounce and **collapses when normalized by `ρU_r²`** across We. Per Section 3.8(b) this doubles as a **runtime assertion**: if injected `p_film` does not collapse on `ρU_r²`, the run is unphysical.

**Rung 6 — restitution (Eq. 13) — the dissipation test (Physics P6).** Confirm `τ_b/τ_o = 0.29 We^(−0.35) + 0.88`. This directly checks that the **sub-grid dissipation** (Section 3.5) is right; an artificially elastic bounce (the failure mode of v1's "add `u·f`") shows up as `τ_b` too long / restitution too high here.

**Rung 7 — S3: oblique bounce (B=0.6), full multi-marker.** Full 3-D, **no symmetry plane** (the only stage that needs the multi-marker machinery, Phases 1–3). This is gated on Phase 1's *driven-together* validation passing (the markers stay distinct under compression, which requires the per-marker-σ / color-repaint fix, MFC P1). Confirm off-center bounce with the film **rotating about the collision center** (Chen & Yang Fig. 11) and drops separating before rupture. Compare qualitatively to Fig. 10 (B=0.9).

**Bouncing-time correlations to report** (Eqs. 10–13): Eq. 10 `τ_b/τ_c = 0.82√We`; Eq. 12 `τ_o = (π/4)√(ρ_l D³/σ)`; Eq. 13 `τ_b/τ_o = 0.29 We^(−0.35) + 0.88`.

---

## 6. Risks, open questions, and decision points

**R1 — Is the injected force necessary, and (reframed) for what? (Physics P2)** The force's job is **load transmission / correct deceleration and prevention of premature numerical merge**, *not* rebound (rebound is the resolved capillary spring). **Decision gate (revised metric):** in Phase 0, run S1 with markers/drop + resolved gas only (force OFF) and measure **whether the drops merge prematurely and whether the deceleration and `h(t)` are correct** — *not* whether `τ_b` is under-predicted (that is a capillary test). Cheap, run early. Honest uncertainty: we do not know a priori how much load the `D/200` cushion already transmits.

**R2 — Continuum lubrication breaks down at the rupture scale.** At `h_crit ≈ 30 nm`, `Kn ≈ 2.2` — transition/free-molecular; the closure is a continuum extrapolation there. **Mitigation:** the *load-bearing* pressure is built at the dimple (`h ~ 0.3 µm`, `Kn ≈ 0.22`), where continuum+slip is defensible; `h_crit` only *gates* merging. So the force is on firmer ground than the final-rupture instant. Still, `f(Kn)` must be sourced properly (R-Kn below).

**R3 — `h_crit` calibration.** Sensitivity is mild above ~50 nm and sharpens below (Mason Fig. 6). Calibrate `h_crit` (default 30 nm) so S2 reproduces both transition boundaries; do **not** tune it per-case (that would hardcode an outcome).

**R4 — Timestep cost of near-contact refinement is real.** Refining the plane to `D/200`–`D/500` sets the global acoustic `dt` (Section 2.3): ~10⁶ steps at `D/200`, ~4×10⁶ at `D/500`. **Coupled with R5 and R-stiff below**, the *effective* cost may be higher: `D/500` may be forced for well-posedness (R5) and/or the film source may need sub-cycling (R-stiff). Stop at the coarsest resolution that passes Rung 2.

**R5 — Film-region representation at `D/200` is a likely failure, not a tail risk (MFC P7).** At max deformation the film (~0.3 µm) is below one `D/200` cell (1.5 µm), so `β` may be 0–1 cells thick and the tangential-gradient stencils in the Reynolds solve (need ≥2 in-plane cells) have no room. "Fall back to the ray integral" does not create resolution. The tension between "stop at `D/200` for `dt`" (R4) and "need an in-plane stencil in the film" is **unresolved by design** and may force `D/500`. Rung 2 decides it; watch for band dropout (film not detected).

**R6 — Compressible-framework subtleties.** (a) **Energy partition (Section 3.5):** the reversible work term keeps load transmission consistent, but the **irreversible Poiseuille dissipation must be deposited as gas internal energy**, not returned to the drop — v1's "add `u·f`" alone makes the bounce artificially elastic. Getting the conservative discretization of this split right is delicate; Eq. 13 (Rung 6) arbitrates it. (b) Gas compressibility is negligible (`ΔP ~ 10² Pa`) — the `ρh` form is harmless over-engineering. (c) The single color function is *attractive* in the film (MFC P1); per-marker-σ / color-repaint is a **prerequisite for S3**, avoided for S1/S2 via the symmetry plane.

**R7 — Merge handoff conservation.** Fusing markers (Phase 3b) is a topology change in a conservative compressible scheme; naive fusion can violate mass/energy conservation or inject acoustic transients. Verify conservation cell-by-cell; fall back to Mason-mode (3c) if fragile. Most numerically delicate step.

**R8 — GPU cost, dominated by the smoothing solver (MFC P7).** The single biggest performance/complexity unknown is the **per-marker smoothing**; v2 specifies an **explicit few-iteration Jacobi/diffusion smoother** (Section 2b) rather than a global implicit solve, precisely to keep it in character with the explicit code. Plus `film_h`/`p_film` fields, masked min-reductions, the elliptic `p_film` solve, and the source kernel. Budget is uncertain and likely **>30–50 %**; profile after Phase 2. All new arrays/kernels/scalars need `GPU_DECLARE`/`GPU_PARALLEL_LOOP`/`GPU_UPDATE` or the run silently falls off the device.

**R9 — Does multi-marker with a single color field hold the drops apart? (MFC P1 — now the central S3 risk.)** No, not intrinsically: single-pressure/single-velocity gives zero relative slip; distinctness survives only while a ≥2–3-cell resolved gas gap persists, and the color force actively *pulls the film shut* (`mpp_lim` in the `num_fluids>2` branch, `m_riemann_solver_hllc.fpp:839`, then fuses on clip). **This is why v2 (i) uses the symmetry plane for S1/S2 to avoid the risk entirely, (ii) makes per-marker-σ / color-repaint a prerequisite for S3, and (iii) replaces the weak passive Phase-1 test with a driven-together test.** Second decision gate; surfaced at Phase-1 validation.

**R-Kn — `f(Kn)` directly moves the S2 thresholds (Physics P7).** `f(Kn)=1+6Kn` is valid only for `Kn ≲ 0.1`; at the dimple `Kn≈0.22`, at rupture `Kn≈2.2`. Because `μ_eff=μ_g/f` rescales both the force and the drainage flux, it *sets* the transition `We`. **Obtain the real `Q(Kn)`/`f(Kn)` (Zhang & Law, or 2025 rarefied dual-VOF/augmented-vdW papers) before S2, and report S2 sensitivity to `f(Kn)` and to `λ` (66 nm vs 100 nm).**

**R-stiff — the binding timestep near contact may be source stiffness, not acoustics (Physics P5; MFC P3).** `p_film ∝ 1/h³` over two decades of `h` is a stiff explicit source that can blow up the momentum/energy RHS before the acoustic CFL bites. **Mitigation:** two-way self-limiting `V_n` (Section 3.8a), a `p_film ≲ C·ρU_r²` runtime guard (3.8b), and **sub-cycling / semi-implicit treatment on a `dt_film` clock** (3.8c). The Section-2 `~10⁶`-step budget assumes the film is sub-cycled off the global clock; without that, the global `dt` collapses.

**R-close — N-fluid closure positivity (MFC P4).** Independent per-marker sharpening can drive `α₃_face < 0` in the film, whose clip *is* the numerical merge. **Mitigation:** renormalize liquid markers to `Σα ≤ 1−α_min` before closure and prove `α₃_face ≥ α_min` (Rung 0-iv). This is the crux robustness problem of Phase 1, not a mechanical index edit.

---

## 7. References

**Primary target (physics + validation):**
- Chen, X. & Yang, V. (2020). *Direct numerical simulation of multiscale flow physics of binary droplet collision.* Phys. Fluids **32**, 062103. (Resolved-film AMR-DNS baseline; Figs. 2, 6, 7, 10, 13, 15–17; **regime map — source of S2 transition We at the target condition**; Eqs. 10–13.)

**Sub-grid film models (design backbone):**
- Musehane, N. M. (2016). *Direct Numerical Simulation of Bubble–Bubble and Droplet–Droplet Interaction Using a Surface Thin Film Model.* MSc, Univ. Cape Town. (Multiple-marker + Reynolds film; Eqs. 3.2–3.10, 3.19–3.20, 3.43–3.69. Local: `…/tool-results/out.txt`.)
- Musehane, Oxtoby & Reddy (2018). *Multi-scale simulation of droplet–droplet interaction and coalescence.* J. Comput. Phys.
- Coyajee & Boersma (2009). Multiple-marker method (credited origin).
- Mason, Stevens & Harvie (2012). *Multi-scale volume of fluid modelling of droplet coalescence.* CSIRO 9th Int. Conf. CFD Minerals & Process Ind., paper 179MAS. (Quasi-1D film + two-way body force; Eqs. 5–10. Local: `…/tool-results/mason2012.pdf`.)
- Harvie & Fletcher (2001); Harvie (1999). (Source of Mason's `dh/dt` closure; not retrieved.)

**Film-drainage theory:**
- Zhang & Law (2011). *An analysis of head-on droplet collision with large deformation in gaseous medium.* Phys. Fluids **23**, 042102. (Deformation-energy + rarefied lubrication + vdW rupture; nonmonotonic transition. **`Q(Kn)`/`f(Kn)` paywalled — must be obtained before S2, R-Kn.**)
- Kwakkel et al. (2013). *J. Comput. Phys.* **258**, 96–125. (Multiple-marker CLSVOF + Zhang–Law drainage; `t_c ≥ t_drain` gate.)
- Rajkotwala et al.; Jiang & James (2007); Yoon et al. (vdW / disjoining-pressure approaches and their documented failure modes at low Ca / final drainage).
- Chen & Wu (2022), arXiv:2205.03604 (rarefaction reframing of Zhang & Law; local PDF).
- Deblais et al. (2024), arXiv:2402.00500 (interferometric jump-to-contact ≈30 nm; `h_crit` calibration; local PDF).
- Sprittles, *Gas Microfilms in Droplet Dynamics: When Do Drops Bounce?*, Annu. Rev. Fluid Mech. (review; paywalled).
- Qian & Law (1997) (regime boundaries; **1.88/10.3 are condition-specific — use only at matching gas/pressure, P8**).
- Dual-VOF w/ rarefied gas + augmented vdW (2025), Int. J. Multiphase Flow S0301932225000850; *Acta Mech. Sinica* (2025) 10.1007/s10409-025-25492-x (modern rarefied + augmented-vdW closures; recommended source for the exact `f(Kn)` and disjoining terms; paywalled).

**MFC implementation hooks (file:subroutine — verified against source in review):**
- `src/simulation/m_thinc.fpp` — `s_compute_mthinc_normals` (207–291), `s_thinc_compression` (295–406), **MTHINC branch 324–356 / per-marker loop 327–355** (int_comp=2; *not* the plain-THINC 357–395 branch), storage (36–38, 199–204). *N-fluid extension + closure positivity + reshaped-index test.*
- `src/simulation/m_weno.fpp:1461`, `src/simulation/m_muscl.fpp:226` — reconstruction callers (directional/reshaped-index correctness).
- `src/simulation/m_riemann_solver_hllc.fpp:839` — `mpp_lim` clip in the `num_fluids>2` branch (silent marker fusion mechanism, R9/R-close).
- `src/common/m_derived_types.fpp:137–155` (`eqn_idx_info`); `src/common/m_global_parameters_common.fpp:161–215` (5-eqn layout), **single color `c` at 296**; `src/pre_process/m_assign_variables.fpp:545` (`c` painted from `patch_icpp%cf_val`). *Index ranges + the coalescence-promoting color-force problem.*
- `src/simulation/m_surface_tension.fpp` — `s_get_capillary` (236–327), `s_compute_capillary_source_flux` (74–233), color grads 255–279, `c_divs(num_dims+1)` at 131. *Per-marker-σ prerequisite for S3.*
- `src/simulation/m_time_steppers.fpp` — **once-per-step** film thickness integration + rupture min-reduction (P2; must NOT live in the per-substage RHS).
- `src/simulation/m_rhs.fpp:828–862` — source hook (per-substage); `m_rhs.fpp:579` — per-substage `s_compute_mthinc_normals` call; template `src/simulation/m_body_forces.fpp:79–137` — **do not copy the 89–99 zero-init**; use `q_prim_qp%vf`/`q_cons_qp%vf`.
- `toolchain/mfc/case_validator.py:427` (int_comp), `:772` (surface_tension). *Guard relaxation.*
- `src/pre_process/m_grid.f90:51–70` / `178–194` (tanh clustering), `~68` (`dx = minval`, global-dt-from-finest-cell).
- `src/common/include/parallel_macros.fpp` — `GPU_PARALLEL_LOOP` (25), `GPU_DECLARE` (114), `GPU_UPDATE`, `GPU_ROUTINE` (57), `GPU_LOOP` (132). *Offload decoration for all new kernels/fields/scalars; note the smoother is an explicit few-iteration solve, not a global implicit one (R8).*