# Extending MFC to Reproduce Droplet Bouncing via a Sub-Grid Gas-Film-Drainage Model

**Status:** Design proposal (v3 — film closure replaced after adversarial verification; energy accounting and citations corrected)
**Author:** Lab assistant (MFC extension)
**Date:** 2026-07-05
**Audience:** MFC domain experts committed to the MFC codebase (compressible, structured-Cartesian, no AMR, diffuse-interface)

**Change log v2→v3 (load-bearing revisions).** (v1→v2's changelog is preserved in git history at commit cda59df9.) Throughout, "M1" now names the first *implementation* milestone; v1/v2/v3 refer only to revisions of this document.

1. **Film closure replaced — v2's system was algebraically degenerate (review finding M1, the headline fix).** Substituting v2's elliptic Reynolds equation into its drainage equation cancels the Poiseuille flux identically, leaving `∂(ρh)/∂t = −2ρV_n`: uniform thinning at twice the kinematic rate, no dimple possible, Rung 0(ii) fails by construction. v3 closes the system with **capillary compliance**: one film mass-conservation PDE plus a *local* Young–Laplace pressure closure, yielding the classical fourth-order thin-film equation (Yiantsios & Davis). The global elliptic band solve **disappears** (Sections 3.3, 3.8, Phase 2c, R8). A non-degeneracy regression test is mandatory at Rung 0.
2. **Dimensional and energy fixes:** the Poiseuille dissipation is `Φ = (h³/12μ_eff)|∇_t p_film|²` — the spurious ρ is dropped (C1); the momentum source is the CSF interface-delta form `f_film = p_film ∇α_θ`, not `p_film n̂` (C11); the energy deposit is defined as the power the force pair *extracts* from the resolved drops, so total energy is conserved by construction rather than gaining Φ from nowhere (M5); a new double-counting honesty flag is added (Section 3.5).
3. **Citations corrected:** the `mpp_lim` marker-fusion mechanism is re-pointed from `m_riemann_solver_hllc.fpp:839` (a non-clipping branch our configuration never enters) to the actual clip sites (F9, F19); Kwakkel volume/pages fixed (M8); the S3 target is the paper's actual off-center case B=0.9/We=48.8 (M3); We = 1.88/10.3 are restored as valid S2 anchors — Chen & Yang attribute them to exactly the target condition (F7); Chen & Yang's own λ ≈ 0.1 µm annotation is acknowledged (F5); `h_min ~ 10⁻⁴D` is qualified as a graphical estimate (F14).
4. **Turning-point claim softened everywhere (C6):** the *viscous-drainage* contribution `∝ V_n` vanishes at maximum compression, but the film continues to transmit the quasi-static capillary load (the Fig. 16 `ρU_r²` plateau). It is never the spring — but `p_film → 0` at the turning point is *not* exact and must not be encoded as a check.
5. **Half-domain bookkeeping box added (C10):** the symmetry-plane path hides factor-of-2/4 traps in mobility, rupture gate, and dissipation; Section 3.9 now fixes the full-gap convention.
6. **Phase-1 gate made coherent (M4):** split into a force-OFF characterization run (expected outcome: merge) and a force-ON gate. **Validation tolerances attached (M9)** as explicit provisional bands. Validator-guard status updated to the committed state (M6). Grid-stretch requirement corrected to origin-anchoring (F16). The 8 MPa estimate recomputed with the rarefaction correction the right way up (C9).
7. **Supportive paper evidence added (F6):** Chen & Yang's own off-center case ran at only ~0.7 µm minimum grid — O(10⁻³)D, MFC-affordable — and a published two-VOF (marker-per-drop) anti-merge precedent exists (He, Xia & Zhang); both strengthen the S3 rationale, with the honest caveat that Fig. 10 is not 15-nm film-resolved DNS.
8. Hygiene: journal names spelled out; broken local reference paths removed; read/paywall status completed (M10); references to the first implementation release renamed M1 (4-doc-hygiene).

---

## 1. Goal and scope

### 1.1 Objective

Reproduce **droplet bouncing** — two liquid drops that approach, deform, and rebound **without ever touching**, cushioned by a thin gas film — inside MFC, a compressible structured-Cartesian diffuse-interface solver with **no AMR**. Chen & Yang (2020) reproduced this in Gerris by *resolving* the film to ~15 nm with octree AMR and balanced-force surface tension (head-on cases; their off-center case is coarser — Section 2.3). That route is infeasible in MFC (Section 2). MFC's route is to **model the film sub-grid** — track its thickness and drainage below the smallest cell, gate coalescence on a physical rupture criterion, and keep the two approaching liquid bodies mechanically separated so the *resolved* capillary spring can rebound them.

One correction of emphasis up front, because it reorders everything downstream: **the resolved physics already contains the restoring force.** The rebound is stored as surface-tension energy in the deformed drop (deformation scale `~D`, well resolved at `D/200`) and released elastically. The sub-grid film's jobs are (i) to **transmit the collision load** so the drops decelerate and deform correctly and (ii) to **gate merging** via `h(t)` against a rupture thickness `h_crit`. The film force is dissipative; its viscous-drainage contribution (`∝ V_n`) vanishes at the turning point, though the film continues to transmit the quasi-static capillary load between the deformed drops throughout the pressed phase (Section 3.2) — what it never does is act as the spring. The whole design is organized around this.

### 1.2 Fluid system and dimensionless targets

Throughout, the physical system is **tetradecane droplets in nitrogen at 1 atm** (Chen & Yang's system): density ratio `ρ_l/ρ_g = 666`, viscosity ratio `μ_l/μ_g = 119`. Properties used for sizing: `ρ_l ≈ 763 kg/m³`, `σ ≈ 0.0265 N/m`, `μ_l ≈ 2.14×10⁻³ Pa·s`, `μ_g ≈ 1.8×10⁻⁵ Pa·s`, liquid sound speed `c_l ≈ 1300 m/s`, gas sound speed `c_g ≈ 349 m/s`, **N₂ mean free path `λ ≈ 66 nm` at 1 atm, 293 K**. (Attribution note, F5: Chen & Yang's own Fig. 5 annotates "mean free path of gas molecules: ~0.1 µm". We use the standard kinetic-theory value 66 nm and treat their ~0.1 µm as a rounded order-of-magnitude annotation; because λ scales every `Kn` and therefore every rarefaction correction and merge threshold, **S2 sensitivity to λ is a first-class reported result** — R-Kn.)

Governing groups (Chen & Yang Eqs. 1–4), with `U_r` the relative (collision) velocity and `D` the initial drop diameter:

```
B  = X/D                     (impact parameter, X = centroid offset ⟂ motion)
We = ρ_l U_r² D / σ          (inertia vs surface tension)
Re = ρ_l U_r D / μ_l
Oh = μ_l/√(ρ_l σ D) = √We/Re  (liquid viscous dissipation)
τ_c = D/U_r                  (convective time)
```

> **Note on the velocity label (do not advertise as exact).** Chen & Yang's Fig. 6 lists `U = 0.97 m/s, D = 306 µm, We = 8.6, Re = 105.9`. Back-solving `U_r = √(We·σ/(ρ_l D))` with tetradecane properties gives `U_r ≈ 0.99 m/s`, i.e. the tabulated `0.97 m/s` is the **relative/collision** velocity `U_r`, with each drop at `U_r/2 ≈ 0.49 m/s`. We adopt `U_r` as the characteristic velocity in `We`, `Re`, `τ_c`. Caveat verified in review: `We=8.6` back-solves to `U≈0.99` but `Re=105.9` back-solves to `U≈0.97` — a **~2 % internal inconsistency**, well within property uncertainty but *not* an exact identity; correlation checks in Section 5 inherit this ~few-% floor.

### 1.3 Staged scope (each stage gates the next)

| Stage | Target (Chen & Yang) | We | Re | B | Purpose |
|---|---|---|---|---|---|
| **S1** | **Fig. 6 head-on bounce** | 8.6 | 105.9 | 0 | First light: one clean head-on bounce; confirm the drops decelerate, deform, and rebound with the film gating merge. Runnable on the **symmetry-plane / single-drop path** (Section 3.9) — no two-marker machinery. |
| **S2** | **Nonmonotonic head-on transition** | sweep (see below) | — | 0 | Reproduce merge→bounce (low-We boundary) and bounce→merge (high-We boundary). Tests the drainage-time / rim-`h_min` / `h_crit` gate. Symmetry-plane path. |
| **S3** | **Off-center bounce, Figs. 10–11** | **48.8** | **260.3** | **0.9** | 3-D, no symmetry plane; rotating film. The only stage that *requires* the general (Musehane) multi-marker film model over the axisymmetric (Mason) one. **Primary case: the paper's only off-center bouncing simulation** — `We=48.8, Re=260.3, B=0.9, D=306 µm, U=2.31 m/s` (Fig. 10; Fig. 11's rotating film is the *same* case). `B=0.6` is retained only as a **secondary regime-map (bounce/no-bounce) check** against Fig. 2 (M3). |

> **S2 transition targets (F7 — v2 over-corrected).** Chen & Yang (p. 062103-12/13) attribute **`We = 1.88` (merge→bounce) and `We = 10.3` (bounce→merge) specifically to head-on tetradecane in nitrogen at 1 atm** — exactly the target condition — citing the experiments of Jiang, Umemura & Law (1992), *Journal of Fluid Mechanics* **234**, 171 (their ref 7 — *not* Qian & Law). **They are therefore valid S2 anchors here.** The general condition-specificity warning stands: Qian & Law's central finding is that the bouncing window is a strong function of ambient gas pressure/composition, so these numbers travel only with the condition. **Circularity note (M10, softened):** Chen & Yang's Fig. 2 solid regime boundaries are reproduced from Qian & Law, while the `1.88/10.3` anchors come from a *different* experiment (Jiang et al. 1992) — so the experimental anchor is independent of the Fig. 2 curves. Chen & Yang still *numerically disabled merging* for their We = 0.2–1.0 runs, so the only independent numerical check in Fig. 2 is their discrete simulation symbols, not the boundary curves; the experimental values remain the primary anchor.

**Diagnostic quartet to reproduce** (Chen & Yang Figs. 13, 15): `h_center` (axis, max film thickness), **`h_min`** (rim of dimple, min film thickness — the primary gate variable), `r_h` (dimple radius), and `p(0,0)` (gauge pressure at collision center), all vs `t/(D/U)`. The **dimpled** film (thicker at center, thinner at the rim), the **nonmonotonic `h_min`**, and the **rim pressure peak / negative rim gauge pressure during bounce-off** are the load-bearing physics. A flat-film closure cannot represent them — and neither could v2's degenerate coupled solve (review finding M1); the v3 capillary-compliance closure (Section 3.3) exists to produce them, and Rung 0(ii) proves it does before any collision run.

**Quantitative pass criteria** (explicit acceptance bands in Section 5, M9): bouncing-time correlation **Eq. 10, `τ_b/τ_c = 0.82√We` (fitted for 0.2 < We < 8.6)**; the **Fig. 16 pressure collapse** (film-center pressure ∝ `ρU_r²`); and the **restitution correlation Eq. 13** (`τ_b/τ_o = 0.29 We^(−0.35) + 0.88`), the direct test that the sub-grid **dissipation** is right (Section 3.5).

### 1.4 Out of scope (M1 — the first implementation milestone)

Shattering; unequal drop sizes; temperature-dependent properties; explicit van der Waals *body* force (we use a fixed `h_crit` cutoff instead, Section 3.6). Reflexive/stretching separation (regions IV/V) are natural follow-ons once merging works but are not M1 targets.

---

## 2. Why model the film instead of resolving it

### 2.1 The two levers Chen & Yang used, and that MFC lacks

Chen & Yang reached ~15 nm (`≈10⁻⁵ D`, level 15) for their head-on cases using **two levers simultaneously**:

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

There is a deeper mismatch: a compressible code pays the acoustic `dt` to resolve *fast liquid acoustics* that are **irrelevant** to a quasi-static lubrication film draining at `O(1 m/s)`. The film physics is essentially incompressible-lubrication; the solver is compressible-acoustic. Resolving the film forces the solver to spend its entire budget on the wrong physics. (Corollary, now load-bearing in Section 3.3: because the physical film pressure jump is only `ΔP ~ 10² Pa`, **gas compressibility in the film is negligible** — the v3 film model therefore treats the film as incompressible; v2's compressible-conservative `ρh` form bought nothing and is dropped.)

### 2.3 What the sub-grid approach buys

Stop refinement at the collision plane at **`D/200`–`D/500` (≈0.6–1.5 µm)** — enough to capture drop deformation and the resolved gas approach — and **model everything below that**. At `dx = D/200 = 1.53 µm`:

```
dt ≈ 0.8·1.53×10⁻⁶/1300 ≈ 9.4×10⁻¹⁰ s
N_steps ≈ 1.5×10⁻³ / 9.4×10⁻¹⁰ ≈ 1.6×10⁶ steps   (~10⁶, verified)
```

**Two decades fewer steps than resolving** (10⁶ vs 10⁸). **Caveat (MFC P3):** this counts *acoustic* steps only. The injected film source is stiff (Section 3.8); the binding constraint on `dt` may not be acoustic at all in the last decade of `h`. The `~10⁶` figure holds only if the film sub-system is integrated on its own (sub-cycled) clock rather than the global one.

The film **must** be sub-grid regardless of how far we refine: at maximum deformation Chen & Yang's film is `z/D ~ 10⁻³` (~0.3 µm at D=306 µm), already **below** `D/500`; at the rim just before rebound `h_min` approaches `O(10⁻⁴)D ≈ 30 nm` — a **graphical estimate** from Fig. 15(a) at We = 1 (the final pre-bounce plunge), not a stated paper number, but consistent with their van der Waals effective distance of 0.03 µm (F14). So across the entire load-bearing phase the film is thinner than the smallest affordable cell. Refinement to `D/200` locates the band and captures deformation; the dimple, rim, and drainage below ~0.6 µm are modeled. **This creates a genuine tension (Section 6, R5): at `D/200` the film is 0–1 cells thick, but the in-plane film stencil needs multiple cells — which may force `D/500` purely to make the film solver well-posed, at 4× the step count.** This is reconciled in Section 5 (Rung 2), not hand-waved.

**Supporting evidence that this resolution class suffices for bouncing (F6).** Chen & Yang's own off-center case (Fig. 10, `We=48.8, B=0.9`) was run with topology-oriented AMR at a minimum grid of only **~0.7 µm — `O(10⁻³)D`**, i.e. inside our `D/200`–`D/500` target band — and still bounced ("positive pressure is established within the gas film and prevents droplet merging"). Two caveats keep this honest: (i) in their method **the merging distance is set by the minimum grid size**, so at 0.7 µm the non-merge outcome is guaranteed by the grid, not by a resolved 15-nm film — **Fig. 10 is not film-resolved DNS**; (ii) the published anti-merge precedent for off-center bouncing at coarse resolution is the **two-VOF-function (one marker per droplet) approach of He, Xia & Zhang** (Chen & Yang's ref 38), which Chen & Yang note "always prompt[s] droplet bouncing" precisely because it has **no rupture criterion** — exactly the gap this design's multi-marker + `h_crit` gate combination closes. Both facts strengthen the S3 rationale (Section 3.1, Rung 7).

---

## 3. Physical model

### 3.1 Model selection: Musehane-general (multiple-marker) as the backbone — but only where needed

Two documented sub-grid architectures exist:

| | **Musehane / multiple-marker** (thesis 2016; journal paper 2018) | **Mason / quasi-1D VOF** (CSIRO 2012) |
|---|---|---|
| Interfaces | N separate markers `α_θ`, each advected + smoothed; per-marker surface tension mechanically forbids merging | Single liquid marker; drops *can* merge numerically |
| Geometry | **General** — not restricted to head-on/equal drops; film detected by marker overlap | **Axisymmetric**, film pinned to an imposed `z=0` symmetry plane |
| Film pressure | Reynolds equation inverted on the band | Solved from reduced NS on the plane |
| Force coupling | One-way / diagnostic (his solver resolves the cushion) | Two-way body force |
| Rupture | `h < h_crit` fuses markers | `h < h_crit` (40 nm) switches a collision-plane ghost BC |
| Validation | Min film thickness for bouncing **matched thickness-based-AMR DNS** | Reproduced Pan et al. 2008 head-on tetradecane |

**Decision:** build the **general 3-D capability on Musehane's multiple-marker mechanics** (per-marker `α_θ`, geometric film detection, `h_crit` merge) **for S3**, because Mason's symmetry-plane assumption cannot represent the *rotating* off-center film (Chen & Yang Fig. 11) that defines region-II bouncing. Committing to Mason would forfeit S3. The two-VOF precedent of He, Xia & Zhang (Section 2.3) shows marker-per-drop anti-merge works for off-center bouncing at MFC-affordable resolution; our `h_crit` gate supplies the rupture physics that precedent lacks. **Honesty note:** the resulting architecture — a dimpled sub-grid film with Mason-style two-way force injection — is a hybrid neither backbone reference validated end-to-end (Musehane's coupling was one-way/diagnostic; Mason's two-way force was flat/quasi-1-D). The validation ladder, not provenance, carries the burden.

**The review's key structural insight (Physics "opportunity"; MFC P1):** for **head-on S1/S2 you do not need the multi-marker machinery at all.** A reflective symmetry plane at the collision midplane physically holds the drops apart, enforces `V_n` symmetry for free, halves the domain, and — critically — **side-steps the entire "does the single color field hold two markers apart" risk** (R9). It also cleanly isolates whether the bounce is real physics (resolved capillary spring + film load transmission) versus an artifact of the marker trick. **The work is therefore staged so the film physics is validated first on the symmetry-plane single-drop path (Section 3.9 / Phase 0), and the multi-marker machinery is only built for S3.**

The multi-marker port still maps onto MFC's machinery — mTHINC is a per-marker extension away from N-fluid, and two liquid markers of the *same* material are exactly what MFC's `cont`/`adv` index ranges support — but with two hard caveats that Section 4.1 treats as the crux (MFC P1, P4).

### 3.2 What the injected film force actually does (Physics P1, P2; softened per C6)

Framing the injected lubrication force as "the rebound force MFC must add to supply the repulsion" is conceptually wrong.

- The **viscous-drainage contribution** to the film pressure scales with the normal approach velocity, so it vanishes at maximum compression (`V_n = 0`). But the film does not go silent there: it **continues to transmit the quasi-static capillary load** between the two deformed drops — Chen & Yang's Fig. 16 shows a sustained `O(ρU_r²)` pressure plateau *through* the turning point, and their peak `p(0,0)` occurs *before* maximum deformation. Do **not** encode `p_film → 0` at maximum compression as an exact statement or a runtime check (C6). What the film never does is act as the spring.
- During separation the drainage contribution reverses sign — the physical negative **rim** gauge pressure of Fig. 15(b) (a local rim feature during t5–t8, reaching the center only at t9, not a global suction): the viscous film *opposes* rebound.
- The **actual restoring force is capillary**: surface-tension energy stored in the deformed drop, at deformation scale `~D`, which MFC **resolves** at `D/200`. Bouncing can therefore emerge in MFC — but because the capillary spring is resolved, *not* because we inject a repulsion.

So the injected force's real roles are: **(a) transmit the collision load** so the resolved drops decelerate and deform on the correct timescale (without it, an under-resolved `D/200` gas cushion passes too little load and the drops interpenetrate/merge numerically before deforming), and **(b) keep the resolved gas gap open** so the markers/interfaces do not numerically fuse before the capillary spring acts. It is a **load-transmission and anti-premature-merge** term, not a rebound term.

**Consequence for R1's decision gate:** the Phase-0 "force off" experiment is judged by *"do the drops merge prematurely / is the deceleration and `h(t)` history correct?"* — **not** by *"is `τ_b` under-predicted / is rebound adequate?"* The rebound test is a *capillary* test and will (if the drop is well resolved) succeed with or without the force; the force's necessity shows up as **premature merge and wrong deceleration**, which is what to measure.

### 3.3 The film model: one mass-conservation PDE closed by capillary compliance (replaces v2's degenerate system — review finding M1)

**What was wrong in v2 (the load-bearing v3 change).** v2 postulated two equations: an elliptic Reynolds equation that determined `p_film` from the imposed squeeze rate, and a drainage equation for `h` that re-added the same Poiseuille flux divergence. The system is **algebraically degenerate**: substituting the elliptic constraint into the drainage equation cancels the flux identically, leaving `∂(ρh)/∂t = −2ρV_n` — uniform thinning at *twice* the kinematic rate, with **no mechanism for the pressure to reshape the film**. No dimple can ever form, the rim-minimum gate has nothing to act on, and the "dimple emerges naturally" claim was false by construction. The missing physics is **interface compliance**: the pressure must be able to feed back on the thickness distribution.

**v3 closure.** The film is treated as **incompressible** (`ΔP ~ 10² Pa` across the film — Section 2.2 — so the gas density varies by ~0.1%; v2's `ρh` conservative form bought nothing and is dropped). `h` is the **full gap** between the two interfaces. Two relations:

1. **Film mass conservation — the only PDE:**
```
∂h/∂t + ∇_t·(h Q) = ∇_t·[ (h³/(12 μ_eff)) ∇_t p_film ]
    Q = (u₁ + u₂)/2          (mean interface velocity, from the resolved field)
    μ_eff = μ_g / f(Kn)       (rarefaction correction, Section 3.4)
    ∇_t(·) ≡ ∇(·) − (∇(·)·n) n   (in-plane gradient), n = unit interface normal
```

2. **Pressure closure — capillary compliance, not mass balance.** `p_film` is *not* determined by mass conservation; it is **slaved to the film shape** by the linearized Young–Laplace normal-stress balance of the two identical deformable interfaces. Each interface deforms from its resolved-kinematics position by `η = (h − h_rigid)/2`, so
```
p_film(x) = p_edge + σ κ₀ − (σ/2) ∇²_t (h − h_rigid)
    ∂h_rigid/∂t = −V_n(x)      (rigid-kinematics gap: the closure rate the resolved interfaces impose)
    κ₀ = resolved interface curvature;  p_edge = resolved gas pressure at the band edge (σ-jump-subtracted, smoothed)
```
The compliance coefficient `σ/2` (each interface's normal-stress balance gives `p_film − p_drop = −σ∇²_t η` with `η = (h − h_rigid)/2`, hence `σ/2` on `h − h_rigid`) is stated as **provisional** pending the Rung-0 benchmark against Yiantsios–Davis solutions.

**Structural win — the global elliptic band solve disappears.** In v2, `p_film` required a global elliptic solve over the band every step. In v3, `p_film` is a **local algebraic operator on `h`** (a Laplacian of `h − h_rigid`): a stencil evaluation, not a solve. Cheaper, simpler, GPU-trivial. Section 3.8, Phase 2c, and R8 are updated accordingly.

**Why the dimple now emerges.** Eliminating `p_film` gives the classical **fourth-order thin-film PDE**:
```
∂h/∂t + ∇_t·(h Q) = −∇_t·[ (h³/(12 μ_eff)) (σ/2) ∇_t(∇²_t (h − h_rigid)) ] + (resolved p_edge, κ₀ gradient terms)
```
Uniform squeezing + capillary compliance is the **textbook drainage mechanism** (Yiantsios & Davis, Journal of Fluid Mechanics **217**, 547 (1990); and their 1991 drop–drop companion): gas squeezed from the center travels farther to escape, pressure builds toward the center, the compliant interfaces bulge there, and `h` develops the center-thick / rim-thin dimple with `h_min` at the rim — exactly the Chen & Yang Fig. 13 structure.

**Non-degeneracy requirement (mandatory Rung-0 regression test).** A uniform-`V_n` squeeze of an initially flat film **must** develop a nonuniform (dimpled) `h` under this closure. v2's system could not, by construction. No collision run before this test passes (Rung 0-ii).

**Boundary conditions** (the fourth-order PDE needs two at the band edge): `h` matches the resolved gap (`h → h_Est`); `p_film` matches the resolved gas pressure, σ-jump-subtracted and smoothed (Musehane Eq. 3.65 — the surviving, correctly-scoped remnant of the earlier R6b, kept to avoid feeding resolved acoustic noise into the film).

The **same `p_film`** still does three jobs, eliminating underdetermination (Physics P3): (1) drives the drainage flux above; (2) is injected as the momentum source (Section 3.7); (3) is the diagnostic pressure compared to Chen & Yang's `p(0,0)` and Figs. 15/16.

**Head-on axisymmetric reduction (S1/S2 sanity check), corrected to the v3 system.** On the symmetry-plane path the solve collapses to the 1-D radial form (fourth-order in `h` once `p` is eliminated) — a cheap, well-posed check of the general solver:
```
∂h/∂t + (1/r) ∂/∂r ( r h Q_r ) = (1/r) ∂/∂r [ r (h³/(12 μ_eff)) ∂p_film/∂r ]
p_film(r) = p_edge + σ κ₀ − (σ/2) (1/r) ∂/∂r ( r ∂(h − h_rigid)/∂r )
```

**Flat-plate Stefan–Reynolds — retained only as the rigid (zero-compliance) limit unit test.** For a *flat* parallel film of lamella radius `a` and gap `h` squeezed at `V_n` with immobile, rigid interfaces:
```
p_flat(r) = 3 μ_eff V_n (a² − r²)/h³ ,   F_lub = 3π μ_eff V_n a⁴/(2 h³) ,   t_drain ∝ μ_eff a⁴/(F h_crit²)
```
Used **only** to (i) unit-test the film solver in the rigid limit (Rung 0-i) and (ii) reason about scalings (`F ∝ 1/h³`, `t_drain ∝ a⁴`). It is **never injected as the force** — it cannot produce a dimple, a rim minimum, or a rim pressure peak.

**Fallback — kinematically slaved `h` (documented decision point, NOT the default).** Drop the compliance entirely: `h ≡ h_rigid` follows resolved kinematics (`∂h/∂t = −V_n`) and `p_film` comes from the rigid-limit elliptic problem. Robust and cheap, but it has **no sub-grid dimple**, so the rupture gate degrades to the kinematic minimum rather than a rim minimum — weakening exactly the high-We bounce→merge mechanism (Section 3.6). Fall back only if the compliance closure proves stiff or uncalibratable, and say so explicitly in reporting.

Van der Waals disjoining pressure is **omitted from the drainage driver** (following Musehane/Malavé et al.: London–vdW forces do not aid dimple formation); vdW enters only through the fixed `h_crit` rupture cutoff (Section 3.6). This is a standard, defensible modeling choice (Mason, Kwakkel); the caveat is that the *merge outcome* is sensitive to `h_crit` below ~50 nm (R3).

### 3.4 Rarefaction correction (μ_eff) — not a harmless placeholder (Physics P7)

At `h → tens of nm`, `Kn = λ/h → O(1)` and continuum lubrication over-predicts resistance; the correction enters through `μ_eff = μ_g/f(Kn)` in the drainage mobility, and therefore in the film pressure and injected force. With `λ ≈ 66 nm` (Section 1.2):
- at the dimple, `h ~ 0.3 µm`: `Kn ≈ 0.22`, first-order slip `f = 1+6Kn ≈ 2.3` (already marginal);
- at rupture, `h ≈ h_crit = 30 nm`: `Kn ≈ 2.2`, `f ≈ 14` — deep in transition/free-molecular, where first-order slip is **meaningless**.

Because `f(Kn)` rescales `μ_g` in *both* the force and the drainage flux, **it directly sets how fast `h` reaches `h_crit`, i.e. it directly sets the S2 transition `We`.** So `f(Kn)=1+6Kn` is **not** a benign placeholder — S2 outcomes are a strong function of this unknown function. **Requirements:**
1. Obtain the actual `Q(Kn)`/`f(Kn)` (Zhang & Law 2011, BGK-derived; or the 2025 rarefied dual-VOF / augmented-vdW papers) **before running S2**.
2. **Report S2 sensitivity to `f(Kn)` and to `λ`** (66 nm vs the paper's ~0.1 µm annotation) as first-class results, not an afterthought.
3. Use it only for gating; the *load-bearing* pressure is built at the dimple (`Kn ≈ 0.22`), where continuum-plus-slip is defensible, so the force is on firmer ground than the final-rupture instant (R2's mitigation, still valid).

### 3.5 Sub-grid dissipation and conservative energy accounting (Physics P6; C1, M5 fixed)

The physical lubrication resistance **dissipates energy inside the sub-grid film**; that dissipation is why bouncing is lossy and why `τ_b/τ_o` deviates from unity (Eq. 13). Two failure modes bracket the design: injecting only reversible work makes the bounce artificially elastic (restitution too high); injecting dissipation *on top of* the work terms with no matching sink creates energy from nowhere (M5). Both are excluded as follows.

**Poiseuille dissipation, per unit film area (C1 — corrected):**
```
Φ = (h³ / (12 μ_eff)) |∇_t p_film|²          [W/m²]
```
equivalently `Φ = q·(−∇_t p_film)` with `q` the volume flux per unit width. (v2 carried a spurious ρ — dimensionally wrong by exactly a factor of density; inherited from the physics review's P6 and caught in verification. In SI the error was only ~15% for N₂, silently skewing Rung 6; in nondimensional code units it could be orders of magnitude.)

**Energy closure — conservative by construction (M5).** The injected momentum source (Section 3.7) forms an **anti-parallel force pair** across the gap; through its consistent work terms it **extracts mechanical power from the resolved drops** at a rate `Ẇ_ext ≈ ∫_band p_film V_n dA`. **Deposit exactly that extracted power as internal energy in the band gas cells.** Total energy is then conserved **by construction** — the deposit is *defined* as whatever the momentum source removed, not an independent source. In quasi-steady lubrication the extracted power equals `∫Φ dA` (the film neither stores nor creates energy). During suction (`V_n < 0`) the sign reverses (a small effect). **Do NOT add Φ as a free-standing energy source on top of the work terms** — that is net energy creation and biases the Rung-6 restitution. The Eq. 13 restitution correlation arbitrates the physical split.

**Honesty flag — double counting (new in v3).** The resolved band gas is still present and already exerts a (smeared, under-resolved) pressure on the interfaces. The injected source should therefore be the **excess** `(p_film − p_gas,resolved,band)` — or, if injected in full, the overlap with the resolved cushion must be listed as an **open calibration point**. It is so listed; Rung 5's `ρU_r²` collapse is the observable that would expose over-forcing.

**Honesty flag (R6/R-diss, unchanged):** the exact *conservative* discretization of "remove mechanical energy from the resolved field, deposit as gas internal energy" in MFC's single-velocity mixture (where liquid and gas share a cell velocity, so the partition is KE ↔ mixture-internal-energy) is genuinely subtle and is an explicit design risk. M1 target: enforce global energy conservation to machine precision and let Eq. 13 arbitrate the physical partition. Getting this term wrong shows up as the wrong restitution, so it is observable.

### 3.6 Rupture criterion and the coalescence gate — track the *rim* (Physics P4)

**Critical thickness `h_crit`** (calibratable, default **30 nm**): Chen & Yang's vdW effective distance is 30 nm; Mason uses 40 nm; Deblais et al. (2024) interferometry measures jump-to-contact at ≈30 nm; Musehane reaches `O(10⁻⁷ m)`.

**Gate on the rim minimum, on the dimpled field:** merge when `min_rim h(t) < h_crit` (the minimum over the film band, which physically sits at the *rim* of the dimple), evaluated on the compliance-closure dimpled `h` — **not** on a single scalar `h`. (Half-domain convention: gate on the full gap `H` — Section 3.9.) When the gate trips, fuse the markers (Section 4, Phase 3); otherwise they stay distinct and the drops **bounce**.

This is the threshold form of **Kwakkel's `t_c ≥ t_drain` gate** (if `h` reaches `h_crit` while the drops are still pressed, they merge; if they separate before, they bounce); we implement the threshold form because `h(x,t)` is tracked directly.

**The We-nonmonotonicity mechanism (Physics P4).** Claiming *both* S2 transitions "emerge from `t_drain ∝ a⁴`" is wrong for the high-We branch:
- **Low-We merge→bounce:** correct as stated — higher We → larger lamella `a` → `t_drain ∝ a⁴` grows → the film survives → bounce.
- **High-We bounce→merge:** does **not** follow from `t_drain ∝ a⁴` (that scaling predicts drops stay *bounced* as We rises — the opposite). The real mechanism is that increased inertia drives the **rim `h_min` below `h_crit` while the drops are still pressed**: the rim thins faster than the wider lamella can compensate. This is a competition of **`h_min(t)` (at the rim) against `h_crit`**, i.e. precisely the dimple/rim structure that a flat-film closure — or v2's degenerate system — cannot represent. State this explicitly in reporting S2, and **confirm the solver actually produces a rim minimum (Rung 0-ii) before claiming S2**.

### 3.7 The injected force: CSF interface-delta form and sign convention (Physics P9; C11 fixed)

The momentum source is the **volumetric CSF interface-delta form**:
```
f_film = p_film ∇α_θ          [N/m³]
```
`∇α_θ` carries both the direction (`n̂_θ = ∇α_θ/|∇α_θ|`, pointing into liquid θ) and the interface delta `|∇α_θ|` (units 1/m) that converts a surface traction into a volumetric source. **`p_film n̂_θ` alone is Pa — dimensionally wrong as a volumetric RHS** and grid-dependent (a fixed Pa added pointwise over an O(4-cell) band does not converge as a surface force; it would break Rung 2's grid-independence and the magnitude-guard calibration — v2's error, C11). This is the corrected transcription of Mason's two-way force. **Unit test:** the band-integrated force must equal `∫ p_film dA` and converge under grid refinement (Rung 0-iv).

**Sign:** `∇α_θ` points into liquid θ, so `p_film > 0` pushes each interface *into* its own drop, away from the film — repulsive. With a gas-film-outward gradient the sign flips and the force would collapse the drops. Because the two markers have **opposite gradients across the gap**, the two injected forces must be verified **anti-parallel** (drops pushed apart) in a unit test before Phase 2. On the symmetry-plane path there is one interface and one gradient, which removes the ambiguity — another reason to validate there first.

### 3.8 Magnitude and stiffness — self-limiting and a hard guard (Physics P5; MFC P3; C9 fixed)

A naive injection over-predicts by orders of magnitude. Plugging nominal numbers into the flat-film estimate at the dimple (`V_n = 0.5 m/s` — a nominal already-decelerated closure rate; under the full-gap convention of Section 3.3 the *initial* closure rate is `U_r ≈ 1 m/s` — `a ≈ 90 µm`, `h = 0.3 µm`): with the **uncorrected** `μ_eff = μ_g = 1.8×10⁻⁵ Pa·s` it gives **~8 MPa**; with the rarefaction correction `f(Kn) ≈ 2.3` at `Kn = 0.22`, `μ_eff = μ_g/2.3` gives **~3.5 MPa**. (The correction *reduces* `μ_eff` — never write `μ_eff ≈ 2×μ_g`; the rarefaction correction goes the other way, C9.) Either way ~3–4 decades above the physical film pressure, which Fig. 16 says is `O(ρU_r²)` — `ρU_r² ≈ 748 Pa` at We = 8.6, with the Fig. 16 plateau sitting at `p/(ρU_r²) ≈ 5`, i.e. a few kPa. The mismatch is diagnostic: the injected closure scales **viscously** (`μ_eff/h³`) while the physical pressure is **inertially capped at `~ρU²`** by drop deformation. The resolution is that by the time `h ~ 0.3 µm` **the drop has already decelerated, so the instantaneous resolved `V_n ≪ U_r/2`** — the pressure never reaches MPa physically. **Mitigations (all M1):**

- **(a) Two-way self-limiting.** Use the **resolved instantaneous `V_n`** (from the resolved interface velocities), so that injecting `f_film` decelerates the drop → `V_n → 0` → the drainage pressure relaxes. The feedback caps the pressure automatically at the inertial scale — *provided the coupling is tight enough within a step*.
- **(b) Runtime magnitude guard.** Enforce `peak p_film ≲ C·ρU_r²` (`C = O(1–10)`) as a **stability guard**, not merely a validation output. This *is* Rung 5 promoted to a runtime assertion: if the injected pressure does not collapse on `ρU_r²`, the run is unphysical and should trip, not silently blow up.
- **(c) Stiffness control — two mechanisms, and the fourth-order term makes it worse.** The mobility's `1/h³` pressure sensitivity spanning two decades of `h` (1.5 µm → 30 nm) was already a stiff explicit source. The v3 closure adds the **fourth-order capillary term, whose explicit stability limit scales as `dt_film ∝ dx_film⁴`** — brutally restrictive on a fine band grid — so **sub-cycling the film sub-system on its own `dt_film` clock and/or semi-implicit treatment of the capillary term is even more necessary than in v2.** (The global elliptic solve is gone — Section 3.3 — so the film integrator *is* the film sub-system's numerical cost.) A **source-based `dt` limiter** is added alongside the acoustic CFL. **Near contact the binding constraint may be source stiffness, not acoustics** — the Section 2 step-count budget assumes the film is sub-cycled off the global clock; if it is not, the global `dt` collapses and the `~10⁶`-step story fails.

### 3.9 The symmetry-plane path for head-on (S1/S2) — recommended, and it de-risks the marker machinery

For head-on collisions the collision midplane is a plane of symmetry: zero normal velocity, free tangential velocity — exactly a **reflective (free-slip) symmetry BC**. Simulate **one drop in the half-domain**; the "film" is the symmetric squeeze film between the drop and its mirror image, and the free-slip plane correctly reproduces the *mobile-interface* film of two identical drops. Benefits:
- **Halves the cost** and enforces `V_n` symmetry for free.
- **Removes the two-marker distinctness problem entirely** for S1/S2 — one drop, one marker, one color normal, so R9 and the coalescence-promoting color force (Section 4.1, MFC P1) do not arise on this path.
- **Isolates the physics from the numerics:** it directly tests whether the bounce comes from the resolved capillary spring + film load transmission, independent of the marker trick.

> **Half-domain bookkeeping box (C10 — factor-of-2/4 traps; read before implementing).** The *resolved-side* equivalence (free-slip mirror) is **exact for symmetric data** — no factor lurks there. The *sub-grid film* bookkeeping is where the traps hide, because the standard mobility `h³/12μ` assumes no-slip at **both** walls of a gap `h`, while the half-domain's natural variable is the drop-to-plane distance (a half gap, with a zero-shear midplane). Convention, fixed here once:
> - **Solve the film in full-gap variables:** `H = 2·h_plane`, with the standard two-wall mobility `H³/(12 μ_eff)` and rigid-kinematics closure `∂H_rigid/∂t = −V_n` (`V_n` = full relative approach rate = 2× the drop's speed toward the plane). Feeding the half gap into the unmodified solver under-predicts drainage mobility 4× and over-predicts film pressure ~4×, directly moving the S2 thresholds.
> - **Gate rupture on the full gap:** `H < h_crit` (`h_crit = 30 nm` is a full drop-to-drop vdW threshold, not a drop-to-plane one).
> - **Deposit `Φ/2`** (the half-domain's share of the film dissipation) in the half-domain band cells, or the mirrored total double-counts the heating and Rung 6 comes out too lossy.
> - **Unit test:** the half-domain film solver must match the full two-sided flat-plate solution (Rung 0-vi).

**Staging consequence:** Phase 0 (Section 4) validates the *film physics* on this single-drop path; the *multi-marker machinery* (Phases 1–3) is built and validated only for S3, where there is no symmetry plane.

---

## 4. MFC implementation plan (phased, with concrete hooks)

All file:line references are to the current MFC tree and were **verified against source in two review rounds** (corrections from the second round are folded in below). Every new device array needs `GPU_DECLARE`; every new kernel loop needs `GPU_PARALLEL_LOOP`/`END_GPU_PARALLEL_LOOP`; every new case-file scalar needs `GPU_UPDATE(device=...)` (macros in `src/common/include/parallel_macros.fpp`).

### Phase 0 — Symmetry-plane head-on film physics (de-risk first; no markers)

Single drop in a half-domain with a reflective symmetry BC at the collision plane, `num_fluids = 2` (tetradecane drop + N₂), plus the film sub-system (Sections 3.3–3.8, with the Section 3.9 full-gap convention). This runs **S1 and S2** and validates the capillary-compliance film solve, the dimple/rim diagnostics, the dissipation/restitution, and the stiffness handling **without any of the multi-marker risk**. Only after Phase 0 passes do the multi-marker phases below (needed solely for S3) become worth their cost. Phase 0 also *is* the R1 experiment: run S1 with the film source **off** and confirm the drop merges into the plane prematurely / decelerates wrongly, then on and confirm correct `h(t)` and gating.

### Phase 1 — N-fluid interface sharpening (for S3): the two hard problems first

**Setup:** `num_fluids = 3` in a 5-equation model: `fluid 1 = tetradecane (drop-A marker)`, `fluid 2 = tetradecane (drop-B marker, identical EOS)`, `fluid 3 = nitrogen (continuous phase, residual)`.

**Reality check the crux — two identical-EOS markers are *not* intrinsically distinct (MFC P1).** In `model_eqns=2` (5-equation) MFC is **single-pressure, single-velocity**. Two identical-EOS liquid markers share the cell's pressure and velocity (zero relative slip) and are distinguished *only* by being separate advected `α` fields *while they occupy different cells*. Once the resolved gas layer `α₃` thins below the numerical diffusion width (~2–3 cells ≈ `D/100`–`D/70` at these resolutions), markers 1 and 2 co-locate and **nothing in the scheme un-mixes them**: mTHINC re-sharpens each marker toward its own `α=1` independently, their sum exceeds 1, and MFC's positivity/limiting machinery silently clips and renormalizes **right where the film matters most** (F9, corrected citations): the `mpp_lim` face-state clip+renormalize in the plain 5-equation Riemann path (`m_riemann_solver_hllc.fpp:1238–1258` — active at **any** `num_fluids`, no `num_fluids>2` qualifier), the volume-fraction correction that mutates `q_cons_vf` itself (`m_pressure_relaxation.fpp:111–138`, `s_correct_volume_fractions`), and the equivalent clamp in `m_variables_conversion.fpp:285–293`. (Nuance, F19: the α renormalization leaves partial densities untouched — mass-conserving — while zeroing a negative `αρ` does alter mass, but only in already-unphysical states.) **Distinctness is therefore contingent on the injected force holding a ≥2-cell resolved gas gap open at every instant** — it is not a property of the representation. This is the central S3 risk; state it as such.

**Worse — the single color-function force is coalescence-promoting in the film (MFC P1).** `eqn_idx%c` is a **single user-painted advected scalar** (`src/common/m_global_parameters_common.fpp:296`, filled from `patch_icpp%cf_val` at `src/pre_process/m_assign_variables.fpp:545`), and surface tension acts on `∇c` (`src/simulation/m_surface_tension.fpp:255–279`). With the natural painting (liquid `c=1`, gas `c=0`) the film is a `1→0→1` trough in `c`, and MFC's CSF capillary force acts to **minimize color-interface area — i.e. it pulls the two drop surfaces together to erase the gas trough.** So the resolved surface-tension model is not merely "unable to express opposing curvature"; **in the film it is an *attractive*, coalescence-promoting force the injected lubrication source must overcome on top of inertia.** This materially raises the required injected force and makes the per-marker-σ refactor **non-optional for S3**, not an escalation path.

**Decision for S3:** promote **per-marker smoothed curvature** `κ_θˢ = −∇·(∇α_θˢ/|∇α_θˢ|)` (Musehane Eq. 3.7) into the S3 critical path, **or** redesign the color painting so `c` traces only one drop and **prove the film capillary force is non-attractive** before relying on it. Do not assume the single color field holds the markers apart.

**1a. Extend mTHINC to N markers — `src/simulation/m_thinc.fpp` — and guarantee closure positivity (MFC P4, P5).**

The module is hardcoded 2-fluid: one normal `mthinc_nhat(1:3,i,j,k)` and one interface position `mthinc_d(i,j,k)` (36–38), normal from the *single* first marker `eqn_idx%adv%beg` (242, 245–258), and four face fields assuming the complement of fluid-1 *is* fluid-2.

- **Storage:** promote `mthinc_nhat → mthinc_nhat(1:3, 1:num_fluids-1, i,j,k)`, `mthinc_d → mthinc_d(1:num_fluids-1, i,j,k)` (one normal + `d` per liquid marker; last fluid is the closure complement). Update allocations (199–204) and the `GPU_DECLARE` at 38.
- **`s_compute_mthinc_normals` (207–291):** wrap in `do f = 1, num_fluids-1`, `eqn_idx%adv%beg → eqn_idx%adv%beg + f - 1`, store into `mthinc_nhat(:,f,…)`, `mthinc_d(f,…)`.
- **`s_thinc_compression` — use the MTHINC branch, not plain-THINC (MFC P5).** For 3-D oriented interfaces the required path is the **MTHINC (`int_comp=2`) branch, lines 324–356**, which uses `f_mthinc_face_average` and has **no** `moncon` monotonicity test and **no** A/B/C THINC algebra (the plain-THINC `else` branch at 357 is the wrong branch for this). The **per-marker loop wraps 327–355**. For each marker compute its own `aTHINC` from its own normal/`d` and write `adv%beg+f-1`, `cont%beg+f-1` (density `rho_f = v_rs_ws(…,cont%beg+f-1)/alpha_f`). Decorate the per-marker inner loop with `GPU_LOOP(parallelism='[seq]')`.
- **Closure with an enforced positivity floor (MFC P4 — the real robustness crux, not a mechanical edit).** The naive running-sum closure `α_N = 1 − Σ_f α_f` is exactly where the scheme breaks: both liquid markers transition in the *same* film cells, mTHINC sharpens each independently (each solving its own `d` in `f_mthinc_solve_d`), and two over-steepened tanh profiles can reconstruct face values with `α₁+α₂ > 1` → `α₃_face < 0` → `ρ₃α₃ < 0` fed to the Riemann solver **right where the gas film matters most**; the clips (F9 sites above) then fire, and clipping in the film *is* the numerical merge we are trying to prevent. **Required strategy:** before forming the closure, **renormalize the liquid markers so `Σ_f α_f ≤ 1 − α_min`** (small floor), then set `α_N = 1 − Σ_f α_f ≥ α_min`, `cont%end = ρ_N·α_N`. **Prove `α₃_face ≥ α_min` at every reconstructed face** as a unit test. Give this top billing in the module changes.
- **Reshaped-vs-physical indexing must be tested in all three directions (MFC P5).** `mthinc_nhat`/`mthinc_d` are computed in **physical** `(j,k,l)` (over `q_prim_qp`) but read inside `s_thinc_compression` using the **reshaped** loop indices, with only a bare comment at line 325 and no actual remap. For `recon_dir=1` reshaped == physical so it works; the per-marker extension must be validated in **all three reconstruction directions** — callers `src/simulation/m_weno.fpp:1461` and `src/simulation/m_muscl.fpp:226` permute the reconstruction axis while the normal components stay in physical x/y/z. Treat directional correctness as an explicit test requirement, not an "applies unchanged."

**1b. Surface tension — per-marker curvature is on the S3 critical path (see the crux above).** `m_surface_tension.fpp` remains `num_fluids`-agnostic for the *outer* interfaces (capillary force from the single color field, grads 255–279, flux routine 74–233, `c_divs(num_dims+1)` at 131), and that is fine away from the film. But because the single color field is *attractive* in the film (MFC P1), S3 requires either per-marker smoothed curvature (`s_get_capillary` 236–327 refactored to per-`α_θˢ`) or a redesigned color painting with a proven non-attractive film force. This is the larger surface-tension refactor; it is **not deferred** for S3 — it is a prerequisite. (S1/S2 avoid it entirely via the symmetry plane, Section 3.9.)

**1c. Validator guards — `toolchain/mfc/case_validator.py` — current status (M6).**
- Line 772 (surface tension): **already relaxed and committed on this feature branch** (commit cda59df9) to `surface_tension and num_fluids not in (2, 3)`, as a **marked prototype edit** whose comment instructs reverting to `num_fluids == 2` for production. Phase 1a supersedes that revert note (N-fluid THINC removes the 2-fluid limitation the comment cites) — update the comment when 1a lands.
- Line 427 (int_comp): **still unrelaxed** — `int_comp != 0 and num_fluids != 2`. This is the guard that blocks N-fluid THINC; Phase 1a relaxes it (`num_fluids not in (2,3)` or a configured max) when the N-fluid extension lands.

**Phase 1 validation — two runs, one gate (M4; v2's single criterion was self-contradictory).** v2 ran the test with the force OFF but passed it on the force machinery "doing its job" — unfalsifiable as written, and its own R1 premise predicts the force-OFF run *must* merge. Split:

- **(i) Characterization run (force OFF):** initialize two markers **driven together at `U_r/2`**, film force OFF. **Expected outcome: merge.** Record the gap-collapse rate and where/when `α₃` clipping fires — this quantifies exactly how much injected force and per-marker-σ correction the machinery must supply. (If the gap survives with no force, that itself falsifies R1's premise and must be reported as such, not silently welcomed.)
- **(ii) The actual gate (force ON):** the same run with the film force + per-marker σ ON. **Pass = the markers stay distinct AND `α₃_face ≥ α_min` at every reconstructed face throughout.** Rung 7 gates on (ii), not (i).

### Phase 2 — Film detection + thickness tracking + compliance closure + source injection

**2a. Film detection — use per-marker `∇α_θ`, not the single `c_divs` (MFC P7).** For each marker compute a smoothed field `α_θˢ` (Section 2b) and detect the thin-film band:
```
β = pos(ξ·δh − h_Est) · pos(α₁ˢ − ε) · pos(α₂ˢ − ε),   ξ ≈ 3–5,  ε ≈ 10⁻³
h_Est = h_Est1 + h_Est2  from the ray integrals  ∇h_Estθ·n_sθ = (1−α_θ)  (Musehane Eqs. 3.46–3.58)
```
There is only **one** `∇c` in the code (`c_divs` holds the single color gradient + magnitude), not `∇α₁` and `∇α₂`; the detection must read the **per-marker gradients** `∇α_θ` from `q_prim_qp%vf(adv%beg+f-1)` (both present within `ξ` cells and anti-parallel), falling back to the full ray integral if too coarse. On the symmetry-plane path (S1/S2) detection is trivial (band = cells within `ξ·δh` of the plane), which is a further reason to validate the solver there first.

**2b. Smoothing solver — explicit, not implicit (MFC P7 / R8).** The per-marker smoothing `α_θˢ − α_θ = ∇²(η_s α_θˢ)` (Musehane Eq. 3.3; `η_s` is a smoothing coefficient with units of length², ~(smoothing width)² — distinct from the interface displacement `η` of Section 3.3) is written as a Helmholtz solve, but a **global implicit solve per marker per substage on GPU would be the single largest performance/complexity unknown** in the design and is out of character for an otherwise local-stencil explicit code. **Specified: an explicit smoother** — a **few (3–5) fixed Jacobi / explicit-diffusion iterations** with a documented iteration count, matrix-free, decorated with `GPU_PARALLEL_LOOP`. If accuracy demands more, escalate to matrix-free multigrid — but do not introduce a global implicit solve without accounting for it as the dominant cost.

**2c. Sub-grid film state and the once-per-step update (MFC P2).** Add device fields `film_h` (full gap `H`), `film_h_rigid`, and stored `p_film` (all `GPU_DECLARE`). **Critical placement correction stands:** the thickness field is a **history-dependent state variable**; its PDE **must be integrated once per physical step in `src/simulation/m_time_steppers.fpp`, NOT inside `s_compute_rhs`.** `s_compute_rhs` (and `s_compute_mthinc_normals` at `m_rhs.fpp:579`) is called **once per Runge–Kutta substage (3× for RK3)**; integrating a history-dependent ODE there would advance drainage **3× per step with the substage `dt`**, corrupting `h(t)`, the `h_crit` gate, and the film pressure. The split:
- **Once per step (time-stepper):** update `film_h_rigid` from the resolved `V_n`; evaluate `p_film` from the **local compliance closure** (Section 3.3 — a Laplacian stencil evaluation, *no elliptic solve*); integrate `film_h` **sub-cycled on `dt_film`** (fourth-order stability limit `∝ dx_film⁴`, Section 3.8c; semi-implicit capillary term if sub-cycling is too costly); compute the rupture min-reduction (Section 3.6); **store `p_film` as a device field.**
- **Per substage (RHS):** the momentum/energy **sources** read the stored `p_film` and are added with `+=` so they get RK-weighted correctly (below). Freezing `p_film` across the substages relaxes strict RK order for the source; this is acceptable for an operator-split source and is the pragmatic default, but flag the order caveat.

On band turn-on initialize `film_h = h_Est`; outside `β`, peg `film_h = h_Est`. Use a one-cell-eroded mask `β_reduced` for the `∇_t` terms (Musehane Eq. 3.63) to avoid a band-edge discontinuity. Local `a = r_h` and `V_n` come from the resolved interface velocities (`mom%beg..mom%end`).

**2d. Lubrication source — RHS hook, and do NOT copy the body-force zero-init (MFC P6).** Add a guarded call in the **post-flux global-source block of `src/simulation/m_rhs.fpp` (829–862)** — the block that already hosts `acoustic_source`, `bubbles_euler`, `bubbles_lagrange`, `chemistry`, `cont_damage`:
```fortran
if (film_drainage) call s_compute_film_drainage_source(q_prim_qp%vf, q_cons_qp%vf, rhs_vf)
```
Implement it **modeled on** `src/simulation/m_body_forces.fpp:79–144` (`s_compute_body_forces_rhs`), the canonical momentum+energy source pattern, **but with two copy-traps called out:**
- **Do not copy the template's zero-init loop (lines 89–99).** The template zeroes `rhs_vf(mom:E)` because body forces are assembled standalone; at the 829–862 hook **`rhs_vf` already holds the flux divergence**, so a literal copy would wipe out all transport. Use `+=` only.
- **Read `q_prim_qp%vf`/`q_cons_qp%vf`** (the hook neighborhood's arrays; pressure at index `E` in `q_prim`, total energy at `E` in `q_cons`), **not** the template's `q_prim_vf`/`q_cons_vf`.

The source, over `do l=0,p / k=0,n / j=0,m` under `GPU_PARALLEL_LOOP`, applies (Sections 3.5, 3.7):
```
rhs_vf(eqn_idx%mom%beg + d)%sf(j,k,l) += p_film · ∂α_θ/∂x_d        ! f_film = p_film ∇α_θ  (CSF interface-delta form, N/m³)
rhs_vf(eqn_idx%E)%sf(j,k,l)          += u_d · (p_film ∂α_θ/∂x_d)   ! consistent work term of the load force
rhs_vf(eqn_idx%E)%sf(j,k,l)          += ẇ_ext(j,k,l)               ! band-gas internal-energy deposit (Section 3.5)
```
where the deposit field `ẇ_ext` distributes, over the band gas cells, **exactly the mechanical power the force pair extracted**: `Σ_band ẇ_ext ΔV ≡ −Σ_band u·f_film ΔV` (`≈ ∫p_film V_n dA ≈ ∫Φ dA` in quasi-steady lubrication) — **never Φ added freely** (M5). **Unit tests:** band-integrated force `= ∫ p_film dA` (C11); global energy conservation to machine precision. All needed fields are in scope and device-resident: velocities (`mom`), pressure (`E` in `q_prim`), volume fractions (`adv`), partial densities (`cont`), color (`c`), geometry (`x_cc,dx,…`), and `film_h`/`p_film`. Enforce the magnitude guard `p_film ≲ C·ρU_r²` (Section 3.8) here.

**2e. New case parameters.** `film_drainage` (logical), `film_compliance` (logical, default T — the capillary-compliance closure; F = the kinematically-slaved-`h` fallback, Section 3.3), `h_crit`, `A_hamaker` (reserved, unused in M1), `xi_film` (ξ), `lambda_mfp` (λ = 66 nm), `f_kn_model`, `p_film_cap` (the `C·ρU_r²` guard), `film_subcycle` (sub-cycle count / `dt_film` policy). Read in `m_global_parameters.fpp` and `GPU_UPDATE(device=[…])` near line 887 (alongside `int_comp, ic_eps, ic_beta`).

### Phase 3 — Rupture criterion + merge handoff

**3a. Rupture test.** Once per step (time-stepper, alongside 2c), reduce `min_rim(film_h)` over each film band (masked min-reduction) and set a merge flag when `< h_crit` (full-gap convention, Section 3.9).

**3b. Merge handoff (Musehane fuse — default).** In flagged cells fuse the markers: `α_new = α₁ + α₂` into marker 1, zero marker 2; recompute partial densities so **mass, momentum, and energy are conserved cell-by-cell** in the fused cells (the delicate step — a topology change in a conservative compressible scheme; see R7). Turn off the injected lubrication source and dissipation in the fused band (single interface now). mTHINC then sharpens the merged interface naturally.

**3c. Alternative (Mason mode).** Instead of fusing, **switch off the lubrication source** in the ruptured band and let mTHINC / the diffuse interface reconnect across the now-thin gap. Cheaper, avoids field surgery, but relies on the auto-merge behaving (and, per MFC P1, the color force *wanting* to merge helps here). Provide both; **default to the fuse** for control, keep Mason mode as a fallback if the fuse proves non-conservative or noisy.

---

## 5. Validation ladder

Run in order; each rung gates the next. Domain `9D×9D×9D` (half-domain for the symmetry-plane path), initial center-to-center `3D`, outflow BCs on the open sides, Courant 0.8. **Head-on rungs use the reflective symmetry plane at `x=0`** (Section 3.9): free-slip BC at the plane, single drop, no markers. **tanh grid clustering at the collision plane** (`src/pre_process/m_grid.f90:51–70` serial / `178–194` parallel): the stretch map is **origin-anchored** — `x=0` is a fixed point of the map and cells grow with `|x|` — so **place the collision plane at `x=0`**; a symmetric full domain or an `x_domain%beg = 0` half domain both work (F16 — v2's "requires `x_domain%beg = −x_domain%end`" was wrong, and would have forbidden the design's own half-domain path). Tune `a_x`/`loops_x` to reach `D/200`–`D/500` at the plane, coarsening to `~D/10` in the bulk. Note (MFC P7): `dx = minval(...)` (`m_grid.f90:~68`) with `s_mpi_reduce_min` sets the global `dt` from the finest cell — clustering buys *locality*, not a larger `dt` (consistent with R4); the only relief from the finest-cell `dt` is sub-cycling the film sub-system (Section 3.8).

**Acceptance bands (M9) — explicit and provisional.** Rungs 0, 1, 2, 5, and 6 carry numeric bands; Rungs 3, 4, and 7 are shape/outcome-based (Rung 3's pass is recovering both transition boundaries in the correct We ordering — a numeric boundary tolerance is deliberately deferred until the f(Kn) sensitivity of Section 3.4 is quantified). They are **provisional**: set before any run, sized to be tight enough to falsify (a flat-film or degenerate closure fails them) yet loose enough not to punish known noise floors (the ~2% We/Re internal inconsistency of Section 1.2; graphical read-off of Figs. 15–17; property uncertainty). Revisit after Rung 1 and report any retuning explicitly — never retune to pass. Note also that Eq. 10 is fitted for `0.2 < We < 8.6`; S1 sits at the endpoint `We = 8.6` and S2 extends beyond it — state explicitly whenever a correlation is applied at or outside its fitted range.

**Rung 0 — film-solver unit tests (before any collision).** (i) **Rigid limit** (`< 1%`): with compliance off, the solver reproduces flat-plate Stefan–Reynolds `p_flat(r)` and `F_lub` to <1%. (ii) **Non-degeneracy / dimple existence (regression test for review finding M1):** a uniform-`V_n` squeeze of an initially flat film **must** develop a nonuniform, center-thick / rim-thin `h` with a rim pressure feature — v2's system fails this identically, by construction; if this fails, S2/Rungs 4–5 cannot pass, stop here. (iii) **Compliance benchmark:** compare drop–drop drainage against Yiantsios–Davis solutions; the `σ/2` compliance coefficient (Section 3.3) is provisional until this passes. (iv) **Force tests (Physics P9, C11):** on two markers, verify `f_film` anti-parallel (drops pushed apart) and band-integrated force `= ∫ p_film dA`, converging under refinement. (v) **Positivity (MFC P4):** `α₃_face ≥ α_min` across all reconstructed faces in all three directions. (vi) **Half-domain equivalence (C10):** the half-domain solver under the Section 3.9 full-gap convention matches the full two-sided flat-plate solution.

**Rung 1 — S1: Fig. 6 head-on bounce (We=8.6), symmetry-plane path.** Confirm the drop **decelerates, deforms, and rebounds** with the film gating merge. Measure `τ_b` as the interval during which `p(0,0) > p_ambient` (Chen & Yang's definition). Targets:
```
τ_c = D/U_r ≈ 306µm / 0.99 ≈ 309 µs
Eq.10: τ_b/τ_c = 0.82√8.6 = 2.40   →  τ_b ≈ 745 µs
Eq.11 form: τ_b = 0.82√(ρ_l D³/σ) ≈ 745 µs   (identical to Eq. 10 by construction when U_r is back-solved
        from We; the ~2% lives entirely in tabulated-U vs back-solved-U — see §1.2 velocity caveat)
τ_o (Eq.12) = (π/4)√(ρ_l D³/σ) ≈ 713 µs
```
**Pass band: measured `τ_b` within ±10% of Eq. 10.** **Also run the R1 experiment here (Section 3.2):** film source OFF → expect premature merge / wrong deceleration; ON → correct `h(t)`. Judge by *merge outcome and deceleration*, not by whether it rebounds (the capillary spring rebounds it either way on this path).

**Rung 2 — near-contact resolution study, reconciled with the film-solver well-posedness tension (MFC P7 / R5).** Rerun S1 at plane resolutions `D/100, D/200, D/500`. **Pass = bounce outcome unchanged and successive-resolution change in `τ_b` < 5%** — the load-bearing verification that the model *works*: unlike Chen & Yang (where outcome is a deliberate function of min grid size), a correct sub-grid model must be **grid-independent**. **But note the real tension:** at `D/200` the max-deformation film (~0.3 µm) is 0–1 cells thick, and the in-plane fourth-order film stencil needs multiple cells — "fall back to the ray integral" does not create resolution that isn't there. So this rung has a second job: **determine whether `D/200` even makes the film solver well-posed, or whether `D/500` (4× the steps) is required** purely for stencil support. Decide `D/200` vs `D/500` here; this may change the cost story from `~10⁶` to `~4×10⁶` steps, and that is an honest possibility, not a failure.

**Rung 3 — S2: nonmonotonic head-on transition, symmetry-plane path.** Sweep `We` across the bouncing window, **anchored to `We = 1.88` (merge→bounce) and `We = 10.3` (bounce→merge)** — valid at the target tetradecane/N₂/1 atm condition per Chen & Yang p. 062103-12/13, citing Jiang, Umemura & Law (1992) (F7). Remember Section 1.3: Chen & Yang's Fig. 2 boundaries *are* Qian & Law's curves and their own low-We runs disabled merging numerically, so the experimental values — from Jiang et al. (1992), a different experiment than the Fig. 2 curves, hence independent of them — are the primary anchor and Fig. 2's simulation symbols the only independent numerical check. Reproduce both boundaries: the lower must emerge from the drainage gate (`t_drain ∝ a⁴`); the **upper must emerge from rim thinning** (`min_rim h(t) < h_crit` while pressed, Section 3.6) — **not** from `t_drain ∝ a⁴` and **not** from grid tuning. **Report S2 sensitivity to `f(Kn)` and to `λ`** (Physics P7) as first-class results. Cross-check the `r_h`↑ / `h_min`-nonmonotonic trends (Chen & Yang Fig. 17).

**Rung 4 — diagnostic quartet.** For S1 and the S2 cases extract `h_center`, `h_min`, `r_h`, `p(0,0)` vs `t/(D/U)` and compare *shapes* to Chen & Yang Figs. 13, 15: `h_min` decreasing faster than `h_center`, the **nonmonotonic `h_min`** (rapid rise then sudden drop at rebound), and the **rim pressure peak / negative rim gauge pressure during bounce-off** (Fig. 15b, t5–t7). These shapes are the physics and are the direct test that the compliance closure produces a real dimple; matching them is stronger evidence than matching `τ_b`. (A flat-film closure — or v2's degenerate system — would fail this rung, which is why the closure was replaced.)

**Rung 5 — Fig. 16 pressure collapse (also the runtime stability guard).** Confirm the film-center pressure is roughly constant during a bounce and **collapses when normalized by `ρU_r²`** across We — **pass band: collapse within ±50% across the We sweep** (the plateau sits at `p/(ρU_r²) ≈ 5`). Per Section 3.8(b) this doubles as a **runtime assertion**: if injected `p_film` does not collapse on `ρU_r²`, the run is unphysical.

**Rung 6 — restitution (Eq. 13) — the dissipation test (Physics P6).** Confirm `τ_b/τ_o = 0.29 We^(−0.35) + 0.88` — **pass band: measured `τ_b/τ_o` within ±5% of Eq. 13.** This directly checks the **sub-grid dissipation** and the conservative energy accounting (Section 3.5); an artificially elastic bounce (work terms without the extraction-deposit pairing) shows up as `τ_b` too long / restitution too high here, and energy creation (M5's failure mode) shows up as the opposite bias.

**Rung 7 — S3: off-center bounce, full multi-marker.** Full 3-D, **no symmetry plane** (the only stage that needs the multi-marker machinery, Phases 1–3). **Primary case (M3): `B = 0.9, We = 48.8, Re = 260.3` — the paper's only off-center bouncing simulation (Fig. 10), whose rotating film is Fig. 11 (same case).** Confirm off-center bounce with the film **rotating about the collision center** and drops separating before rupture; compare the film-rotation and separation sequence against Figs. 10–11. **Secondary check: `B = 0.6`** as a regime-map bounce/no-bounce outcome only (Fig. 2). Gated on Phase 1's gate **(ii)** — the force-ON driven-together run — passing (M4). Recall Fig. 10's own resolution caveat (Section 2.3): it is an `O(10⁻³)D` simulation whose non-merging is grid-guaranteed, so shape/rotation comparison is the meaningful content, not film-scale quantities.

**Bouncing-time correlations to report** (Eqs. 10–13): Eq. 10 `τ_b/τ_c = 0.82√We`; Eq. 12 `τ_o = (π/4)√(ρ_l D³/σ)`; Eq. 13 `τ_b/τ_o = 0.29 We^(−0.35) + 0.88`.

---

## 6. Risks, open questions, and decision points

**R1 — Is the injected force necessary, and for what? (Physics P2)** The force's job is **load transmission / correct deceleration and prevention of premature numerical merge**, *not* rebound (rebound is the resolved capillary spring). **Decision gate:** in Phase 0, run S1 with drop + resolved gas only (force OFF) and measure **whether the drops merge prematurely and whether the deceleration and `h(t)` are correct** — *not* whether `τ_b` is under-predicted (that is a capillary test). Cheap, run early. Honest uncertainty: we do not know a priori how much load the `D/200` cushion already transmits — which is also why the double-counting flag (Section 3.5) exists.

**R2 — Continuum lubrication breaks down at the rupture scale.** At `h_crit ≈ 30 nm`, `Kn ≈ 2.2` — transition/free-molecular; the closure is a continuum extrapolation there. **Mitigation:** the *load-bearing* pressure is built at the dimple (`h ~ 0.3 µm`, `Kn ≈ 0.22`), where continuum+slip is defensible; `h_crit` only *gates* merging. So the force is on firmer ground than the final-rupture instant. Still, `f(Kn)` must be sourced properly (R-Kn below).

**R3 — `h_crit` calibration.** Sensitivity is mild above ~50 nm and sharpens below (Mason Fig. 6). Calibrate `h_crit` (default 30 nm) so S2 reproduces both transition boundaries; do **not** tune it per-case (that would hardcode an outcome).

**R4 — Timestep cost of near-contact refinement is real.** Refining the plane to `D/200`–`D/500` sets the global acoustic `dt` (Section 2.3): ~10⁶ steps at `D/200`, ~4×10⁶ at `D/500`. **Coupled with R5 and R-stiff below**, the *effective* cost may be higher: `D/500` may be forced for well-posedness (R5) and/or the film source may need sub-cycling (R-stiff). Stop at the coarsest resolution that passes Rung 2.

**R5 — Film-region representation at `D/200` is a likely failure, not a tail risk (MFC P7).** At max deformation the film (~0.3 µm) is below one `D/200` cell (1.5 µm), so `β` may be 0–1 cells thick — and the v3 fourth-order in-plane operator needs an even wider stencil (≥2–3 cells each side) than v2's second-order one. "Fall back to the ray integral" does not create resolution. The tension between "stop at `D/200` for `dt`" (R4) and "need an in-plane stencil in the film" is **unresolved by design** and may force `D/500`. Rung 2 decides it; watch for band dropout (film not detected). If the compliance closure itself proves unusable on the affordable band, the documented fallback is the kinematically-slaved-`h` mode (Section 3.3) — with its stated cost to the rim rupture gate.

**R6 — Compressible-framework subtleties.** (a) **Energy accounting (Section 3.5):** the extraction-deposit pairing makes total energy conservative *by construction*; the remaining delicacy is the KE↔internal partition inside MFC's single-velocity mixture, and Eq. 13 (Rung 6) arbitrates it. The double-counting overlap with the resolved cushion is an open calibration point. (b) Gas compressibility in the film is negligible (`ΔP ~ 10² Pa`) — the v3 film model is incompressible by construction. (c) The single color function is *attractive* in the film (MFC P1); per-marker-σ / color-repaint is a **prerequisite for S3**, avoided for S1/S2 via the symmetry plane.

**R7 — Merge handoff conservation.** Fusing markers (Phase 3b) is a topology change in a conservative compressible scheme; naive fusion can violate mass/energy conservation or inject acoustic transients. Verify conservation cell-by-cell; fall back to Mason-mode (3c) if fragile. Most numerically delicate step.

**R8 — GPU cost (updated: the elliptic solve is gone).** v3's compliance closure makes `p_film` a **local stencil evaluation** — the global elliptic band solve that dominated v2's cost projection **disappears** (Section 3.3) — though note the kinematically-slaved *fallback* mode reintroduces a rigid-limit elliptic problem, so its cost story differs. The remaining unknowns are the **per-marker explicit smoother** (Section 2b — kept explicit precisely to stay in character with the code) and the **sub-cycled fourth-order film integrator** (Section 3.8c), plus `film_h`/`p_film` fields, masked min-reductions, and the source kernel. Budget is uncertain; profile after Phase 2. All new arrays/kernels/scalars need `GPU_DECLARE`/`GPU_PARALLEL_LOOP`/`GPU_UPDATE` or the run silently falls off the device.

**R9 — Does multi-marker with a single color field hold the drops apart? (MFC P1 — the central S3 risk.)** No, not intrinsically: single-pressure/single-velocity gives zero relative slip; distinctness survives only while a ≥2–3-cell resolved gas gap persists, and the color force actively *pulls the film shut*. When the gap closes, the positivity clips silently redistribute the markers (F9: `m_riemann_solver_hllc.fpp:1238–1258` face-state clip+renormalize at any `num_fluids`; `m_pressure_relaxation.fpp:111–138` mutating `q_cons_vf`; `m_variables_conversion.fpp:285–293`) — clipping in the film *is* the numerical merge. **This is why (i) the symmetry plane carries S1/S2 to avoid the risk entirely, (ii) per-marker-σ / color-repaint is a prerequisite for S3, and (iii) Phase 1's validation is a driven-together pair — characterization with the force off, gate with it on (M4).** Second decision gate; surfaced at Phase-1 validation.

**R-Kn — `f(Kn)` directly moves the S2 thresholds (Physics P7).** `f(Kn)=1+6Kn` is valid only for `Kn ≲ 0.1`; at the dimple `Kn≈0.22`, at rupture `Kn≈2.2`. Because `μ_eff=μ_g/f` rescales both the force and the drainage flux, it *sets* the transition `We`. **Obtain the real `Q(Kn)`/`f(Kn)` (Zhang & Law, or the 2025 rarefied dual-VOF/augmented-vdW papers) before S2, and report S2 sensitivity to `f(Kn)` and to `λ` (66 nm vs the paper's ~0.1 µm annotation).**

**R-stiff — the binding timestep near contact may be source stiffness, not acoustics (Physics P5; MFC P3).** Two stiff mechanisms: the mobility's `1/h³` pressure sensitivity over two decades of `h`, and the v3 fourth-order capillary term with its `dt_film ∝ dx_film⁴` explicit limit. **Mitigation:** two-way self-limiting `V_n` (Section 3.8a), a `p_film ≲ C·ρU_r²` runtime guard (3.8b), and **sub-cycling / semi-implicit treatment on a `dt_film` clock** (3.8c). The Section-2 `~10⁶`-step budget assumes the film is sub-cycled off the global clock; without that, the global `dt` collapses.

**R-close — N-fluid closure positivity (MFC P4).** Independent per-marker sharpening can drive `α₃_face < 0` in the film, whose clip *is* the numerical merge. **Mitigation:** renormalize liquid markers to `Σα ≤ 1−α_min` before closure and prove `α₃_face ≥ α_min` (Rung 0-v). This is the crux robustness problem of Phase 1, not a mechanical index edit.

---

## 7. References

**Primary target (physics + validation):**
- Chen, X. & Yang, V. (2020). *Direct numerical simulation of multiscale flow physics of binary droplet collision.* Physics of Fluids **32**, 062103. (Resolved-film AMR-DNS baseline for head-on; off-center case at ~0.7 µm minimum grid with grid-size-determined merging distance — Section 2.3. Figs. 2, 6, 7, 10–11 (off-center primary S3 case), 13, 15–17; Eqs. 10–13; `We = 1.88/10.3` anchors at the target condition. PDF archived at `~/wpi/demulsification/gas_film_drainage/chen_yang_2020.pdf`.)

**Sub-grid film models (design backbone):**
- Musehane, N. M. (2016). *Direct Numerical Simulation of Bubble–Bubble and Droplet–Droplet Interaction Using a Surface Thin Film Model.* MSc thesis, University of Cape Town. (Multiple-marker + Reynolds film; Eqs. 3.2–3.10, 3.19–3.20, 3.43–3.69.)
- Musehane, Oxtoby & Reddy (2018). *Multi-scale simulation of droplet–droplet interaction and coalescence.* Journal of Computational Physics.
- Coyajee & Boersma (2009). Multiple-marker method (credited origin). Journal of Computational Physics **228**, 4444–4467. (Not read directly.)
- Mason, Stevens & Harvie (2012). *Multi-scale volume of fluid modelling of droplet coalescence.* Ninth International Conference on CFD in the Minerals and Process Industries, CSIRO, paper 179MAS. (Quasi-1D film + two-way body force; Eqs. 5–10.)
- Harvie & Fletcher (2001); Harvie (1999). (Source of Mason's `dh/dt` closure; not retrieved.)

**Film-drainage theory:**
- Yiantsios, S. G. & Davis, R. H. (1990). *On the buoyancy-driven motion of a drop towards a rigid surface or a deformable interface.* Journal of Fluid Mechanics **217**, 547–573. (**The capillary-compliance drainage/dimple mechanism the v3 closure implements; Rung 0-iii benchmark.**)
- Yiantsios & Davis (1991). Drop–drop companion paper, Journal of Colloid and Interface Science **144**, 412. (Not read directly; benchmark solutions for the two-drop compliance coefficient.)
- Zhang & Law (2011). *An analysis of head-on droplet collision with large deformation in gaseous medium.* Physics of Fluids **23**, 042102. (Deformation-energy + rarefied lubrication + vdW rupture; nonmonotonic transition. **`Q(Kn)`/`f(Kn)` paywalled — must be obtained before S2, R-Kn.**)
- Kwakkel, Breugem & Boersma (2013). *Extension of a CLSVOF method for droplet-laden flows with a coalescence/breakup model.* Journal of Computational Physics **253**, 166–188. (Multiple-marker CLSVOF + Zhang–Law drainage; `t_c ≥ t_drain` gate.) [Volume/pages corrected in v3 — M8.]
- He, C., Xia, X. & Zhang, P. (2020). Physics of Fluids **32**, 032004. (Two-VOF-function, marker-per-drop anti-merge precedent for off-center bouncing — Chen & Yang's ref 38; no rupture criterion, so it "always prompts bouncing"; not read directly, characterized via Chen & Yang. Section 2.3.)
- Rajkotwala et al. (local front reconstruction / bouncing droplets — citation incomplete: year/venue not pinned; not read directly); Jiang & James (2007); Yoon et al. (vdW / disjoining-pressure approaches — citation incomplete; documented failure modes at low Ca / final drainage).
- Chen & Wu (2022), arXiv:2205.03604 (rarefaction reframing of Zhang & Law; local PDF).
- Deblais et al. (2024), arXiv:2402.00500 (interferometric jump-to-contact ≈30 nm; `h_crit` calibration; local PDF).
- Sprittles, *Gas Microfilms in Droplet Dynamics: When Do Drops Bounce?*, Annual Review of Fluid Mechanics (review; paywalled, not read).
- Qian & Law (1997). Journal of Fluid Mechanics **331**, 59. (**Not read directly** — source of Chen & Yang's Fig. 2 solid regime-boundary curves and of the gas-pressure/composition sensitivity of the bouncing window. **Not** the source of the `We = 1.88/10.3` anchors — those are Jiang, Umemura & Law (1992), below. Circularity note (M10, softened): Chen & Yang's Fig. 2 solid boundaries *are* Qian & Law's curves, so the independent S2 check against Fig. 2 is Chen & Yang's simulation symbols, not the curves; the `1.88/10.3` experimental anchors, being from a different experiment, are independent of the curves.)
- Jiang, Y. J., Umemura, A. & Law, C. K. (1992). *An experimental investigation on the collision behaviour of hydrocarbon droplets.* Journal of Fluid Mechanics **234**, 171–190. (**Not read directly** — primary experimental source of the `We = 1.88` (merge→bounce) and `We = 10.3` (bounce→merge) S2 anchors for head-on tetradecane in N₂ at 1 atm, taken via Chen & Yang's citation — their ref 7, p. 062103-12 — F7.)
- Dual-VOF with rarefied gas + augmented vdW (2025), International Journal of Multiphase Flow S0301932225000850; *Acta Mechanica Sinica* (2025) 10.1007/s10409-025-25492-x (modern rarefied + augmented-vdW closures; recommended source for the exact `f(Kn)` and disjoining terms; paywalled).

**MFC implementation hooks (file:subroutine — verified against source in review; F9/F8 corrections applied):**
- `src/simulation/m_thinc.fpp` — `s_compute_mthinc_normals` (207–291), `s_thinc_compression` (295–406), **MTHINC branch 324–356 / per-marker loop 327–355** (int_comp=2; *not* the plain-THINC 357–395 branch), storage (36–38, 199–204). *N-fluid extension + closure positivity + reshaped-index test.*
- `src/simulation/m_weno.fpp:1461`, `src/simulation/m_muscl.fpp:226` — reconstruction callers (directional/reshaped-index correctness).
- **Positivity-clip sites (the silent marker-fusion mechanism, R9/R-close — F9 corrected):** `src/simulation/m_riemann_solver_hllc.fpp:1238–1258` (`mpp_lim` face-state clip+renormalize in the plain 5-equation branch, any `num_fluids`); `src/simulation/m_pressure_relaxation.fpp:111–138` (`s_correct_volume_fractions`, mutates `q_cons_vf`); `src/common/m_variables_conversion.fpp:285–293`. (v2 cited `m_riemann_solver_hllc.fpp:839` — that branch does no clipping and only runs with `bubbles_euler`.)
- `src/common/m_derived_types.fpp:137–155` (`eqn_idx_info`); `src/common/m_global_parameters_common.fpp:161–215` (5-eqn layout), **single color `c` at 296**; `src/pre_process/m_assign_variables.fpp:545` (`c` painted from `patch_icpp%cf_val`). *Index ranges + the coalescence-promoting color-force problem.*
- `src/simulation/m_surface_tension.fpp` — `s_get_capillary` (236–327), `s_compute_capillary_source_flux` (74–233), color grads 255–279, `c_divs(num_dims+1)` at 131. *Per-marker-σ prerequisite for S3.*
- `src/simulation/m_time_steppers.fpp` — **once-per-step** film state update + compliance-closure `p_film` evaluation + rupture min-reduction (MFC P2; must NOT live in the per-substage RHS).
- `src/simulation/m_rhs.fpp:829–862` — post-flux source hook (per-substage); `m_rhs.fpp:579` — per-substage `s_compute_mthinc_normals` call; template `src/simulation/m_body_forces.fpp:79–144` — **do not copy the 89–99 zero-init**; use `q_prim_qp%vf`/`q_cons_qp%vf`.
- `toolchain/mfc/case_validator.py:427` (int_comp — still unrelaxed), `:772` (surface_tension — relaxed, committed at cda59df9). *Guard status per Phase 1c (M6).*
- `src/pre_process/m_grid.f90:51–70` / `178–194` (tanh clustering — origin-anchored map, F16), `~68` (`dx = minval`, global-dt-from-finest-cell).
- `src/common/include/parallel_macros.fpp` — `GPU_PARALLEL_LOOP` (25), `GPU_DECLARE` (114), `GPU_UPDATE`, `GPU_ROUTINE` (57), `GPU_LOOP` (132). *Offload decoration for all new kernels/fields/scalars; the smoother is an explicit few-iteration solve and the compliance closure is a local stencil — no global implicit solves (R8).*
