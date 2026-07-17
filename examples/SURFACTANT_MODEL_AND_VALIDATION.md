# Surfactant transport in MFC: the model, its discretization, and what the examples actually establish

This document is the explanation to give when someone asks "how do you know this is right?" It runs in a
straight line: the equation we are solving → why it cannot be coded literally → the reformulation that makes
it work → how it is discretized → what each example does and does not prove.

**The one-sentence summary.** MFC transports a *volumetric* surfactant density `Γ̃ = Γ|∇c|` conservatively;
this reproduces the surface transport equation without ever building a surface mesh; the surface-diffusion
operator is validated quantitatively against an exact analytic solution (to ~1.4%), and the flow coupling is
validated qualitatively against three published drop benchmarks.

---

## 1. The physics we are modeling

Surfactant is a surface-active contaminant that sits *on* a drop interface. It does two things:

1. It **lowers the surface tension** `σ` wherever it collects.
2. Because the flow **sweeps it around the interface**, `σ` becomes *non-uniform* — and a tension gradient
   pulls tangentially along the surface. That tangential pull is the **Marangoni stress**.

This matters for droplet coalescence: surfactant swept into the gap between two approaching drops raises the
local tension gradient and resists drainage, changing whether and when they merge.

## 2. The equation being solved (sharp-interface form)

For an insoluble surfactant of surface concentration `Γ` (mass per unit *area*) on an interface moving with
velocity `u`, the transport equation is (Stone 1990):

```
∂Γ/∂t  +  ∇ₛ·(Γ uₛ)  +  Γ κ (u·n)   =   D_s ∇ₛ²Γ
   (1)        (2)            (3)              (4)
```

- **(2) Surface convection** — surfactant slides *along* the interface with the tangential velocity `uₛ`.
- **(3) Dilatation** — when the interface **stretches**, the same surfactant spreads over more area, so `Γ`
  drops. `κ` is the curvature; `u·n` the normal velocity.
- **(4) Surface diffusion** — `∇ₛ²` is the Laplace–Beltrami operator: diffusion *constrained to the surface*.
  `D_s` is the surface diffusivity. The surface Péclet number is `Pe = γ̇R²/D_s`.

Every symbol with a subscript `s` lives *on the surface*. That is the whole problem, as the next section
explains.

## 3. Why that equation cannot be coded literally in MFC

MFC is a **diffuse-interface** solver. There is no surface mesh, no list of surface elements, no `n` or `κ`
as primary data. The interface is represented by a **color function** `c` that goes smoothly from 0 to 1 over
a few cells. There is no surface to integrate `∇ₛ·` over.

So terms (2), (3), (4) cannot be evaluated directly. Something else is needed.

## 4. The reformulation: transport a volume density instead

The trick (Teigen 2009/2011; James & Lowengrub 2004; Jain 2024) is to stop tracking `Γ` and instead track

```
Γ̃  =  Γ |∇c|          ("surfactant per unit VOLUME", smeared across the interface band)
```

`|∇c|` acts as a **surface delta function**: it is ~zero away from the interface and peaks on it. It
converts a surface integral into a volume integral:

```
∫ Γ dS   =   ∫ Γ|∇c| dV   =   ∫ Γ̃ dV
```

**The key result**, and the reason the whole scheme works:

> `Γ̃` obeys a *plain conservative advection equation in the volume*:
> ```
> ∂Γ̃/∂t + ∇·(Γ̃ u) = 0
> ```

Terms (2) **and** (3) come out automatically. There is **no curvature term and no `∇ₛ·uₛ` term anywhere in
the code**. This works because `|∇c|` itself changes as the interface stretches: stretch the interface, and
`|∇c|` falls in exactly the way needed to make `Γ = Γ̃/|∇c|` dilute correctly. The stretching is captured
*structurally*, not by a term someone wrote down.

The physical concentration is recovered by dividing back:

```
Γ  =  Γ̃ / |∇c|        on the band where |∇c| > 1e-6 (capillary_cutoff)
```

This recovery is used only where `Γ` is *needed* — i.e. in the `σ(Γ)` closure.

## 5. How it is discretized

**Advection.** `Γ̃` is carried as one extra **conserved scalar** (`eqn_idx%surf`). It is not density-weighted:
its conservative and primitive forms are identical. It is reconstructed with **WENO5** and fluxed by the same
**HLLC** Riemann solver as the rest of the state, advected at the **contact-wave speed** `s_S`:

```
F_surf = (upwinded Γ̃) · s_S                        [m_riemann_solver_hllc.fpp]
rhs(surf) += ( F_{j-1} − F_j ) / dx                [m_rhs.fpp]
```

The update is a **pure flux difference** — a discrete divergence. This has a consequence that must be stated
plainly:

> **Because the update is a flux difference, the total `∑Γ̃` can only change through the domain boundaries.
> Conservation to machine precision is therefore guaranteed by construction. It confirms the implementation
> is correct; it is NOT evidence that the physics is right.** Any conservative scheme, right or wrong, does
> this. Do not present it as validation.

The color function `c` is advected like a volume fraction (with the usual non-conservative correction).

**Surface diffusion.** This was the hardest piece and the one place a real bug was found and fixed:

| attempt | form | outcome |
|---|---|---|
| naive projection | `D_s ∇·((I − n⊗n)∇Γ̃)` | **leaked** `Γ̃` off curved interfaces — band width grew 0.044→0.168 |
| Rätz–Voigt | `D_s ∇·(|∇c|(I − n⊗n)∇(Γ̃/|∇c|))` | leak gone, but the `Γ̃/|∇c|` division is fragile → rate only 0.58× exact |
| **Jain 2024, Eq. 6** *(implemented)* | `D_s( ∇Γ̃ − 2(0.5 − c)·n·Γ̃/ε )`, `ε ≈ Δx` | **correct** — 0.986× exact |

The implemented flux has two parts:

- `D_s ∇Γ̃` — ordinary **isotropic** diffusion of the smeared density.
- `− 2 D_s (0.5 − c) n Γ̃ / ε` — a **sharpening flux**. Note `(0.5 − c)` is `+0.5` outside the drop and `−0.5`
  inside, so this term pushes `Γ̃` *back onto the `c = 0.5` contour from both sides*.

Isotropic diffusion would smear surfactant off the interface; the sharpening flux re-confines it. Their sum
reproduces true Laplace–Beltrami surface diffusion **with no projection tensor and no division by `|∇c|`** —
which is precisely why it is robust where the first two attempts were not. It is written in divergence form,
so it also conserves `Γ̃`.

**The `σ(Γ)` closure.**

| `sigma_model` | closure |
|---|---|
| 0 | `σ` constant |
| 1 | `σ + (dσ/dT)(T − T_ref)` — thermal Marangoni |
| 2 | `σ + (dσ/dΓ)·Γ` — linear solutocapillary |
| 3 | `σ₀(1 + E·ln(1 − Γ/Γ∞))` — **nonlinear Langmuir** (`sigma_El` = `E`, `surf_max` = `Γ∞`) |

`σ` is floored at `1e-3·σ₀` so that Langmuir saturation (`Γ → Γ∞` drives `σ → 0`) cannot crash the capillary
force.

**Where the Marangoni force comes from — this is the elegant part.** There is **no Marangoni source term in
the code**. The capillary stress tensor is

```
T = σ (I − n⊗n) |∇c|                (the solver stores its flux form Ω = −T; Schmidmayer et al. 2017)
```

Taking `∇·T` with a *spatially varying* `σ` produces **both** the normal (Laplace) force **and** the
tangential (Marangoni) force automatically. Making `σ` a function of `Γ` is therefore the *entire* coupling:
nothing else had to be added.

---

## 6. What each example establishes — tiered by strength of evidence

The evidence is **not all the same quality**. Presenting it as though it were is the single easiest way to
lose credibility. Use these three tiers.

### Tier 1 — validated against an exact analytic solution *(the strongest evidence)*

**`1D_/2D_/3D_solutocapillary_diffusion`** — Laplace–Beltrami eigenmode decay.

A surfactant pattern that is an eigenmode of the surface Laplacian decays at a rate known in closed form:

| geometry | pattern | exact decay rate | measured |
|---|---|---|---|
| flat | `cos(kx)` | `D_s k²` | **−0.1%** |
| circle | `cos(mθ)` | `m² D_s / R²` | **0.986×** |
| sphere | `Y_lm` | `l(l+1) D_s / R²` | **0.949×** |

**This is the only quantitative validation in the set, and it is the one to lead with.** It is what proves the
surface-diffusion operator is *the operator we claim it is* — and it is what caught the original leaking
implementation.

### Tier 2 — validated against published benchmarks, qualitatively

These reproduce the **directions and orderings** the benchmark papers establish — which is what those papers'
acceptance gates ask for — **not their numbers**. Section 7 explains why not.

| example | benchmark | result |
|---|---|---|
| **M1** `2D_Xu2006_surfactant_shear` | Xu et al. 2006, Table 3.1 | `D` rises with coverage (0.563 → 0.566 → 0.572); surfactant sweeps to the tips; `σ` minimal there |
| **M2** `2D_Xu2012_surfactant_sweep` | Xu et al. 2012, Table 3.2 | `Ca↑→D↑` (0.455/0.608/0.706); `λ↑→D↓` (0.612/0.608/0.582); `Re↑→D↑` (0.608/0.652); `Pe↑→` surfactant more tip-concentrated (2.91/3.38/3.46) |
| **M3** `2D_PimentaOliveira_rheology` | Pimenta & Oliveira 2021, Table 3.3 | `[η_m] = 0` exactly for a clean drop; grows with coverage (0 → 0.018 → 0.054); `N₁ > 0`; `[η] = [η_c] + [η_m]` |

The sharpest single result here is **M3's `[η_m] = 0` for a clean drop**: with no surfactant there are no
tension gradients, so the Marangoni stress must vanish identically — and it does, with no tuning.

### Tier 3 — consistency and characterization, **no reference solution**

**`2D_solutocapillary_transport` (M0)** — this is *not* a validation, and its README says so explicitly.
It provides:

- mass conserved to 1.00000 (self-consistency — see the warning in §5),
- numerical diffusion *quantified* (+11% band broadening per lap; 10.3% `m=1` decay per rotation at `R/dx≈26`),
- transport-follows-flow agreement to ~1% (a real cross-check: momentum and surfactant are independently
  evolved fields, and `Γ` moves at exactly the rate the velocity field says it should),
- a documented **negative result**: solid-body rotation fails as an exact reference (the *flow* decays 13.5%
  per rotation, boundary-driven; the surfactant rides it faithfully). See that README for why this is
  fundamental rather than a tuning problem.

**Not validation at all:** `2D_/3D_solutocapillary_marangoni` are qualitative demos.

---

## 7. The limits — state these before you are asked

All four are consequences of MFC being an **explicit compressible** solver. None is a bug.

1. **We cannot reach the Stokes limit.** `Re → 0` requires kinematic viscosity `ν = μ/ρ → ∞`, and the explicit
   viscous time step `~Δx²/ν` then vanishes. Every drop benchmark (M1/M2/M3) runs at **finite `Re ≈ 1`**,
   where the papers use `Re ≈ 0`. This is the main reason we match trends, not numbers.
2. **We cannot script a velocity field.** MFC always solves the full momentum equation, so the classic
   imposed-velocity transport benchmarks (Xu & Zhao 2003; Jain 2024 §6.2) are a *different problem*, not a
   harder one. This is why M0 has no reference.
3. **2D diffuse-interface drops over-deform**, and Langmuir saturation destabilizes under-resolved sharp
   tips — so `Ca` is capped near 0.3 rather than Xu's 0.7.
4. **Low density blows the acoustic CFL** (`c ~ 1/√ρ`), which bounds how `Re` can be varied.

---

## 8. How to say it — honest phrasings

The results are good. They get *undermined* by overclaiming, so use the right-hand column.

| result | ❌ do not say | ✅ say |
|---|---|---|
| mass = 1.00000 | "conservation validates the transport" | "confirms we implemented the conservative form correctly — it's guaranteed by construction, so it's a code check, not physics" |
| diffusion rate 0.986× | "the model is validated" | "**the surface-diffusion operator reproduces the exact Laplace–Beltrami rate to 1.4%** — this is our quantitative validation, and it's what caught the original leaking operator" |
| `D` rises with coverage (M1) | "matches Xu 2006" | "reproduces Xu's *trend*; we're at finite Re, so we match orderings, not their deformation values" |
| M2 sweeps | "validated against Xu 2012" | "all four of Xu 2012's qualitative response gates come out in the right direction" |
| `[η_m] = 0` clean (M3) | "matches Pimenta" | "the Marangoni stress vanishes identically without surfactant and grows with coverage — the decomposition behaves structurally as Pimenta & Oliveira require" |
| M0 transport | "validates transport" | "consistency checks and a numerical-diffusion characterization; **there is no reference solution — MFC can't impose a velocity field**" |
| overall | "the surfactant model is validated" | "the surface-diffusion operator is validated quantitatively against theory; the flow coupling is validated qualitatively against three published benchmarks; the transport is self-consistent but has no reference" |

**A finding worth volunteering, not hiding:** the guide's shorthand pairs "Pe↑" with "more uniform
surfactant." That is backwards — higher `Pe` means *weaker* surface diffusion and therefore a *more*
tip-concentrated interface. We measured it (2.91 → 3.38 → 3.46 as `Pe` goes 1 → 10 → 100) rather than
asserting a direction. Volunteering this demonstrates the sweeps were actually run.

## 9. The paragraph for a committee

> Surfactant is transported in MFC as a conserved volumetric density `Γ̃ = Γ|∇c|`, following Teigen and
> Jain. Advecting `Γ̃` conservatively reproduces both surface convection and interface stretching without any
> explicit curvature term, and the physical concentration is recovered as `Γ = Γ̃/|∇c|` on the interface band.
> Surface diffusion uses Jain's (2024) sharpening-flux formulation, which we adopted after finding that the
> naive projected operator leaks surfactant off curved interfaces; the corrected operator reproduces the exact
> Laplace–Beltrami decay rate to within 1.4% in 2D and 5% on a sphere. Coupling to the flow requires no source
> term: making `σ` depend on `Γ` inside the capillary stress tensor produces the Marangoni force
> automatically. The coupled model reproduces the qualitative gates of three published benchmarks — Xu et al.
> (2006, 2012) and Pimenta & Oliveira (2021) — including the requirement that the Marangoni contribution to
> the bulk stress vanish identically for a clean drop. Because MFC is an explicit compressible solver, these
> comparisons are at finite Reynolds number rather than in the Stokes limit, so they establish trends and
> orderings rather than the benchmarks' quantitative values.

## 10. References

**The full, CrossRef-verified reference list lives in the
[1D README](1D_solutocapillary_diffusion/README.md#references)** — cite from there, not from here. The key
attributions:

- **Transport equation** — Stone (1990), *Phys. Fluids A* **2**(1), 111–112.
- **Diffuse-interface surfactant (`Γ̃ = Γ|∇c|`, conserve-then-recover)** — Teigen et al. (2009), *Commun.
  Math. Sci.* **7**(4); Teigen et al. (2011), *JCP* **230**(2), 375–393 (closest precedent); James &
  Lowengrub (2004), *JCP* **201**(2), 685–722. *Note: our exact recipe is a consolidation of these two, not
  lifted verbatim — say so if asked about provenance.*
- **Surface diffusion / sharpening flux** — Jain (2024), *JCP* **515**, 113277, Eq. (6); Rätz & Voigt (2006),
  *Commun. Math. Sci.* **4**(3).
- **CSF force** — Brackbill, Kothe & Zemach (1992), *JCP* **100**(2), 335–354; the capillary-stress-tensor
  implementation in MFC follows Schmidmayer et al. (2017).
- **`σ(Γ)` Langmuir EOS** — Pawar & Stebe (1996), *Phys. Fluids* **8**(7), 1738–1751.
- **Benchmarks** — Xu, Li, Lowengrub & Zhao (2006), *JCP* **212**(2), 590–616 (M1); Xu, Yang & Lowengrub
  (2012) (M2); Pimenta & Oliveira (2021) (M3); Stone & Leal (1990), *JFM* **220**, 161–186.
- **Eigenmode-decay validation practice** — Dziuk & Elliott (2013), *Acta Numerica* **22**, 289–396;
  Macdonald, Brandman & Ruuth (2011), *JCP* **230**(22).
