# Surfactant transport in MFC: the model, its discretization, and what the examples actually establish

This document is the explanation to give when someone asks "how do you know this is right?" It runs in a
straight line: the equation we solve → why it cannot be coded literally → the reformulation that makes it
work → how it is discretized → what each example does and does not prove.

**The one-sentence summary.** MFC transports a *volumetric* surfactant density $\tilde\Gamma = \Gamma\lvert\nabla c\rvert$
conservatively; this reproduces the surface transport equation without ever building a surface mesh; the
surface-diffusion operator is validated quantitatively against an exact analytic solution (to $\sim 1.4\%$),
and the flow coupling is validated qualitatively against three published drop benchmarks.

---

## 1. The physics we are modeling

Surfactant is a surface-active contaminant that sits *on* a drop interface. It does two things:

1. It **lowers the surface tension** $\sigma$ wherever it collects.
2. Because the flow **sweeps it around the interface**, $\sigma$ becomes *non-uniform* — and a tension
   gradient pulls tangentially along the surface. That tangential pull is the **Marangoni stress**.

This matters for droplet coalescence: surfactant swept into the gap between two approaching drops raises the
local tension gradient and resists drainage, changing whether and when they merge.

## 2. The equation being solved (sharp-interface form)

For an insoluble surfactant of surface concentration $\Gamma$ (mass per unit *area*) on an interface moving
with velocity $\mathbf{u}$, the transport equation is (Stone 1990):

$$
\underbrace{\frac{\partial \Gamma}{\partial t}}_{\text{(1) rate of change}}
\;+\; \underbrace{\nabla_s\cdot\left(\Gamma\,\mathbf{u}_s\right)}_{\text{(2) surface convection}}
\;+\; \underbrace{\Gamma\,\kappa\,\left(\mathbf{u}\cdot\mathbf{n}\right)}_{\text{(3) dilatation}}
\;=\; \underbrace{D_s\,\nabla_s^2\,\Gamma}_{\text{(4) surface diffusion}}
$$

- **(2) Surface convection** — surfactant slides *along* the interface with the tangential velocity
  $\mathbf{u}_s$.
- **(3) Dilatation** — when the interface **stretches**, the same surfactant spreads over more area, so
  $\Gamma$ drops. $\kappa$ is the curvature, $\mathbf{u}\cdot\mathbf{n}$ the normal velocity.
- **(4) Surface diffusion** — $\nabla_s^2$ is the Laplace–Beltrami operator: diffusion *constrained to the
  surface*. $D_s$ is the surface diffusivity, and the surface Péclet number is
  $\mathrm{Pe} = \dot\gamma R^2 / D_s$.

Every operator carrying a subscript $s$ lives *on the surface*. That is the whole problem, as the next
section explains.

## 3. Why that equation cannot be coded literally in MFC

MFC is a **diffuse-interface** solver. There is no surface mesh, no list of surface elements, and neither
$\mathbf{n}$ nor $\kappa$ exists as primary data. The interface is represented by a **color function** $c$
that varies smoothly from $0$ to $1$ over a few cells. There is no surface on which to evaluate
$\nabla_s\cdot$.

So terms (2), (3) and (4) cannot be evaluated directly. Something else is needed.

## 4. The reformulation: transport a volume density instead

The trick (Teigen 2009/2011; James & Lowengrub 2004; Jain 2024) is to stop tracking $\Gamma$ and instead
track

$$\tilde\Gamma \;=\; \Gamma\,\lvert\nabla c\rvert$$

that is, surfactant per unit **volume**, smeared across the interface band.

Here $\lvert\nabla c\rvert$ acts as a **surface delta function**: it is $\approx 0$ away from the interface
and peaks on it. It converts a surface integral into a volume integral,

$$\int_S \Gamma \, \mathrm{d}S \;=\; \int_V \Gamma\,\lvert\nabla c\rvert \, \mathrm{d}V \;=\; \int_V \tilde\Gamma \, \mathrm{d}V$$

**The key result**, and the reason the whole scheme works:

> $\tilde\Gamma$ obeys a *plain conservative advection equation in the volume*:
> $$\frac{\partial \tilde\Gamma}{\partial t} \;+\; \nabla\cdot\left(\tilde\Gamma\,\mathbf{u}\right) \;=\; 0$$

Terms (2) **and** (3) come out automatically. There is **no curvature term and no $\nabla_s\cdot\mathbf{u}_s$
term anywhere in the code.** This works because $\lvert\nabla c\rvert$ itself changes as the interface
stretches: stretch the interface and $\lvert\nabla c\rvert$ falls in exactly the way needed to make
$\Gamma = \tilde\Gamma/\lvert\nabla c\rvert$ dilute correctly. The stretching is captured *structurally*, not
by a term someone wrote down.

The physical concentration is recovered by dividing back,

$$\Gamma \;=\; \frac{\tilde\Gamma}{\lvert\nabla c\rvert} \qquad \text{on the band where } \lvert\nabla c\rvert > 10^{-6}$$

(the threshold is `capillary_cutoff`), and this recovery is used only where $\Gamma$ is actually *needed* —
i.e. inside the $\sigma(\Gamma)$ closure.

## 5. How it is discretized

**Advection.** $\tilde\Gamma$ is carried as one extra **conserved scalar** (`eqn_idx%surf`). It is not
density-weighted, so its conservative and primitive forms are identical. It is reconstructed with **WENO5**
and fluxed by the same **HLLC** Riemann solver as the rest of the state, advected at the **contact-wave
speed** $s_S$:

$$F_{\mathrm{surf}} \;=\; \tilde\Gamma^{\mathrm{up}}\, s_S \qquad\qquad
\frac{\mathrm{d}\tilde\Gamma_j}{\mathrm{d}t} \;=\; \frac{F_{j-1} - F_j}{\Delta x}$$

The update is a **pure flux difference** — a discrete divergence. This has a consequence that must be stated
plainly:

> **Because the update is a flux difference, the total $\sum \tilde\Gamma$ can change only through the domain
> boundaries. Conservation to machine precision is therefore guaranteed by construction. It confirms the
> implementation is correct; it is NOT evidence that the physics is right.** Any conservative scheme, right
> or wrong, does this. Do not present it as validation.

The color function $c$ is advected like a volume fraction (with the usual non-conservative correction).

**Surface diffusion.** This was the hardest piece, and the one place a real bug was found and fixed:

| attempt | form | outcome |
|---|---|---|
| naive projection | $D_s \nabla\cdot\left[(\mathbf{I} - \mathbf{n}\otimes\mathbf{n})\nabla\tilde\Gamma\right]$ | **leaked** $\tilde\Gamma$ off curved interfaces — band width grew $0.044 \to 0.168$ |
| Rätz–Voigt | $D_s \nabla\cdot\left[\lvert\nabla c\rvert(\mathbf{I} - \mathbf{n}\otimes\mathbf{n})\nabla\left(\tilde\Gamma/\lvert\nabla c\rvert\right)\right]$ | leak gone, but the $\tilde\Gamma/\lvert\nabla c\rvert$ division is fragile → rate only $0.58\times$ exact |
| **Jain 2024, Eq. (6)** *(implemented)* | $D_s\left(\nabla\tilde\Gamma - 2\left(\tfrac{1}{2} - c\right)\dfrac{\tilde\Gamma}{\varepsilon}\,\mathbf{n}\right)$, $\;\varepsilon \approx \Delta x$ | **correct** — $0.986\times$ exact |

The implemented flux has two parts:

$$\mathbf{F}_s \;=\; \underbrace{D_s \nabla\tilde\Gamma}_{\text{isotropic diffusion}} \;-\; \underbrace{2 D_s\left(\tfrac{1}{2} - c\right)\frac{\tilde\Gamma}{\varepsilon}\,\mathbf{n}}_{\text{sharpening flux}}$$

Note that $\left(\tfrac{1}{2} - c\right)$ is $+\tfrac{1}{2}$ outside the drop and $-\tfrac{1}{2}$ inside, so
the second term pushes $\tilde\Gamma$ *back onto the $c = \tfrac{1}{2}$ contour from both sides*.

Isotropic diffusion alone would smear surfactant off the interface; the sharpening flux re-confines it. Their
sum reproduces true Laplace–Beltrami surface diffusion **with no projection tensor and no division by
$\lvert\nabla c\rvert$** — which is precisely why it is robust where the first two attempts were not. It is
written in divergence form, so it also conserves $\tilde\Gamma$.

**The $\sigma(\Gamma)$ closure.**

| `sigma_model` | closure |
|---|---|
| 0 | $\sigma = \text{const}$ |
| 1 | $\sigma + \dfrac{\mathrm{d}\sigma}{\mathrm{d}T}\left(T - T_{\mathrm{ref}}\right)$ — thermal Marangoni |
| 2 | $\sigma + \dfrac{\mathrm{d}\sigma}{\mathrm{d}\Gamma}\,\Gamma$ — linear solutocapillary |
| 3 | $\sigma_0\left[1 + E\,\ln\!\left(1 - \Gamma/\Gamma_\infty\right)\right]$ — **nonlinear Langmuir** (`sigma_El` $= E$, `surf_max` $= \Gamma_\infty$) |

$\sigma$ is floored at $10^{-3}\sigma_0$ so that Langmuir saturation ($\Gamma \to \Gamma_\infty$ drives
$\sigma \to 0$) cannot crash the capillary force.

**Where the Marangoni force comes from — this is the elegant part.** There is **no Marangoni source term in
the code.** The capillary stress tensor is

$$\mathbf{T} \;=\; \sigma\,\left(\mathbf{I} - \mathbf{n}\otimes\mathbf{n}\right)\lvert\nabla c\rvert$$

(the solver stores its flux form $\boldsymbol{\Omega} = -\mathbf{T}$; Schmidmayer et al. 2017). Taking
$\nabla\cdot\mathbf{T}$ with a *spatially varying* $\sigma$ produces **both** the normal (Laplace) force
**and** the tangential (Marangoni) force automatically. Making $\sigma$ a function of $\Gamma$ is therefore
the *entire* coupling: nothing else had to be added.

---

## 6. What each example establishes — tiered by strength of evidence

The evidence is **not all of the same quality.** Presenting it as though it were is the single easiest way to
lose credibility. Use these three tiers.

### Tier 1 — validated against an exact analytic solution *(the strongest evidence)*

**`1D_/2D_/3D_solutocapillary_diffusion`** — Laplace–Beltrami eigenmode decay.

A surfactant pattern that is an eigenmode of the surface Laplacian decays at a rate known in closed form:

| geometry | pattern | exact decay rate | measured |
|---|---|---|---|
| flat | $\cos(kx)$ | $D_s k^2$ | $-0.1\%$ |
| circle | $\cos(m\theta)$ | $m^2 D_s / R^2$ | $0.986\times$ |
| sphere | $Y_{lm}$ | $l(l+1)\,D_s / R^2$ | $0.949\times$ |

**This is the only quantitative validation in the set, and it is the one to lead with.** It proves the
surface-diffusion operator is *the operator we claim it is* — and it is what caught the original leaking
implementation.

### Tier 2 — validated against published benchmarks, qualitatively

These reproduce the **directions and orderings** the benchmark papers establish — which is what those papers'
acceptance gates ask for — **not their numbers.** Section 7 explains why not.

| example | benchmark | result |
|---|---|---|
| **M1** `2D_Xu2006_surfactant_shear` | Xu et al. 2006, Table 3.1 | $D$ rises with coverage ($0.563 \to 0.566 \to 0.572$); surfactant sweeps to the tips; $\sigma$ minimal there |
| **M2** `2D_Xu2012_surfactant_sweep` | Xu et al. 2012, Table 3.2 | $\mathrm{Ca}\uparrow \Rightarrow D\uparrow$ ($0.455/0.608/0.706$); $\lambda\uparrow \Rightarrow D\downarrow$ ($0.612/0.608/0.582$); $\mathrm{Re}\uparrow \Rightarrow D\uparrow$ ($0.608/0.652$); $\mathrm{Pe}\uparrow \Rightarrow$ more tip-concentrated ($2.91/3.38/3.46$) |
| **M3** `2D_PimentaOliveira_rheology` | Pimenta & Oliveira 2021, Table 3.3 | $[\eta_m] = 0$ exactly for a clean drop; grows with coverage ($0 \to 0.018 \to 0.054$); $N_1 > 0$; $[\eta] = [\eta_c] + [\eta_m]$ |

For M3 the bulk interfacial stress is the volume average of the capillary stress tensor,

$$\Sigma_{ij} \;=\; \frac{1}{V}\int_V \sigma(\Gamma)\left(\delta_{ij} - n_i n_j\right)\lvert\nabla c\rvert\,\mathrm{d}V,
\qquad [\eta] = \frac{\Sigma_{xy}}{\mu\,\dot\gamma\,\phi},
\qquad N_1 = \Sigma_{xx} - \Sigma_{yy},
\qquad \phi = \frac{\pi R^2}{WH}$$

Splitting $\sigma$ into its band-mean and fluctuation splits $[\eta]$ into capillary and Marangoni parts
exactly. The sharpest single result in this tier is **$[\eta_m] = 0$ for a clean drop**: with no surfactant
there are no tension gradients, so the Marangoni stress must vanish identically — and it does, with no
tuning.

### Tier 3 — consistency and characterization, **no reference solution**

**`2D_solutocapillary_transport` (M0)** — this is *not* a validation, and its README says so explicitly.
It provides:

- mass conserved to $1.00000$ (self-consistency — see the warning in §5),
- numerical diffusion *quantified* ($+11\%$ band broadening per lap; $10.3\%$ $m=1$ decay per rotation at
  $R/\Delta x \approx 26$),
- transport-follows-flow agreement to ${\sim}1\%$ (a real cross-check: momentum and surfactant are
  independently evolved fields, and $\Gamma$ moves at exactly the rate the velocity field says it should),
- a documented **negative result**: solid-body rotation fails as an exact reference (the *flow* decays
  $13.5\%$ per rotation, boundary-driven; the surfactant rides it faithfully). See that README for why this
  is fundamental rather than a tuning problem.

**Not validation at all:** `2D_/3D_solutocapillary_marangoni` are qualitative demos.

---

## 7. The limits — state these before you are asked

All four are consequences of MFC being an **explicit compressible** solver. None is a bug.

1. **We cannot reach the Stokes limit.** $\mathrm{Re}\to 0$ requires kinematic viscosity
   $\nu = \mu/\rho \to \infty$, and the explicit viscous time step $\Delta t \sim \Delta x^2/\nu$ then
   vanishes. Every drop benchmark (M1/M2/M3) runs at **finite $\mathrm{Re}\approx 1$**, where the papers use
   $\mathrm{Re}\approx 0$. This is the main reason we match trends, not numbers.
2. **We cannot script a velocity field.** MFC always solves the full momentum equation, so the classic
   imposed-velocity transport benchmarks (Xu & Zhao 2003; Jain 2024 §6.2) are a *different problem*, not a
   harder one. This is why M0 has no reference.
3. **2D diffuse-interface drops over-deform**, and Langmuir saturation destabilizes under-resolved sharp
   tips — so $\mathrm{Ca}$ is capped near $0.3$ rather than Xu's $0.7$.
4. **Low density blows the acoustic CFL** ($c \sim 1/\sqrt{\rho}$), which bounds how $\mathrm{Re}$ can be
   varied.

Non-dimensional groups used throughout:

$$\mathrm{Ca} = \frac{\mu\dot\gamma R}{\sigma_0}, \qquad
\mathrm{Re} = \frac{\rho\dot\gamma R^2}{\mu}, \qquad
\mathrm{Pe} = \frac{\dot\gamma R^2}{D_s}, \qquad
\lambda = \frac{\mu_{\mathrm{drop}}}{\mu_{\mathrm{matrix}}}$$

---

## 8. How to say it — honest phrasings

The results are good. They get *undermined* by overclaiming, so use the right-hand column.

| result | ❌ do not say | ✅ say |
|---|---|---|
| mass $= 1.00000$ | "conservation validates the transport" | "confirms we implemented the conservative form correctly — it's guaranteed by construction, so it's a code check, not physics" |
| diffusion rate $0.986\times$ | "the model is validated" | "**the surface-diffusion operator reproduces the exact Laplace–Beltrami rate to $1.4\%$** — this is our quantitative validation, and it's what caught the original leaking operator" |
| $D$ rises with coverage (M1) | "matches Xu 2006" | "reproduces Xu's *trend*; we're at finite $\mathrm{Re}$, so we match orderings, not their deformation values" |
| M2 sweeps | "validated against Xu 2012" | "all four of Xu 2012's qualitative response gates come out in the right direction" |
| $[\eta_m] = 0$ clean (M3) | "matches Pimenta" | "the Marangoni stress vanishes identically without surfactant and grows with coverage — the decomposition behaves structurally as Pimenta & Oliveira require" |
| M0 transport | "validates transport" | "consistency checks and a numerical-diffusion characterization; **there is no reference solution — MFC can't impose a velocity field**" |
| overall | "the surfactant model is validated" | "the surface-diffusion operator is validated quantitatively against theory; the flow coupling is validated qualitatively against three published benchmarks; the transport is self-consistent but has no reference" |

**A finding worth volunteering, not hiding:** the guide's shorthand pairs "$\mathrm{Pe}\uparrow$" with "more
uniform surfactant." That is backwards — higher $\mathrm{Pe}$ means *weaker* surface diffusion and therefore
a *more* tip-concentrated interface. We measured it ($2.91 \to 3.38 \to 3.46$ as $\mathrm{Pe}$ goes
$1 \to 10 \to 100$) rather than asserting a direction. Volunteering this demonstrates the sweeps were
actually run.

## 9. The paragraph for a committee

> Surfactant is transported in MFC as a conserved volumetric density $\tilde\Gamma = \Gamma\lvert\nabla c\rvert$,
> following Teigen and Jain. Advecting $\tilde\Gamma$ conservatively reproduces both surface convection and
> interface stretching without any explicit curvature term, and the physical concentration is recovered as
> $\Gamma = \tilde\Gamma/\lvert\nabla c\rvert$ on the interface band. Surface diffusion uses Jain's (2024)
> sharpening-flux formulation, which we adopted after finding that the naive projected operator leaks
> surfactant off curved interfaces; the corrected operator reproduces the exact Laplace–Beltrami decay rate to
> within $1.4\%$ in 2D and $5\%$ on a sphere. Coupling to the flow requires no source term: making $\sigma$
> depend on $\Gamma$ inside the capillary stress tensor produces the Marangoni force automatically. The
> coupled model reproduces the qualitative gates of three published benchmarks — Xu et al. (2006, 2012) and
> Pimenta & Oliveira (2021) — including the requirement that the Marangoni contribution to the bulk stress
> vanish identically for a clean drop. Because MFC is an explicit compressible solver, these comparisons are
> at finite Reynolds number rather than in the Stokes limit, so they establish trends and orderings rather
> than the benchmarks' quantitative values.

## 10. References

**The full, CrossRef-verified reference list lives in the
[1D README](1D_solutocapillary_diffusion/README.md#references)** — cite from there, not from here. The key
attributions:

- **Transport equation** — Stone (1990), *Phys. Fluids A* **2**(1), 111–112.
- **Diffuse-interface surfactant** ($\tilde\Gamma = \Gamma\lvert\nabla c\rvert$, conserve-then-recover) —
  Teigen et al. (2009), *Commun. Math. Sci.* **7**(4); Teigen et al. (2011), *JCP* **230**(2), 375–393
  (closest precedent); James & Lowengrub (2004), *JCP* **201**(2), 685–722. *Note: our exact recipe is a
  consolidation of these two, not lifted verbatim — say so if asked about provenance.*
- **Surface diffusion / sharpening flux** — Jain (2024), *JCP* **515**, 113277, Eq. (6); Rätz & Voigt (2006),
  *Commun. Math. Sci.* **4**(3).
- **CSF force** — Brackbill, Kothe & Zemach (1992), *JCP* **100**(2), 335–354; the capillary-stress-tensor
  implementation in MFC follows Schmidmayer et al. (2017).
- **$\sigma(\Gamma)$ Langmuir EOS** — Pawar & Stebe (1996), *Phys. Fluids* **8**(7), 1738–1751.
- **Benchmarks** — Xu, Li, Lowengrub & Zhao (2006), *JCP* **212**(2), 590–616 (M1); Xu, Yang & Lowengrub
  (2012) (M2); Pimenta & Oliveira (2021) (M3); Stone & Leal (1990), *JFM* **220**, 161–186.
- **Eigenmode-decay validation practice** — Dziuk & Elliott (2013), *Acta Numerica* **22**, 289–396;
  Macdonald, Brandman & Ruuth (2011), *JCP* **230**(22).
