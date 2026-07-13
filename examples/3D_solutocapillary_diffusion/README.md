# 3D surfactant surface diffusion on a sphere — a curved-interface stress test

This case tries to check that surfactant spreads along a **curved** droplet surface at the physically
correct rate. **It does not pass:** the surface-diffusion operator leaks surfactant off a curved
interface. The defect is diagnosed cleanly in the [2D-circle companion](../2D_solutocapillary_diffusion)
(with a controlled on/off experiment); the numbers here are kept to show *why* the earlier
"convergence" reading was a measurement artifact.

## In plain terms

Surfactant on a droplet's surface can spread along it — *surface diffusion*. Coat a sphere unevenly and
let it spread; theory gives the exact rate the pattern fades, set by the radius `R` and diffusivity
`D_s`. On a **flat**, grid-aligned interface MFC gets this exactly right (the
[1D case](../1D_solutocapillary_diffusion), −0.1%). On a **curved** interface it does not.

## What actually happens

A whole-field moment `M₁ = Σ Γ̃ z`, fit to an exponential, gives a rate that *rises* toward the exact
value as the mesh is refined:

| cells across the radius (`R/Δx`) | moment rate ÷ exact |
|---|---|
| 5 | 0.53 |
| 11 | 0.79 |
| 16 | 0.88 |

At first glance that looks like convergence. **It is not.** The
[2D-circle companion](../2D_solutocapillary_diffusion) resolves the interface far more finely and runs
the surface diffusion on-vs-off. It shows that on a curved interface the surfactant **diffuses off the
interface into the bulk** — the operator's tangential projection fails to confine diffusion to the
surface. Because the concentration is being drained as well as spread, *no* single decay rate exists:
on the circle the very same field reads anywhere from 0.55× to 3.1× exact depending on how you weight
it, and several of those readings move *away* from exact as the mesh is refined. The `Σ Γ̃ z` moment
here happens to rise, but it is one biased reading of a leaking field — not the operator converging.

What *is* real: the total surfactant is conserved to round-off at every resolution.

## Where this leaves the operator

- **Flat, grid-aligned interfaces:** exact (1D, −0.1% across the spectrum).
- **Curved interfaces:** the operator leaks surfactant normally and does **not** reproduce the exact
  rate. It needs a curvature-aware formulation (see the 2D companion) before curved-interface surface
  diffusion can be trusted.
- **Unaffected:** `surf_diff` defaults to `0`. The core surfactant advection and the σ(Γ) Marangoni
  coupling do not use surface diffusion and are not touched by this defect.

## Details (the math and how to run it)

- **Setup:** a passive surfactant (`sigma_model = 0`, constant σ, nothing flows — the drop stays put),
  seeded with the `l = 1` spherical-harmonic mode `Γ = Γ₀(1 + ε z/R)` (`z/R = cos θ`).
- **Exact rate:** a spherical-harmonic mode of degree `l` decays as `exp(−l(l+1) D_s/R² · t)`. For
  `l = 1` that is `2 D_s/R²` — the *sphere* eigenvalue, **not** the flat/circle value `D_s k²`.
- **Measured** as the `z`-moment `M₁ = Σ Γ̃ z` — a biased estimator of a leaking field, kept only to
  show the artifact.
- **Run the sweep** — one build serves every resolution:
  ```
  examples/3D_solutocapillary_diffusion/run_convergence.sh   # R/Δx ≈ 5, 11, 16 (release build, multi-rank)
  ```

| `R/Δx` | band `w/R` ≈ | moment rate | exact `2 D_s/R²` | rate / exact | surfactant drift |
|---|---|---|---|---|---|
| 5.3 | 0.37 | 0.845 | 1.600 | 0.53 | 0.000 % |
| 10.7 | 0.19 | 1.259 | 1.600 | 0.79 | 0.000 % |
| 16.0 | 0.12 | 1.400 | 1.600 | 0.88 | 0.000 % |

The rising ratio is **not** evidence of convergence — it is one biased reading of a field that is losing
surfactant off the interface (diagnosed in the [2D companion](../2D_solutocapillary_diffusion)). Mass
conservation to round-off is the only solid result here.

## References

The prior work backing the surfactant model and the eigenmode-decay validation is listed in the
[1D README](../1D_solutocapillary_diffusion#references). Most relevant here: spherical harmonics are
eigenfunctions of the surface Laplacian with eigenvalue `l(l+1)/R²`, so the `l = 1` mode decays at
`2D_s/R²` — the sphere instance of the standard Laplace–Beltrami eigenmode-decay test (Dziuk & Elliott
2013; Macdonald, Brandman & Ruuth 2011). James & Lowengrub (2004) is the closest single precedent using
both a surface-diffusion-on-a-drop check and a coupled surfactant benchmark.
