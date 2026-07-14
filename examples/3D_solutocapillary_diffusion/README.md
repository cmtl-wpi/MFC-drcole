# 3D surfactant surface diffusion on a sphere — curved-interface check

The 3D analog of the [2D-circle](../2D_solutocapillary_diffusion) validation: does surfactant spread
along a **curved** droplet surface at the physically correct rate? This case first exposed a defect in
the original operator (it leaked surfactant off curved interfaces) and now runs on the **fixed**
operator — Jain's (2024) interface-confined scalar flux, validated to ~1.4% on the circle.

## In plain terms

Surfactant on a droplet's surface can spread along it — *surface diffusion*. Coat a sphere unevenly and
let it spread; theory gives the exact rate the pattern fades, set by the radius `R` and diffusivity
`D_s`. On a **flat**, grid-aligned interface MFC always got this right (the
[1D case](../1D_solutocapillary_diffusion), −0.1%). On a **curved** interface the *original* operator
did not — it leaked the surfactant off the surface — which the fix corrects.

## The defect and the fix

The clean diagnosis is in the [2D-circle companion](../2D_solutocapillary_diffusion): the original
projected operator `D_s·(I−n⊗n)∇Γ̃` diffused surfactant *normally*, off the interface into the bulk. A
tell-tale of that leaking field was that different ways of measuring the decay rate disagreed wildly
(0.55×–3.1× exact) and drifted *away* from exact with resolution. The old 3D numbers below came from a
`Σ Γ̃ z` moment of that leaking field, so their apparent "rise toward exact" was a measurement artifact,
not real convergence.

The fix replaces the operator with Jain's interface-confined scalar flux (isotropic diffusion + a
sharpening flux that re-confines `Γ̃` to the interface; *J. Comput. Phys.* 515 (2024) 113277, Eq. 6).
On the circle it recovers the exact rate to ~1.4% with a clean single-exponential decay and no leakage.
Mass is conserved to round-off throughout (true both before and after the fix).

## Old, superseded numbers (original *leaking* operator)

Kept only to show the artifact the 2D companion diagnosed — **not** a validation:

| `R/Δx` | moment rate | exact `2 D_s/R²` | rate / exact (leaking op) |
|---|---|---|---|
| 5.3 | 0.845 | 1.600 | 0.53 |
| 10.7 | 1.259 | 1.600 | 0.79 |
| 16.0 | 1.400 | 1.600 | 0.88 |

## Details (the math and how to run it)

- **Setup:** a passive surfactant (`sigma_model = 0`, constant σ, nothing flows — the drop stays put),
  seeded with the `l = 1` spherical-harmonic mode `Γ = Γ₀(1 + ε z/R)` (`z/R = cos θ`).
- **Exact rate:** a spherical-harmonic mode of degree `l` decays as `exp(−l(l+1) D_s/R² · t)`. For
  `l = 1` that is `2 D_s/R²` — the *sphere* eigenvalue, **not** the flat/circle value `D_s k²`.
- **Run the sweep** (fixed operator) — one build serves every resolution:
  ```
  examples/3D_solutocapillary_diffusion/run_convergence.sh   # R/Δx ≈ 5, 11, 16 (release build, multi-rank)
  ```

The rigorous curved-interface validation lives in the finely-resolved
[2D-circle companion](../2D_solutocapillary_diffusion) (0.986× exact); this sphere is the 3D
confirmation that the same operator carries over.

## References

The prior work backing the surfactant model and the eigenmode-decay validation is listed in the
[1D README](../1D_solutocapillary_diffusion#references). Most relevant here: spherical harmonics are
eigenfunctions of the surface Laplacian with eigenvalue `l(l+1)/R²`, so the `l = 1` mode decays at
`2D_s/R²` — the sphere instance of the standard Laplace–Beltrami eigenmode-decay test (Dziuk & Elliott
2013; Macdonald, Brandman & Ruuth 2011). James & Lowengrub (2004) is the closest single precedent using
both a surface-diffusion-on-a-drop check and a coupled surfactant benchmark.
