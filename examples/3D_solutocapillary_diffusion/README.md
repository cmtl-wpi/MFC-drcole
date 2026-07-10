# 3D surfactant surface diffusion on a sphere — convergence test

Checks that surfactant spreads along a **curved** droplet surface at the physically correct speed, and
that MFC's answer approaches the exact one as the mesh is refined.

## In plain terms

Surfactant sitting on a droplet's surface can spread out along it — *surface diffusion*. Coat a sphere
unevenly (more surfactant near the top, less near the bottom) and let it spread: the uneven pattern
smooths out toward a uniform coating. Theory gives the **exact rate** at which that pattern fades, set
only by the sphere's radius `R` and the surfactant's diffusivity `D_s`.

This case does exactly that — an uneven coating on a sphere, nothing else pushing on it — and measures
how fast MFC smooths it out, versus the exact rate.

Why a sphere, and not just the flat [1D case](../1D_solutocapillary_diffusion)? On a flat interface
lined up with the grid, the diffusion runs straight along a mesh line — the easy case. A sphere is
curved, so the diffusion has to follow a surface that cuts across the mesh at every angle. That's the
harder, realistic case, and the one that actually stresses the operator.

## The result

MFC's rate climbs toward the exact value as the interface is better resolved:

![convergence](figures/convergence.png)

| cells across the radius (`R/Δx`) | MFC rate ÷ exact |
|---|---|
| 5 | 0.53 |
| 11 | 0.79 |
| 16 | 0.88 |

On a coarse mesh MFC is too slow (about half the exact rate) because the interface is "fuzzy" — spread
over a few cells, which at low resolution is a big fraction of the whole sphere. As the mesh is refined
the fuzz shrinks and the rate rises steadily toward the exact answer, halving the error each time the
mesh doubles. The total amount of surfactant is conserved to round-off the whole time.

**The takeaway:** the operator is **correct** — it converges to the exact answer. The gap on a coarse
mesh is just a resolution effect, not a bug. The approach is gradual (first order), so getting within a
few percent needs a well-resolved interface (~30 cells across the radius) — the same resolution the
harder surfactant-laden-drop benchmarks (Stone & Leal) would need.

## Details (the math and how to run it)

- **Setup:** a passive surfactant (`sigma_model = 0`, so surface tension is constant and nothing flows —
  the drop stays put), seeded with the `l = 1` spherical-harmonic mode `Γ = Γ₀(1 + ε z/R)` (`z/R = cos θ`).
- **Exact rate:** a spherical-harmonic mode of degree `l` decays as `exp(−l(l+1) D_s/R² · t)`. For
  `l = 1` that is `2 D_s/R²` — the *sphere* eigenvalue, **not** the flat/circle value `D_s k²`.
- **Measured** as the `z`-moment `M₁ = Σ Γ̃ z` (a clean surface quantity because the stored density
  `Γ̃ = Γ·|∇c|` sits on the interface), fit to an exponential.
- **Run the sweep** — one build serves every resolution, since the `z/R` initial condition does not
  depend on the mesh:
  ```
  examples/3D_solutocapillary_diffusion/run_convergence.sh   # R/Δx ≈ 5, 11, 16 (release build, multi-rank)
  ```

| `R/Δx` | band `w/R` ≈ | measured rate | exact `2 D_s/R²` | rate / exact | surfactant drift |
|---|---|---|---|---|---|
| 5.3 | 0.37 | 0.845 | 1.600 | 0.53 | 0.000 % |
| 10.7 | 0.19 | 1.259 | 1.600 | 0.79 | 0.000 % |
| 16.0 | 0.12 | 1.400 | 1.600 | 0.88 | 0.000 % |

The steadily rising ratio (first order in the band thickness `w/R`) is exactly the behaviour expected
of a diffuse-interface surface operator, and it confirms the discretization is correct.
