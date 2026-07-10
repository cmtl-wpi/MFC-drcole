# 1D solutocapillary surface diffusion — mode-decay validation

Validates MFC's tangential (interfacial) surfactant diffusion operator against an exact analytic
solution: the decay of a surface-concentration mode under pure surface diffusion.

**Scope (what this does and does not test).** The interface here is flat and grid-aligned, so the
surface diffusion is effectively 1D (along `x`) and the surface Laplacian reduces to `∂²/∂x²`. This
validates the operator's flux/divergence **scaling** (that the measured rate equals `D_s k²`), and it
is exact to −0.1% below. It does **not** exercise the tangential projection `(I − n⊗n)`: with the
normal `n = (0,1)` along a grid axis, the projection's cross-derivative terms are multiplied by
`n_x = 0` and contribute nothing. Validating the projection needs an interface whose tangent is not
grid-aligned (a tilted flat interface, or a resolved curved interface — see the note at the end).

## Problem

An insoluble surfactant on an interface is transported along it and diffuses along it with surface
diffusivity `D_s`. On a **flat** interface the surface-transport equation reduces, for a passive
surfactant on a static interface, to 1-D diffusion along the interface,

```
∂Γ/∂t = D_s ∂²Γ/∂s²    (s = arclength along the interface)
```

so a Fourier mode `Γ = Γ₀(1 + ε cos k x)` decays exactly as

```
Γ(x,t) = Γ₀(1 + ε e^{−D_s k² t} cos k x)   ⇒   amplitude A(t) = A(0) e^{−D_s k² t}.
```

The measured decay rate must equal `D_s k²`.

## Setup (`case.py`)

- Flat, grid-aligned interface at `y = 0` between two identical fluids (`48²`, periodic in `x`).
- Passive surfactant (`sigma_dGamma` unset ⇒ `σ` constant, no Marangoni, interface stays flat and
  static), seeded with the `m = 1` mode `Γ = Γ₀(1 + ε cos k x)`, `k = 2π/Lx`, `ε = 0.3`.
- `surf_diff = 0.2`. Time step is the min of the acoustic and explicit surface-diffusion CFL limits.

The surfactant is stored as the smeared area-density `Γ̃ = Γ·|∇c|`, so it is concentrated on the
interface; the diffusion operator's tangential projection `(I − n⊗n)∇Γ̃` confines diffusion to the
interface, and the total `∫Γ̃` is conserved.

## Result

```
./mfc.sh run case.py -n 1 -t pre_process simulation
python3 measure.py
```

![mode decay](figures/decay.png)

| quantity | value |
|---|---|
| exact rate `D_s k²` | 1.9739 |
| MFC measured rate | 1.9711 |
| relative error | **−0.1 %** |
| total surfactant drift | 0.000 % |

The amplitude tracks the exact exponential to within 0.1 % on this grid, and the total surfactant is
conserved to round-off — the surface-diffusion operator reproduces the analytic Laplace–Beltrami rate.

## Interface field

The surfactant density `Γ̃` on the interface, rendered with the built-in viewer:

```
./mfc.sh run case.py -n 1 -t pre_process simulation post_process
./mfc.sh viz . --var surfactant --step 0 --png --vmin 0 --vmax 16 --cmap magma
./mfc.sh viz . --var surfactant --step 1600 --png --vmin 0 --vmax 16 --cmap magma
```

| initial (`t = 0`) | final (`t ≈ 2τ`) |
|---|---|
| ![initial](figures/surfactant_field_initial.png) | ![final](figures/surfactant_field_final.png) |

The surfactant stays concentrated on the interface at `y = 0` (dark bulk everywhere else — the
tangential projection leaks no surfactant off the interface), while the along-interface `cos kx`
modulation — bright at the crest, dim at the troughs at `t = 0` — homogenizes into a uniform band as
it diffuses. Same colour scale in both frames.

## Note on curved / non-aligned interfaces (not yet validated)

On a curved interface with the tangent *not* grid-aligned, the projection cross-terms are active and
the finite band thickness matters. A `Nx=64`, `R=0.5` drop (`R/dx ≈ 10.6`, `w/R ≈ 0.2`) underestimates
`D_s/R²` by a factor of ~2.8. **This 1D test does not isolate that error**: because the flat interface
is grid-aligned, the projection contributes nothing here, so the sphere discrepancy could be finite
band thickness (a resolution effect that converges), a projection cross-term error, or both — this
example cannot distinguish them.

Two checks close that gap and are not yet done:
- **Tilted flat interface** (tangent at 45°, exact `D_s k²`) — activates the projection cross-terms
  with no curvature confound, so it validates the projection cleanly.
- **Sphere spherical-harmonic mode decay** (`Γ ~ Y_lm`, rate `l(l+1)D_s/R²`) at increasing `R/dx` —
  the canonical surface-diffusion benchmark; a convergence study separates resolution from projection.

Quantitative surfactant-laden-drop benchmarks (e.g. Stone & Leal extensional flow) additionally
require a resolved interface plus an imposed strain, and are future work.
