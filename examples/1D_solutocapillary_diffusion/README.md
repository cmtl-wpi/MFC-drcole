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

## Dispersion — the full eigenvalue spectrum

The canonical surface-diffusion benchmark is not a single mode but the whole spectrum: every mode `k`
must decay at its Laplace–Beltrami eigenvalue rate `D_s k²` (Xu, Li, Lowengrub & Zhao, *J. Comput.
Phys.* 2006; the surface-FEM literature, Dziuk & Elliott). Because the operator is linear, seeding a
superposition of modes and letting them decay independently recovers the entire `rate(k)` curve in one
run:

```
./mfc.sh run sweep.py -n 1 -t pre_process simulation
python3 measure.py sweep
```

![dispersion](figures/dispersion.png)

| mode | `k` | measured rate | exact `D_s k²` | error |
|---|---|---|---|---|
| `n=1` | 3.14 | 1.971 | 1.974 | −0.14 % |
| `n=2` | 6.28 | 7.851 | 7.896 | −0.57 % |
| `n=3` | 9.42 | 17.538 | 17.765 | −1.28 % |

The measured rates lie on `D_s k²` across all three modes; the error grows mildly with `k` because the
higher modes have fewer cells per wavelength (`48/n`), the expected spatial-resolution trend. This is
the operator matching the canonical mode-decay benchmark, not just one eigenvalue.

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

## Note on curved / non-aligned interfaces

Because this flat interface is grid-aligned, the projection `(I − n⊗n)` contributes nothing here, so
the 1D test alone cannot tell whether a curved-interface error is finite band thickness (converges) or
a projection bug. The **[3D sphere convergence study](../3D_solutocapillary_diffusion)** speaks to that:
the `l = 1` mode should decay at the sphere eigenvalue `l(l+1)D_s/R² = 2D_s/R²` (note: the *sphere*
value, not the flat/circle `D_s k²`), and the measured rate rises toward it as the interface is
resolved (`rate/exact ≈ 0.53 → 0.79 → 0.88` for `R/Δx = 5.3 → 10.7 → 16`), with surfactant conserved to
round-off. That is consistent with a resolution effect rather than a bug — **but** those numbers come
from a whole-field moment that is only an approximate estimator (a 2D-circle cross-check reads the same
moment as either 0.6× or 1.8× exact depending on masking, bracketing the true value), so the *rate* of
convergence is bracketed, not pinned. Read the caveat in the 3D README.

Still open, and needed to make the curved-interface convergence quantitative: a proper interfacial
measurement (recover `Γ = Γ̃/|∇c|` on the band and project onto the mode); a **tilted flat interface**
(tangent at 45°, exact `D_s k²`) to isolate the projection from curvature; and the coupled
surfactant-laden-drop benchmarks (e.g. Stone & Leal), which also need a resolved interface plus an
imposed strain — all future work.
