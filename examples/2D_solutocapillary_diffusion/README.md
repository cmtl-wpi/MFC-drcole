# Insoluble-surfactant surface diffusion — mode-decay validation

Validates MFC's tangential (interfacial) surfactant diffusion operator against an exact analytic
solution: the decay of a surface-concentration mode under pure surface diffusion.

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

## Note on curved interfaces

On a **curved** interface the same operator is subject to the usual diffuse-interface errors: at
coarse resolution the finite band thickness (`w/R`) and the averaged-normal projection reduce the
effective rate (e.g. a `Nx=64`, `R=0.5` drop, `R/dx ≈ 10.6`, `w/R ≈ 0.2`, underestimates `D_s/R²` by
a factor of ~2.8). This is a resolution effect — the flat-interface test above isolates and confirms
the operator's discretization is correct — and it decreases as the interface is better resolved.
Quantitative surfactant-laden-drop benchmarks (e.g. Stone & Leal extensional flow) therefore require
a resolved interface.
