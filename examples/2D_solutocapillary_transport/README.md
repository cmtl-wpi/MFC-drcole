# 2D surfactant transport under a prescribed flow (M0)

Validates that an insoluble interfacial surfactant is **transported correctly along a moving/deforming
interface** — the advection half of the model, complementary to the
[surface-diffusion validation](../2D_solutocapillary_diffusion). Two prescribed-velocity tests
(guide milestone "M0"): a translating drop (interface **confinement**) and a drop in extensional flow
(the interface-**stretching** term). Both decouple the surfactant transport from the flow solver by
imposing a known velocity, so any error is in the transport, not the hydrodynamics.

## Why prescribed velocity, and why these two flows

MFC is a full compressible solver; the only velocity field it maintains *exactly* is **uniform
translation** (Galilean invariance). An **extensional** flow `u=(εx,−εy)` is divergence-free, so it is
low-Mach-friendly and holds well enough over a short advective time to drive a clean redistribution.
Together they exercise the two things interfacial transport must get right:

- **Translation** → does the surfactant stay *confined* to the interface as it moves through the mesh,
  or does numerical diffusion smear it into the bulk?
- **Extensional strain** → when the interface stretches, does the surfactant *redistribute* correctly
  (dilute where stretched, concentrate where the surface flow sweeps it)?

## Test 1 — confinement under translation

A drop with a `cosθ` surfactant pattern translates one full lap around the periodic domain
(`surf_diff=0`, `sigma_model=0`, so pure passive advection). The surfactant should ride rigidly with the
interface.

![confinement](figures/confinement.png)

| metric | result |
|---|---|
| interfacial mass | **1.00000** (conserved to machine precision) |
| `cosθ` pattern amplitude | preserved (0.51 → 0.54) |
| band width (radial spread) | 0.058 → 0.065 — **+11% over one full lap** (8 drop-diameters of travel) |
| off-band leakage | ~2.4% → 3% (mostly the initial tanh tail) |

The ring stays a well-defined ring — MFC's WENO interface capturing keeps `Γ̃` fairly sharp (a naively
advected scalar would smear far more). There is a *mild* numerical-diffusion broadening (~11% per lap);
it is small for the short trajectories of a coalescence run, and it is what a Jain-style advection
sharpening flux would remove — at the cost of adding artificial surface diffusion (positivity requires
`Pe_c = Δx|u|/D ≤ 1`), which is undesirable in the high-Péclet regime. Left as a documented option.

```
MFC_NX=128 ./mfc.sh run examples/2D_solutocapillary_transport/case_translate.py --no-debug -n 4 -t pre_process simulation
python3 examples/2D_solutocapillary_transport/measure_translate.py examples/2D_solutocapillary_transport
```

## Test 2 — the stretching term under extensional strain

A drop with a **uniform** surfactant coating sits in `u=(εx,−εy)` (`ε=1`, `Ca=με R/σ=0.5`). The flow
elongates the drop along `x`; the surface flow and interface stretching redistribute the surfactant.

![strain](figures/strain.png)

Left: uniform coating at `t=0`. Right (`t=0.10`): the drop has elongated along `x` and the surfactant is
visibly **concentrated at the x-tips** (the elongation axis) and depleted at the y-poles — the classic
Stone & Leal trend. The `m=2` mode amplitude `a₂` rises positive (to ~0.10), and **total interfacial
surfactant is conserved to machine precision** the whole time. This is the interface-stretching term
`−Γ ∇ₛ·uₛ` plus surface convection working; the conservative `Γ̃ = Γ|∇c|` transport captures it
structurally (no explicit stretching term is coded). The imposed strain is not maintained indefinitely
(MFC relaxes it), so the deformation plateaus at `D≈0.065` — long enough to demonstrate the redistribution.

```
MFC_NX=128 ./mfc.sh run examples/2D_solutocapillary_transport/case_strain.py --no-debug -n 4 -t pre_process simulation
python3 examples/2D_solutocapillary_transport/measure_strain.py examples/2D_solutocapillary_transport
```

## What this establishes

- Interfacial surfactant **mass is conserved to machine precision** under both translation and strain.
- The surfactant **stays confined** to the interface under advection (mild, bounded numerical broadening).
- The **stretching term is correct** — surfactant redistributes to the elongation tips as theory predicts.

Together with the [surface-diffusion validation](../2D_solutocapillary_diffusion) (exact rate to ~1.4%),
this covers the transport side of the insoluble-surfactant model. The remaining validation is the coupled
surfactant-laden-drop-in-shear benchmark (Xu et al. 2006 / Stone & Leal), where the surfactant feeds
back on the flow through `σ(Γ)`.

## References

See the [1D README](../1D_solutocapillary_diffusion#references). The extensional-flow redistribution
trend is Stone & Leal (1990), *J. Fluid Mech.* 220, 161; the interface-confined transport formulation is
Jain (2024), *J. Comput. Phys.* 515, 113277.
