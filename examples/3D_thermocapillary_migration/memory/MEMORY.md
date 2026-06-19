# Example memory — 3D_thermocapillary_migration

Memory scoped to this example. One line per memory; full content in the linked file. Format
and conventions are in `../../CLAUDE.md` (Memory section). Project-wide facts go in
`../../memory/`.

- [Why thermal_scalar, not the density proxy](why-thermal-scalar.md) — the u_YGB validation uses the decoupled scalar; case.py's density proxy advects and reverses the gradient
- [Confinement → 1.0 validation logic](confinement-to-one.md) — recovering u_YGB is a 3-axis convergence (confinement/Ma/grid → 0); the cube sweep + 1/W→0 extrapolation is the headline; cost corner
