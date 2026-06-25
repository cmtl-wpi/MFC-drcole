# Example memory — 3D_thermocapillary_migration

Memory scoped to this example. One line per memory; full content in the linked file. Format
and conventions are in `../../CLAUDE.md` (Memory section). Project-wide facts go in
`../../memory/`.

- [Grid convergence to u_YGB = 1.0](confinement-to-one.md) — the real axis is GRID (W=10 → 1.009≈1.0); confinement is negligible (force-free drop, never swept); 0.95 is NOT a limit but Samareh's finite-res value, reproduced only as a matched-grid anchor point
- [Temperature via density proxy](../../2D_thermocapillary_migration/memory/temperature-via-density-proxy.md) — thermal_scalar removed; T is EOS-derived, imposed via the (per-fluid) density proxy + conduction (full note in the 2D example)
