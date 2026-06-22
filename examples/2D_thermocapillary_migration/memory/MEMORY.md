# Example memory — 2D_thermocapillary_migration

Memory scoped to this example. One line per memory; full content in the linked file. Format
and conventions are in `../../CLAUDE.md` (Memory section). Project-wide facts go in
`../../memory/`.

- [fig ↔ tc ↔ case mapping](fig-tc-case-mapping.md) — run.py figure targets, run-dir axis, and source case files (run.py docstring is stale)
- [TC1 droop is numerical interface diffusion](tc1-droop-is-numerical-interface-diffusion.md) — late-time velocity sag = color-band smearing (×1.4, vol/shape conserved, −0.98 corr); NOT conduction/grid/deformation; sharper interface helps
- [Temperature via density proxy](temperature-via-density-proxy.md) — thermal_scalar removed; T is EOS-derived, imposed via the (per-fluid) density proxy + conduction; distinct fluids work without a scalar; it's a constructed IC — show its limits, don't oversell
