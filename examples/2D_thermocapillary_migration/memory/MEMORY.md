# Example memory — 2D_thermocapillary_migration

Memory scoped to this example. One line per memory; full content in the linked file. Format
and conventions are in `../../CLAUDE.md` (Memory section). Project-wide facts go in
`../../memory/`.

- [fig ↔ tc ↔ case mapping](fig-tc-case-mapping.md) — run.py figure targets, run-dir axis, and source case files (run.py docstring is stale)
- [TC1 droop is numerical interface diffusion](tc1-droop-is-numerical-interface-diffusion.md) — late-time velocity sag = color-band smearing (×1.4, vol/shape conserved, −0.98 corr); NOT conduction/grid/deformation; sharper interface helps
- [Temperature via density proxy](temperature-via-density-proxy.md) — thermal_scalar removed; T is EOS-derived, imposed via the (per-fluid) density proxy + conduction; distinct fluids work without a scalar; it's a constructed IC — show its limits, don't oversell
- [case_Ma_20 NaN = 6-eq pressure-relaxation 0/0 (fixed)](thermocap-nan-pressure-relaxation.md) — bulk phase drains to alpha_rho_1=0 inside the drop → 0/0 in `s_equilibrate_pressure` (t_step 14870); NOT a GPU bug / wall corner; fixed by `max(rho_K_s,sgm_eps)`; verified to completion (45999). GPU≡CPU peak U*/U_r=0.134 @ step 11500 still holds.
- [Working mode: understanding over results](working-mode-understanding-over-results.md) — user is exploring to understand the physics/numerics here, not to ship a figure; explain mechanisms, treat "messy" runs as experiments, don't rush to deliverables
