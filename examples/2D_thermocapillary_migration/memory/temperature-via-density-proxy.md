---
name: temperature-via-density-proxy
description: thermal_scalar was REMOVED; temperature is EOS-derived and imposed via the (per-fluid) density proxy + conduction. It is a constructed IC, not an independently solved field — show its limitations, don't oversell.
metadata:
  type: decision
---

The `thermal_scalar` feature (independent advected `T_s`) was **removed** (2026-06-22). MFC is
compressible: temperature is EOS-derived, `T = (p+p_inf)/((gam-1)*rho*cv)`
([[mfc_temperature_is_derived]]), with no separate temperature unknown. The imposed linear gradient
is set up through the **density proxy** and the closures (σ(T), μ(T), conduction) read the EOS
temperature via `f_compute_mixture_temperature`. Post-process flag `T_wrt` writes that field as
`temperature` (replaced `T_s_wrt`).

**Key result — distinct fluids work without a decoupled scalar.** Use a PER-FLUID density proxy:
`rho_i(y) = (p+p_inf_i)/((gam-1)*cv_i*T(y))`. The density RATIO is height-independent (`T(y)`
cancels), so the prescribed fluid ratio is preserved by tuning each fluid's stiffening `p_inf_i`; and
the mixture EOS recovers `T_mix(y)=T(y)` exactly across the interface (volume fractions cancel) — no
jump. Verified on `case_Ma_20.py` (Nas-Tryggvason ratio 0.5): IC T linear to machine precision,
drop/bulk density ratio 0.51, u=0 at t=0. This is why all 5 cases (incl. the distinct-fluid Fig 7 /
Fig 8/13) were kept rather than deleted. One analytic patch only (a 2nd patch's analytic density
leaks globally — [[mfc-analytic-patch-density-global]]); `alpha`/`alpha_rho`/`pres`/`cf` all analytic
off one smooth `eta(x,y)`; `-ffree-line-length-none` lets the long expressions compile.

**Honest stance (the user's, carried over from the old `thermal-scalar-not-physical-temperature`
note).** The imposed temperature is a *constructed IC*, not an energy-equation solution independent
of the setup. Don't oversell it: (1) without conduction a density proxy is FROZEN and advects with
the drop ([[frozen-t-proxy-advects-not-frozen]]) — only a `t/t_r<=2` anchor; conduction + isothermal
walls sustain it for finite Ma. (2) Each fluid's absolute density stratifies `~1/T` — a
compressibility artifact absent in Samareh's incompressible reference (`~dT/T`, shrunk by a large
`T0`/`T_base`). So equal-density low-Ma cases match quantitatively (Fig 5 → 0.80, Fig 6 → 0.95) but
distinct-fluid finite-Ma cases are qualitative (right overshoot-then-settle shape, compressible-fidelity
limited peak). Treat Claude-authored code/docs skeptically ([[treat-codebase-skeptically]]).
