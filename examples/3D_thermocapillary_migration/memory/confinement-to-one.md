---
name: confinement-to-one
description: "Recovering u_YGB validates variable surface tension as a 3-axis convergence (confinement/Ma/grid → 0 ⇒ ratio → 1.0); cube sweep + 1/W→0 Richardson is the headline; cost corner is large-box fine-grid"
metadata:
  type: project
---

`case_ygb.py` validates MFC's σ(T) Marangoni physics by recovering the Young-Goldstein-Block
analytic terminal velocity. "Recover u_YGB" is a **convergence claim, not a single number**:
`v_t/u_YGB → 1.0` as three deficits vanish, each its own sweep axis (driven by `run_ygb.py` /
reduced by `validate_ygb.py`):

- **confinement** — finite box drags the drop; cube geometry, drop centered, sweep `YGB_W` ∈
  {6,8,10,12}, fit ratio vs `1/W`, extrapolate `1/W→0` (unbounded). **This is the headline.**
- **finite Ma** — conduction is needed for a 3D plateau, but Ma>0 distorts the interfacial gradient;
  sweep `YGB_MA`, extrapolate Ma→0. Cheap: dt is acoustic-limited until Ma≈0.1.
- **grid** — finite dx; sweep `YGB_NX`, Richardson in dx.

Re_M = ρ·v_YGB·D/μ ≈ 0.018 is already deep Stokes, so NO Reynolds sweep. Converged corner is
cube/W10/Nx80/Ma0.5 (shared by all three reductions). The `samareh` geometry (offset 5D×7.5D box)
is the de-risking **anchor**: it should reproduce Samareh Fig 6 ≈0.95 in the confined box before the
cube sweep extrapolates past it toward 1.0. Hitting 1.0 unbounded is a *stronger* statement than
matching Samareh's confined 0.95.

**Cost corner (the real constraint):** 3D thermocapillary + conduction + σ(T) + WENO5 is expensive
and acoustic-dt-limited. A 40³ serial smoke ran ~4.5 s/step (40³ can't MPI-decompose — below the
weno5 ≥25-cell/block floor; smoke must be serial or ≥64³). Production cube runs (Nx 80-96, ~3 t_r)
are multi-hour each even at 8-27 ranks; the W12/Nx128 grid-check is the ~1-day corner. The full
sweep is a multi-day background/batch campaign — launch heavy `run_ygb.py` selectors under nohup,
not interactively. Rank rule (`ranks_for`): cube `(Nx//32)³` keeps blocks ≥32 cells. See
[[why-thermal-scalar]] and project memory [[mfc-sigmaT-3d-drift-no-conduction]].
