# AMR + surface-tension containment experiment

Tests whether keeping a diffuse interface **contained** well inside a static 2:1 AMR
block (so the color-function gradient is ~0 at the coarse/fine seam) avoids the spurious
seam current that made upstream PR 1628 gate `surface_tension` under AMR. See `HANDOFF.md`
for the full plan and `results/FINDINGS.md` for what was actually found.

## Provenance
- PR https://github.com/MFlowCode/MFC/pull/1628, branch `sbryngelson:up/mega`,
  **pinned SHA `ace2285a`**. Worktree at `~/repos/MFC-amr-st` (branch `amr-st-experiment`).
- Both AMR-vs-surface_tension gates lifted: `src/simulation/m_checker.fpp` (Fortran) and
  `toolchain/mfc/case_validator.py` (Python). Intel/ifx build.

## Layout
```
cases/      case_laplace.py     # Phase-1 static-drop Laplace test (one case, --variant/--c-l/...)
analysis/   mfc_read.py         # shared lustre reader (level-0 field, model-aware indices)
            seam_analysis.py    # per-run metrics + raw-field figures + RESULT_JSON
            compare_phase1.py   # trio -> verdict + overlay + summary.json
scripts/    build.sh run_case.sh run_durable.sh   # build / one-off run / durable screen run
runs/       <run_id>/           # gitignored run dirs (case.py, *.inp, restart_data, manifest)
results/    FINDINGS.md summary.json figures/      # tracked deliverables
```

## Reduced-config decision (why this deviates from the handoff spec)
Lifting the gate revealed that **AMR + surface_tension only runs at np=1** (a multi-rank
fine-halo bug — see FINDINGS) and there costs ~2.25 s/step, making the handoff's 10τ /
c_l=100 / R50 spec ~26 days per run. To answer the containment question feasibly, the
Phase-1 runs trim three feasibility knobs that do **not** affect the seam mechanism (which
is sound-speed- and domain-independent), while preserving the geometry that tests it
(interface at R, seam at `block_R·R`, margin ≫ seam band):

| knob | handoff | reduced | why safe |
|---|---|---|---|
| `c_l` | 100 m/s | 12 m/s | dt ∝ 1/c → 8× fewer steps; Ma_parasitic still ~0.03 (incompressible); seam mechanism is c-independent |
| domain | ±4R | ±3.2R | drop stays ≫ buff from the outflow; only trims empty gas |
| block | ±1.75R | ±1.5R | interface→seam margin 25 cells > seam band (12); still fully contained |
| duration | 10τ | 2τ | the 27–540× growing-failure signature appears within ≪1τ; 2τ shows plateau + no-growth |

The full-spec values remain the case defaults (`--c-l 100 --domain-R 4 --block-R 1.75
--n-periods 10`) for reference/reproduction.

Numerics also necessarily differ from the production `../case.py`: **WENO5** (AMR requires
WENO; case.py uses MUSCL) with the shipped `2D_laplace_pressure_jump` ST-safe WENO
settings, and a **2D Laplace jump σ/R** (case.py's 2σ/R is the 3D-sphere value and would
leave a static drop unbalanced). Validated: reconstructed jump matches σ/R to 0.68%.
