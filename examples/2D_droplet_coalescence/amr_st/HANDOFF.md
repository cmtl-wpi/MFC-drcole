# Handoff — AMR + surface tension containment experiment (upstream PR 1628)

**Goal:** test whether *interface containment* eliminates the coarse/fine seam-current
failure that made upstream PR 1628 gate `surface_tension` under AMR. If it does, droplet
cases can use the PR's 2:1 block AMR with a static (later: guarded dynamic) block, and the
result is actionable feedback on the PR.

**Hypothesis:** the seam failure requires `|∇c| > capillary_cutoff` (= 1e-6,
`src/common/m_constants.fpp:45`) within a stencil width of the block boundary. Away from
the diffuse interface the color function is exactly 0 or 1, ∇c ≡ 0, and neither level
computes any capillary force — the inconsistent-normal mechanism has nothing to act on.
Keep the interface ≥ margin inside the block and the failure mode is absent by
construction.

**Provenance:** PR https://github.com/MFlowCode/MFC/pull/1628, branch
`sbryngelson:up/mega`, head SHA `ace2285a7e72fdce9dc3fbc3a5629e1b9d1a89b7` (2026-07-08).
The PR is open and moving — **pin this SHA** for all runs and record it in every manifest.

**Read first:** `.claude/RESEARCH_WORKFLOWS.md` (run/archive/plot contracts — WF-BUILD,
WF-RUN-LOCAL, WF-PLOT, WF-ARCHIVE apply throughout).

---

## Background (verified facts, with pointers)

- **The gate:** one `@:PROHIBIT(surface_tension, ...)` inside the `if (amr)` section of
  `src/simulation/m_checker.fpp` on the PR branch. Two *other* ST gates (`active_box`,
  `hybrid_riemann`) exist in the same file — leave those alone.
- **The mechanism** (PR's `docs/documentation/case.md` §7.1 "AMR + surface tension
  (unsupported)"): the capillary stress is Ω ∝ ∇c⊗∇c/|∇c|
  (`src/simulation/include/inline_capillary.fpp`) — normalized by |∇c|, so it depends on
  the interface-normal *direction*. Conservative-linear prolongation of the color into
  fine ghost cells can't reproduce the coarse normal across the 2:1 seam → unbalanced
  force → growing spurious seam current. Conservation is NOT the problem (capillary flux
  is refluxed; defects at machine precision).
- **Their three attempted fixes all had the interface crossing/near the seam** and failed
  at ~540× / ~27× / ~556× the baseline seam velocity. The contained case
  (interface far from the seam) was apparently never characterized. That's the gap this
  experiment fills.
- **Verified at `ace2285a`:** `src/simulation/m_amr.fpp` contains no
  surface-tension/capillary special-casing (grep for `surface_tension|capillary|color`
  returns nothing) — the attempted-fix machinery was removed; with the gate lifted, the
  color function rides the generic prolong/restrict/reflux path. That generic path is
  exactly what we want to measure.
- **Our case already satisfies the other AMR requirements:** `model_eqns=2`,
  `num_fluids=2`, `mpp_lim='T'`, WENO (`recon_type=1`), RK3 (`time_stepper=3`).
- **AMR shape:** two-level, 2:1 only. Static block = `amr_regrid_int=0` +
  `amr_block_beg/end(i)` in level-0 cell indices. Parameter table:
  PR `docs/documentation/amr.md` §Parameters. Constraint to respect: the block may cover
  at most ~half of any rank's subdomain per dimension.

---

## Phase 0 — branch, gate removal, build (~1 h)

```bash
cd ~/repos/MFC
git remote add sbryngelson https://github.com/sbryngelson/MFC.git 2>/dev/null || true
git fetch sbryngelson up/mega
git worktree add ../MFC-amr-st ace2285a7e72fdce9dc3fbc3a5629e1b9d1a89b7 -b amr-st-experiment
```

1. In `../MFC-amr-st/src/simulation/m_checker.fpp`, delete **only** the
   `@:PROHIBIT(surface_tension, ...)` in the AMR section (grep `amr does not support
   surface_tension`). Do not touch the `active_box`/`hybrid_riemann` ST gates.
2. Build with the Intel toolchain (ifx ~2× faster than gfortran here). Source the full
   Intel env atop the build command; `setvars.sh` returns 3 on re-init, so under `set -e`
   use `|| true`. Plain build, **no** `--case-optimization` (dev iteration; also the
   coalescence case has analytic ICs → case-specific slug; never `--no-build` after an IC
   edit).
   ```bash
   cd ~/repos/MFC-amr-st && ./mfc.sh build -t pre_process simulation post_process -j 32
   ```
3. **Pre-flight greps** (resolve before designing cases — each is a possible plan change):
   - Is `bc = -8` gated under AMR? (`amr.md` gates `bc=-4` only; confirm in `m_checker.fpp`.)
   - Is `cfl_dt`/`cfl_adap` allowed under AMR? (AMR restart code references `cfl_dt`, so
     likely yes; if gated, use fixed `dt` matched across all runs in a trio.)
   - Is `amr_max_blocks = 1` accepted? (Slots preallocate at max size; 1 saves memory for
     a single static block.)
   - Actual `buff_size` with `surface_tension=T` (read `s_configure_coordinate_bounds` in
     `src/common/m_helper_basic.fpp`) — sets the seam-band width used in the audits below.

Keep this main repo checkout untouched; all experiment code/cases/runs live under
`examples/2D_droplet_coalescence/amr_st/` (this dir), binaries come from the worktree.
Run with absolute case paths (output lands next to `case.py`). `mpirun` is blocked in the
sandboxed shell — disable sandbox or exec the binary directly; watch `MFC.out`, not stdout.

---

## Phase 1 — static-drop Laplace seam test (~half a day)

The decisive, cheap test. Directly reproduces the PR authors' seam-velocity metric with
the interface *contained* instead of crossing the seam.

**Case** (`cases/case_laplace.py`, new): one static drop, radius R and fluids/EOS taken
from case f of `../case.py` (keep the cavitation-fix sound speeds), zero background flow,
pressure initialized with the Laplace jump (the balanced IC already used in the
coalescence cases). Uniform grid, square domain ±4R, coarse dx = R/50 (≙ D/100)
→ 400×400 cells. No stretching. Run ≥ 10 capillary periods τ = sqrt(ρ_l R³/σ).

**Three runs:**
| id | grid | AMR |
|---|---|---|
| `laplace__uniform_coarse` | 400², dx=R/50 | off |
| `laplace__uniform_fine` | 800², dx=R/100 | off |
| `laplace__amr_static` | 400² + static block | block x,y ∈ ±1.75R → cells 112–288 (176 ≤ 200 = half ✓) |

Interface at R, seam at 1.75R → 0.75R ≈ 37 coarse cells of separation ≫ any stencil.

**Metrics** (script in `analysis/`, reads `simulation.inp` — never hardcode; show raw
|u| pcolormesh at final time *before* extracting scalars):
- max|u|(t) domain-wide (parasitic current level) for all three runs
- max|u|(t) restricted to a band of ±2·buff_size cells around the block boundary (AMR run)
- pressure jump vs σ/R; drop-center drift
- containment audit: α in the seam band stays within 1e-9 of {0,1} at every saved step

**Success:** AMR max|u| within ~2× of `uniform_coarse` and **non-growing** over the run.
The §7.1 failure signature is 27–540× and growing — the two outcomes are unambiguous.
**If this fails, stop.** The hypothesis is dead; write up the measured growth rate vs
seam distance and skip to Phase 3 (that's still useful PR feedback).

---

## Phase 2 — coalescence trio (~a day)

**Case** (`cases/case_amr.py`, derived from `../case.py --case f`): keep the physics,
EOS, balanced IC, and per-case R_um exactly; replace the stretched grid with uniform to
keep stretch×AMR interaction out of the first test:
- Domain ±4D × ±3D, uniform D/100 → 800×600 (the ±4D also honors the case-p lesson on
  outflow distance). `bc = -8` as in the original (pending Phase 0 pre-flight).
- t_end ≥ 3 ms (house rule — drop is mid-stretch at 1 ms).

**Three runs:**
| id | grid | AMR |
|---|---|---|
| `coal_f__ref_d200` | 1600×1200, D/200 | off |
| `coal_f__coarse_d100` | 800×600, D/100 | off |
| `coal_f__amr_static` | 800×600 + block | x ∈ ±1.8D → cells 220–580 (360 ≤ 400 ✓), y ∈ ±1.4D → cells 160–440 (280 ≤ 300 ✓) |

Block sizing rationale: drop centers start at ±1.18D (edges ±1.68D) and advect inward;
case f (B=0.55, oblique) merged-blob excursion must stay inside ±1.4D in y — check
against the archived D/200 case-f frames before launching, widen if it gets within
~10 coarse cells of the block edge.

**Metrics:**
- bridge radius vs t: AMR should track `ref_d200` markedly better than `coarse_d100`
  (reuse the extraction in `../analysis/` per the experiment-comparison recipe)
- V_liq(t) and domain max|u|(t) sanity traces (the case-p failure signatures)
- |u| pcolormesh at merge time and late time, raw field first
- same seam-band velocity + α containment audit as Phase 1; **a containment violation
  voids the run** — rerun with a bigger block, don't hand-wave it

MPI note: ≥ ~25 cells per rank-block per split dimension, and the block ≤ ~half of each
rank's subdomain — pick the decomposition to satisfy both (or run modest rank counts;
these are cheap 2D runs on the EPYC).

---

## Phase 3 — write-up

- `results/summary.json` (metric archive: parasitic-current ratios, growth rates, bridge
  RMS error vs reference, containment-audit pass/fail) + figures.
- Every run gets a `run_manifest.json` (git SHA = `ace2285a…`, invocation, slug, ranks,
  host, wall time) — a run without a manifest is not reproducible.
- If Phase 1+2 are positive: draft a PR 1628 comment — "contained-interface case
  characterized: seam current at parasitic baseline; suggest relaxing the hard gate to a
  containment guard (the PR already ships this pattern for EL bubble clouds and moving
  IBs)". **Do not post to GitHub without Davis's approval.**
- Update project memory (status of this handoff, key numbers).

## Phase 4 — contingent, out of scope here

Dynamic-regrid guard (only if Phases 1–2 pass and the PR merges or Davis wants it):
c-based tagging dilated by the containment margin + per-stage abort guard mirroring the
EL-cloud/moving-IB containment pattern already in the PR. Sketch only; do not start
without discussion.

---

## Caveats / framing

- **The payoff is capped:** two-level 2:1 AMR ≈ 2–4× at best in 3D, and it does NOT buy
  film-scale resolution at rupture. This experiment is about feasibility + upstream
  feedback, not an immediate production speedup for the coalescence study.
- PR unmerged and large; results are against `ace2285a` and may need re-validation.
- Never regenerate golden files; never claim success while a containment audit fails;
  never post upstream without approval.

**Estimated total:** ~2 days wall, mostly 2D runs.
