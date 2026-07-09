# AMR + surface-tension containment experiment — findings

Provenance: PR 1628 (`sbryngelson:up/mega`), pinned SHA `ace2285a`, worktree
`~/repos/MFC-amr-st`. Intel/ifx build. Host `ubuntu` (EPYC). Started 2026-07-08.

## Phase 0 — gate removal + build (COMPLETE)

The `surface_tension`-under-AMR gate exists in **two** places (the handoff named only the
first):
1. `src/simulation/m_checker.fpp:118` — the Fortran `@:PROHIBIT(surface_tension, ...)`
   inside `if (amr)`. Lifted.
2. `toolchain/mfc/case_validator.py:1430` — a Python `prohibit(surface_tension or
   hyperelasticity, ...)`. Lifted the `surface_tension` half; kept `hyperelasticity`.

Both edits keep the `active_box`/`hybrid_riemann` ST gates and the AMR `hyperelasticity`
gate intact. Build (ifx, all 3 targets) succeeds.

Pre-flight greps (resolved against `ace2285a`):
- `bc=-8` (`BC_CHAR_NR_SUB_OUTFLOW`) is NOT AMR-gated (only `bc=-4` is). OK.
- `cfl_adap_dt` works under AMR **lockstep**; `amr_subcycle=T` would require fixed dt
  (`m_checker.fpp:201`). We run lockstep (`amr_subcycle=F`, the default).
- `amr_max_blocks >= 1` accepted; but see the multi-rank finding below.
- `buff_size` with WENO5+viscous = `2*weno_polyn+2 = 6` → seam band ±12 cells.
- **AMR requires WENO** (`recon_type=1`); the production `case.py` uses MUSCL
  (`recon_type=2`). Cases here use WENO5 with the shipped Laplace example's ST-safe
  WENO settings — a deliberate, required deviation.

## Blocking discovery — AMR + surface_tension is broken at multi-rank

Lifting the gate does **not** yield a runnable-at-scale configuration. Measured on the
static-drop Laplace case (`cases/case_laplace.py --variant amr`, 400² base + static 2:1
block, cells 112–288):

| config | surface_tension | result | s/step |
|---|---|---|---|
| AMR, np=1  | **on**  | RUNS, stable (ICFL controlled) | ~2.25 |
| AMR, np=2  | **on**  | **hangs at step 0** (fine-halo deadlock) | — |
| AMR, np=16 | **on**  | **`MPI_Sendrecv` message-truncated abort** at step 0 | — |
| AMR, np=4  | off     | runs to completion (106 steps) — control | ~0.56 |
| AMR, np=16 | off     | **runs to completion (106 steps)** — control | **0.16** |

**The bug is surface-tension-specific, verified by controls.** The *identical* AMR case
(same 400²+block, same 16-rank 4×4 decomposition) runs to completion **without**
surface_tension but truncation-aborts **with** it. So it is not a general
decomposition/face-size bug — enabling the ST color function (which adds a variable to
`sys_size`) breaks the parallel fine-level halo exchange
(`s_mpi_sendrecv_amr_fine_halo`, `src/simulation/m_mpi_proxy.fpp:82`). This path was never
exercised because ST was gated under AMR.

The AMR core math is correct at every config (init diagnostics):
`prolong linear-reproduction err = 2.2e-16`, `restrict independent-integral err = 4.3e-16`,
`restrict-prolong conservation err = 2.2e-16`. AMR itself also parallelizes well (0.16
s/step at np=16 without ST) — so the ~2.25 s/step at np=1-with-ST is purely the absence of
parallelism, forced by the halo bug.

### Root cause (CONFIRMED) and fix
Not the fine-fine halo (an earlier guess; `s_mpi_sendrecv_amr_fine_halo` is dead code).
The real cause: the AMR fine-block advance runs the full solver RHS
(`s_amr_fine_stage_advance → s_compute_rhs`, `m_amr.fpp:2532`). For surface_tension that
RHS calls `s_get_capillary → s_populate_capillary_buffers` (`m_surface_tension.fpp:275`),
which does a **coarse-decomposition MPI halo** on the color-gradient `c_divs`
(`s_mpi_sendrecv_variables_buffers`). But only the block-owning rank(s) enter the fine
advance, so that MPI has no matching partner on the other ranks → **deadlock (np=2) /
truncation (np≥ with an uneven split)**. Surface_tension is the **only** physics that does
an MPI exchange *inside* `s_compute_rhs`, which is exactly why the no-ST controls run.
The IGR solver already guards its analogous populate with `if (.not. amr_in_fine_advance)`
(`m_igr.fpp:306`); surface_tension was simply missed — it was gated under AMR, so the path
was never exercised.

**Fix 1 (deadlock) — DONE & VERIFIED** (branch `fix/amr-st-multirank-halo`): skip the
c_divs MPI halo during the fine advance (the fine block is interior; its c_divs ghosts are
filled locally, as at np=1 where the halo is a no-op). Guard the populate with
`amr_in_fine_advance`, mirroring `m_igr.fpp:306`; provably a no-op for non-AMR runs.
**Verified: np=2 no longer deadlocks at step 0** — it runs all 260 steps.

**BUT a second, deeper problem remains — np=2 BLOWS UP** (not fixed). With the deadlock
gone, np=2 diverges from np=1 and blows up (alpha1→4.5, pres→2.5e5, |u|→2.9 vs np=1's
alpha1∈[0,1], pres~285, |u|~0.16). Localized to **(x/R≈0, y/R≈±1) — the interface∩rank-
boundary**. Diagnosis (data-backed):
- The coarse solve is **bit-identical** np1/np2 away from the block (diff 0). The
  divergence is entirely **inside the fine block**, at the interface.
- With `amr_max_blocks=1 < num_procs`, the block is **whole-owned by one rank**, which must
  **gather** the other rank's coarse cells to prolong its fine ghosts. The multi-rank
  restrict-prolong carries a **7.3e-12 conservation defect** (vs 2.2e-16 at np=1). For a
  sharp interface with surface tension this tiny gather/prolong inconsistency **seeds a
  growing spurious capillary current** at the seam → blowup. For no-ST (smooth solution)
  the defect is negligible, which is why the no-ST controls stay sane.
- The blowup is **identical** across three attempted capillary fixes (c_divs ghost-extrap;
  c_divs recomputed from the prolonged color; wiring in the dead-code rank-seam fine halo
  `s_mpi_sendrecv_amr_fine_halo`). It is **independent of the capillary-path handling** —
  the seed is upstream in the AMR gather/prolong. `amr_max_blocks=2` (block tiled 2-ways)
  **crashes at init (SIGABRT)** — the tiling+ST path is separately broken.

**Conclusion:** the deadlock is a clean, isolated, fixable bug (fixed). Full multi-rank
AMR+ST **correctness** needs bit-consistency in the PR's core multi-rank gather/prolong
(and a working tiled path) — deep, unfinished PR-internal work, amplified by ST's sharp
interface. This is precisely why the PR gated ST under AMR; it is not a bounded fix.

### Feasibility consequence

At np=1's ~2.25 s/step, the handoff's spec (≥10τ ≈ 1.07M steps at R/50, c_l=100) is
~26 days **per run**; even 1τ is ~2.8 days. The uniform baselines parallelize freely and
are cheap, but the AMR run is pinned to np=1 by the halo bug and is the bottleneck.

The experiment as written assumed AMR+ST would simply run once the gate was lifted. It
does not. Chosen path: (a) reduced-scale np=1 containment test (this document) + this
characterization as PR feedback.

## Phase 1 — static-drop containment test (PRELIMINARY: hypothesis supported)

Reduced config (README): c_l=12, ±3.2R, R/50, 2τ, 100 saves. Trio: uniform coarse (320²),
uniform fine (640²), AMR (320² + static block cells 85–235 ≈ ±1.5R; interface at R stays
25 coarse cells inside → cf within 2.7e-11 of {0,1} across the seam). All at np=1/32r;
the AMR run is np=1 (halo bug).

**Rigorous metric:** AMR **seam-band** max|u| vs uniform-coarse max|u| in the *same* cell
band (control) — isolates any current the seam ADDS from the ambient parasitic flow that
reaches ±1.5R anyway. Coarse and AMR share the identical 320² base grid, so the band is
identical in cell indices.

The AMR **seam-band max|u| plateaus flat** almost immediately and stays there
(t/τ = 0.06 / 0.12 / 0.18 / 0.24 / 0.30 / 0.34 → 0.018 / 0.020 / 0.024 / 0.022 / 0.022 /
0.021 m/s) — no growth. Meanwhile the uniform-coarse parasitic flow keeps developing, so
the AMR-seam / coarse-same-band ratio *falls* over time:

| quantity | t/τ=0.16 | t/τ=0.34 |
|---|---|---|
| AMR seam-band max\|u\| | 0.023 m/s | 0.021 m/s (**flat**) |
| uniform-coarse **same-band** max\|u\| (control) | 0.032 m/s | 0.042 m/s (rising) |
| **AMR / coarse-same-band** | 0.71× | **0.49×** |
| AMR seam growth (2nd/1st half) | 0.95 | 1.06 (**non-growing**) |
| containment audit (cf dist to {0,1} in band) | 2.7e-11 | 2.7e-11 (**pass**) |

(uniform-coarse domain parasitic current reaches ~1.0 U_σ ≈ 0.48 m/s at 2τ, plateauing;
Laplace jump error 0.07% coarse / 0.7% amr.)

**Reduced-config caveat (fine run).** The uniform-fine (R/100) run developed a
near-cavitation spot (ρ→~0 spikes the viscous CFL) and its adaptive dt collapsed to ~4e-10
at t/τ≈1.6, so it was stopped there (context-only; 79 saves kept). This is a low-c_l=12 +
low-pressure artifact at R/100, not a containment issue — the **coarse** run (R/50, which
the verdict is built on) completed a clean 2τ, and the **AMR** run's dt is coarse-level-
controlled (R/50) so it runs stably at 2.6e-8. Watch the AMR fine block (also R/100) for
the same spot; so far (t/τ=0.5) it is healthy and the containment audit holds.

**The AMR block introduces no seam current.** Its seam-band velocity is a flat plateau
at/below the uniform run's ambient parasitic flow in the same region (the fine block
resolves the interface better, so *less* parasitic flow spreads to ±1.5R), non-growing,
with cf exactly {0,1} at the seam — the opposite of the §7.1 27–540× growing failure. The
containment hypothesis is **supported**. The flat plateau is established by t/τ≈0.06 and
holds through 0.34 (verified); the full 2τ run (grinding ~5h at np=1) extends it. Verdict
artifact: `results/summary.json`; overlay `results/figures/phase1_overlay.png`.
