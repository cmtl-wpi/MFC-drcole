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

### The likely mechanism (`m_mpi_proxy.fpp:181`)
The fine-halo `MPI_SENDRECV` uses the **same `cnt` for send and receive**, assuming the
neighbor's matching fine face has identical size and that the neighbor participates. Under
the mirror decomposition adjacent ranks can own different-sized block pieces, and a rank
owning no block cells returns early (line 92) — so a poster can lack a partner (deadlock)
or a size can mismatch (truncation). Why this only bites with ST is not yet pinned down
(the count scales uniformly with `sys_size`); a plausible cause is a separate
color-function/curvature halo on the fine level. Diagnosing the exact trigger is deferred.

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

Result (AMR at t/τ=0.16 so far; coarse/finished):

| quantity | value |
|---|---|
| uniform-coarse domain parasitic current (2τ) | ~1.0 U_σ ≈ 0.48 m/s, plateauing (growth 2.1→1.58) |
| AMR seam-band max\|u\| | 0.023 m/s |
| uniform-coarse **same-band** max\|u\| (control) | 0.032 m/s |
| **AMR / coarse-same-band** | **0.71×** (AMR is *lower*) |
| AMR seam growth (2nd/1st half) | 0.95 (**non-growing**) |
| containment audit (cf dist to {0,1} in band) | 2.7e-11 (**pass**) |
| Laplace jump error | 0.07% (coarse), 0.7% (amr) |

**The AMR block introduces no seam current.** Its seam-band velocity is at/below the
uniform run's ambient parasitic flow in the same region (the fine block resolves the
interface better, so *less* parasitic flow spreads to ±1.5R), non-growing, with cf exactly
{0,1} at the seam — the opposite of the §7.1 27–540× growing failure. This supports the
containment hypothesis. **Caveat:** the AMR run is at t/τ=0.16 of 2τ; the full run (grinding
~6h at np=1) is needed to confirm the ratio stays O(1) through the parasitic plateau.
Verdict artifact: `results/summary.json`; overlay `results/figures/phase1_overlay.png`.
