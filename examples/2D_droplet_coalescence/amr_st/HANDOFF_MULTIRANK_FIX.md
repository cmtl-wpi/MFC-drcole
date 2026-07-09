# Handoff — implementing multi-rank AMR + surface_tension (upstream PR 1628)

**Goal:** make `surface_tension` work under the PR-1628 block AMR at **np ≥ 2** (it works
at np=1 today after the deadlock fix). This finishes the PR authors' unfinished
multi-rank machinery; it is a **contribution-to-the-PR project**, not a fix for the
coalescence runs (see "Reality check" at the end).

**Provenance / pins.** PR https://github.com/MFlowCode/MFC/pull/1628, branch
`sbryngelson:up/mega`, base SHA `ace2285a`. Worktree `~/repos/MFC-amr-st`, branch
`fix/amr-st-multirank-halo` (off a `test-enabler` commit that lifts the ST-under-AMR gate
in `m_checker.fpp` + `case_validator.py`). The PR is moving — pin `ace2285a` and record it
in every manifest. **Read first:** `results/FINDINGS.md` (the full diagnosis this builds
on) and `.claude/RESEARCH_WORKFLOWS.md`.

---

## State on entry (what is already true)

- **Deadlock: FIXED** (commit `b8c17198`). The capillary buffer populate
  (`s_get_capillary → s_populate_capillary_buffers`, `m_surface_tension.fpp`) ran a
  coarse-decomposition MPI halo on `c_divs` *inside* the AMR fine-block RHS; only
  block-owning ranks reach the fine advance → no MPI partner → deadlock. Fix: during the
  fine advance (`amr_in_fine_advance`, `m_global_parameters.fpp:278`) recompute `c_divs`
  over the ghost shell from the prolonged color and skip the populate — the exact pattern
  IGR already uses (`m_igr.fpp:306`). np=1 unregressed. **This stays.**
- **Containment (the §7.1 seam-normal problem): SOLVED at np=1, with data.** Keeping the
  interface ≥ (stencil+margin) cells inside the block edge removes the 27–540× c/f-seam
  current (`results/summary.json`: seam-band max|u| 0.3–0.7× the parasitic baseline,
  non-growing). This is Barrier 2 below.
- **Blowup at np ≥ 2: OPEN.** With the deadlock gone, np=2 blows up
  (alpha1→4.5, pres→2.5e5) at the interface∩rank-boundary. Root cause traced (below).

---

## The three barriers (verified facts, with pointers)

**Barrier 1 — capillary MPI deadlock.** Done (`b8c17198`). Nothing left.

**Barrier 2 — coarse/fine seam-normal (§7.1).** At the block *edge*, fine ghost color is
prolonged from coarse → inconsistent capillary normal → seam current. **Resolved by
containment** (keep the interface off the block edge). Needs to become a *guard* (Phase 2).

**Barrier 3 — multi-rank consistency (the blowup).** Two sub-cases:
- **3a. Whole-owned block** (`amr_max_blocks < num_procs`; the DEFAULT). One rank
  `s_amr_gather_coarse_patch` (`m_amr.fpp:493`) gathers the *other* ranks' coarse cells.
  The gather values are bit-exact, BUT the round-trip `restrict(prolong(·))`
  (`s_amr_conservation_check`, `m_amr.fpp:1750`, prints `[amr] restrict-prolong
  conservation err`) is **2.2e-16 at np=1 but 7.3e-12 at np=2**. The 7.3e-12 lives in
  **total energy at a cell on the interface** (largest prolongation slope → catastrophic
  FP cancellation of `u0 ± 0.25·slope` in the restrict sum). ST amplifies this ~1e-11 seed
  into the blowup. The prolong (`s_amr_fill_fine_ghosts` :2066, `s_prolong_one_var` :942)
  and restrict (`s_restrict_one_var` :1162) are **coordinate-free, index-based, and exactly
  conservative in exact arithmetic** — this is a floating-point *reproducibility* gap
  between the np=1 and multi-rank paths, not an operator bug. Exact ULP source not yet
  isolated (see Phase 0).
- **3b. Distributed/split block** (`amr_max_blocks ≥ num_procs`). Fine ghosts at rank/tile
  seams are still *prolonged from coarse*, not exchanged with the neighbour's actual fine
  cells (`s_amr_fine_stage_fill` comment, `m_amr.fpp:2488` — "final ghost state EXCEPT at
  faces shared with an adjacent sub-block", i.e. only tile seams via `s_amr_fine_fine_halo`
  :2391; rank seams are unhandled). The intended fix, `s_mpi_sendrecv_amr_fine_halo`
  (`m_mpi_proxy.fpp:82`), **exists but is DEAD CODE (never called)**. Wiring it naively into
  Phase 2 **crashed at init (SIGABRT)** with a tiled block — it needs proper integration.

---

## The feasible path (do in order; each gates the next)

### Phase 0 — nail the exact ULP divergence (½–1 day, do FIRST)
Before writing any fix, isolate WHERE the np=1 vs np=2 floating-point paths first diverge.
Do NOT debug blind.
- Instrument `s_amr_gather_coarse_patch` / `s_amr_conservation_check` to dump `amr_cg` and
  the freshly-prolonged fine block **bit-for-bit** (hex or full-precision) at init, np=1 vs
  np=2, for `var=5` (total energy) around the interface cell the check flags. Compare.
- Distinguish the two hypotheses: (a) a genuine 1-ULP *data* difference in the gather/MPI
  path, vs (b) a *compiler* FMA/vectorization difference on the multi-rank code path (test
  by rebuilding with `-fp-model=precise`/`-no-fma` for ifx and re-checking the 7.3e-12).
- `gdb` attach is blocked here (`ptrace_scope=1`, no sudo) — use source prints to stderr
  (`write(0,*)`, unbuffered) as in this session, or ask to set `ptrace_scope=0`.
- **Decision gate:** if (b) compiler FMA, the "fix" may be a build-flag/pragma on the
  prolong/restrict kernels (cheap). If (a) data, proceed to Phase 3a.

### Phase 1 — keep the deadlock fix, add regression coverage
`b8c17198` is in. Add a test that AMR+ST at np=1 runs and matches a golden (the reduced
Laplace case, `cases/case_laplace.py --variant amr`). Confirm non-AMR ST is byte-identical
(the fix only touches the `amr_in_fine_advance` branch).

### Phase 2 — containment guard (Barrier 2; the containment result already proves it)
Relax the hard `@:PROHIBIT(surface_tension)` (currently deleted on the test-enabler commit)
to a **guard**: permit ST under AMR only when the tagged/blocked region keeps the diffuse
interface ≥ (WENO stencil + margin) cells inside every block edge; abort otherwise. Static
block: check `amr_block_beg/end` vs the interface at init. Dynamic regrid: dilate the
`amr_tag_eps` color-gradient tags by the margin (mirror the EL-cloud / moving-IB
contain-and-guard pattern the PR already ships). Put runtime checks in
`src/simulation/m_checker.fpp` (static geometry) and a per-regrid guard in `m_amr.fpp`.

### Phase 3 — multi-rank seam consistency (Barrier 3; the hard, bounded part)
**3a — distribute blocks.** Set `amr_max_blocks ≥ num_procs` (case param) so each rank owns
its LOCAL fine piece; no cross-rank coarse gather for the interior → removes the whole-owned
FP-gather seed. (If Phase 0 shows the gather itself is fixable bit-consistently, do that
instead so the default `amr_max_blocks` also works.)

**3b — finish `s_mpi_sendrecv_amr_fine_halo`** so fine ghosts at rank/tile seams carry the
neighbour's actual fine interior (not prolonged coarse). It is written for the mirror
decomposition (`amr_isect_lo/hi`, `amr_region_lo/hi`, self-guards on
`amr_rank_owns_block`). Wire it into Phase 2 of the step loop
(`m_time_steppers.fpp:509–528`, alongside `s_amr_fine_fine_halo`), per-owned-block. My naive
wiring crashed at init with a tiled block (`amr_max_blocks=2`) — likely an
extent/buffer-sizing or ownership-model (SFC `amr_block_owner` vs mirror
`amr_rank_owns_block`) mismatch, or the routine reads stage-entry data at the wrong phase.
Debug against the SIGABRT: build `--debug` (bounds-checking), reproduce at np=2
`amr_max_blocks=2`, read the backtrace. Once seams carry real fine data, the FP-amplified
blowup has nothing to feed on.

---

## Verification contract (a fix is not done until all pass)
1. **np-independence (the decisive test):** np=2 result matches np=1 to FP round-off
   (~1e-10 rel), not a blowup. Harness: run `cases/case_laplace.py --variant amr` at n=1 and
   n=2 from the same IC, compare the final level-0 checkpoint with
   `analysis/seam_analysis.py` / a field-diff (the exact script used this session:
   read both `restart_data/lustre_<N>.dat`, reshape `(8, n+1, m+1)`, compare rho/speed/
   alpha1/cf/pres). A CORRECT fix → alpha1∈[0,1], pres~285, |u|~0.16 at both n.
2. **Conservation:** `restrict-prolong conservation err` back to ~1e-16 at np=2 (or
   understood/accepted per Phase 0).
3. **No-ST unregressed:** AMR without ST still runs at np=16 (it does today).
4. **Non-AMR ST byte-identical** to base.
5. **np=16 and `amr_max_blocks≥num_procs`** run to completion (not just np=2).
6. Never regenerate golden files to hide a diff; a wrong seam fix is a *silent* wrong
   answer — the np-independence test is the guard.

## Tooling (what works here)
- Build: `scripts/build.sh` (Intel ifx; source setvars atop every build/run; `|| true` on
  setvars under `set -e`). Incremental: `./mfc.sh build -t simulation -j 32`.
- Run: **standalone-case pattern** — `scripts/gen_case.py <run_id> <dir> -- <case args>`
  writes a case.py with args BAKED into `sys.argv` (mfc.sh `-- <args>` forwarding is
  unreliable in detached/zsh shells; the tool shell here is zsh — `$VAR` does NOT
  word-split). Then `taskset -c <cores> ./mfc.sh run <dir>/case.py -n <ranks> -t
  pre_process simulation`. `mpirun` needs the sandbox disabled; watch stderr, not stdout.
- Reduced config for fast np=1/np=2 iteration (does NOT change the seam mechanism):
  `--variant amr --block-R 1.5 --c-l 12 --domain-R 3.2 --n-periods 0.02 --n-saves 4`
  (320²+block, ~256 steps, ~5 min at np=2). np=1 stays sane; np=2 blows up by save 4 today.
- Durable long runs: `scripts/run_durable.sh` (detached screen + self-healing resume loop,
  survives the session reaper).
- `amr_max_blocks` is a case param (`--amr-max-blocks N` in `case_laplace.py`).

## Reality check (read before committing weeks to this)
Two-level 2:1 AMR is ~2–4× at best in 3D and does **not** buy film-scale resolution at
rupture — so a fully-working ST+AMR would not materially help the coalescence study.
Treat this as upstream contribution, not a production speedup. The high-value, low-effort
alternative: send the PR authors (1) the deadlock fix, (2) the Barrier-3 diagnosis + the
7.3e-12 conservation-defect repro, (3) the containment result — that's exactly what lets
*them* finish the seam halo. **Do not post to GitHub without Davis's approval.**
