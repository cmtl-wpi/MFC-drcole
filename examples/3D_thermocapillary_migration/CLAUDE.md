# Run organization — 3D_thermocapillary_migration

How to keep simulation runs in this example organized. Physics and validation status live in
`README.md`; broader MFC run/data-lifecycle conventions live in `../../.claude/RESEARCH_WORKFLOWS.md`.
This file is only about *where output goes and what it's named*. It mirrors the 2D example's
`../2D_thermocapillary_migration/CLAUDE.md`, with this example's own sweep axes.

## Golden rule: never run in-place

MFC writes its output (`restart_data/`, `silo_hdf5/`, `D/`, `*.inp`, `MFC.out`, …) into the
directory it's launched from. **Always run a case into its own `runs/...` leaf, never directly in
the example root.** The root holds source only: `case*.py`, the helper scripts (`run_ygb.py`,
`measure.py`, `validate_ygb.py`, `fields_ygb.py`), `README.md`, and the `figures/` + `results/`
outputs. `run_ygb.py` does this for you — it copies the case into the leaf and launches there.

## The three case files

- `case_ygb.py` — the **u_YGB validation** (this directory's purpose). Decoupled `thermal_scalar`
  drives σ(T); uniform density. Two geometry modes (`YGB_GEOM`): `cube` (centered drop, the
  confinement sweep) and `samareh` (offset 5D×7.5D box, the ≈0.95 anchor). See [[why-thermal-scalar]].
- `case.py` — legacy Samareh Fig 6 baseline using the **density-proxy** temperature (no
  `thermal_scalar`). Kept for reference; not the validation path.
- `case_Ma_1723.py` — the experimental large-Ma (LMS Shuttle) case with SI units + μ(T) Arrhenius.

## Run-dir naming

`case_ygb.py` is parameterized by env vars (`YGB_GEOM/YGB_W/YGB_NX/YGB_MA/YGB_TR`), so a run's
identity is its sweep coordinates, encoded in the leaf path (one subdir per axis):

```
runs/ygb/<geom>/<W>/<grid>/<Ma>/
         cube     w8     nx080  ma0p5
         samareh  w10    nx128  ma0p25
```

- `<geom>` — `cube` (unbounded sweep) or `samareh` (confined anchor); separates them so a cube W5
  and the samareh W5 anchor never collide.
- `<W>` — cube box width in D: `w6`, `w8`, `w10`, `w12` (`p` = decimal point, e.g. `w7p5`).
- `<grid>` — cells per box width: `nx064`, `nx080`, `nx128`.
- `<Ma>` — Marangoni number: `ma1`, `ma0p5`, `ma0p25`, `ma0p1`.

A leaf dir IS one run (`restart_data/`, `*.inp`, the `case_ygb.py` copy). `validate_ygb.py` parses
the coordinates back out of the path, so the encoding is the source of truth — no separate manifest.

## Pipeline

```
run_ygb.py {smoke|anchor|confinement|grid|ma|all} [--force]   # run sims into runs/, measure each
measure.py <run_dir>                                          # rise velocity -> RESULT_JSON + PNG
validate_ygb.py {confinement|grid|ma|all}                     # convergence fits + figures
fields_ygb.py <run_dir> [step]                                # T_s / color midplane sanity field
```

**`run_ygb.py` does NOT rmtree a populated leaf** (the opposite of the 2D `run.py`). It SKIPS a leaf
that already has `restart_data/`, so a partial multi-day sweep is safe to re-invoke — you won't
clobber checkpoints. Use `--force` to deliberately rerun (it deletes first). The heavy selectors are
multi-hour-per-run; run them under `nohup` / in the background.

## Precious vs regenerable (inside a leaf)

- **Precious — never delete:** `restart_data/` (sole source of truth), the `case_ygb.py` copy, the
  `*.inp` inputs.
- **Regenerable bloat — safe to prune for disk:** `silo_hdf5/`, `viz/`, `run_time.inf`, `MFC.out`,
  `D/`. All rebuildable from `restart_data/` + inputs.

## Keep the root clean

If a case is ever run in-place here, remove the leftovers afterward (all gitignored, none tracked):
`restart_data/`, `D/`, `silo_hdf5/`, `viz/`, `*.inp`, `*.dat`, `MFC.out`, `MFC.sh`, `run_time.inf`,
`smoke_run.log`, `__pycache__/`. Don't commit them; don't leave them.

## Memory

This example has its own memory in `memory/`, co-located with this CLAUDE.md. Read `memory/MEMORY.md`
before working here, and keep it current. Example-specific facts go here; project-wide facts go in
`../../memory/`. Format and conventions are in `../../CLAUDE.md` (Memory section).
