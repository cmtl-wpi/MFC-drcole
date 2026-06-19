# Run organization — 2D_thermocapillary_migration

How to keep simulation runs in this example directory organized. Physics, validation
status, and figure provenance live in `README.md`; broader MFC run/data-lifecycle
conventions live in `../../.claude/RESEARCH_WORKFLOWS.md`. This file is only about *where
output goes and what it's named*.

## Golden rule: never run in-place

MFC writes its output (`restart_data/`, `silo_hdf5/`, `D/`, `*.inp`, `MFC.out`, …) into the
directory it's launched from. **Always run a case into its own `runs/<name>/` dir, never
directly in the example root.** The root holds source only: `case_*.py`, the helper scripts
(`run.py`, `run_sc099.py`, `measure.py`, `plot.py`, `plot_sc099.py`, `compare_tc3_visc.py`),
`README.md`, `samareh2014.pdf`, and the `figures/` + `animations/` outputs.

`run.py` already does this for you — it copies the case into `runs/<name>/`, rewrites the
`Nx = <n>` line for the grid, and launches there. Prefer it over manual `mfc.sh run`.

## Run-dir naming

General grammar: `tc<N>_<descriptor>_w<grid>[_ma<value>][_sc<pct>]`
- `tc<N>` — Samareh case: tc1 (Ma→0 rise), tc2 (low-Ma), tc3 (3D large-Ma, in the 3D dir)
- `w<grid>` — cells per drop diameter (`w064`, `w128`, `w256`)
- `_ma<value>` — Marangoni number when swept (`ma0p1`, `ma0p01`, `ma0p001`)
- `_sc<pct>` — interface smoothing `smooth_coeff` when studied (`sc050`, `sc099`)
- `<descriptor>` — `cond` (bulk conduction + independent T_s) or `frozen` (frozen-T density proxy)

⚠️ **Two naming schemes currently coexist — reconcile before adding runs:**
- `run.py` emits short names with Ma implied by the case file: `tc1_w064/128/256` (fig5,
  from `case_Ma_0p001.py`), `tc2_w064/128` (fig7, from `case_Ma_20.py`).
- The runs actually on disk and the names `plot.py`/`plot_sc099.py` read use the long form:
  `tc1_frozen_w064`, `tc1_cond_w064_ma0p1`, `tc1_cond_w128_ma0p1`, `tc1_cond_w064_ma0p1_sc099`, …

These do **not** match. A run only feeds a figure if its dir name matches what the plot
script hardcodes — grep the script (`grep tc1_ plot.py`) before naming a new run, or you'll
produce output the plotter can't find. When in doubt, follow the long form (it encodes Ma
and smoothing explicitly) and make `run.py`'s variant table agree with it.

`2D_thermocapillary_migration-<timestamp>` is an `--archive` snapshot — leave those alone.

## Pipeline

```
run.py {fig5|fig7|tc3|all} [run|remeasure]   # run sims into runs/, then measure + replot
measure.py <run_dir> {fig5|fig7|tc3}         # extract rise velocity -> RESULT_JSON + PNG
plot.py {samareh|ma|fields|clean}            # overlays / Ma sweep / field maps / orphan cleanup
```

- `remeasure` re-measures existing `runs/` and replots **without** re-running any sim.
- `run.py` `shutil.rmtree`s the run dir before each run — **do not** use it to resume; you'll
  delete the checkpoints. Re-running a sim also wipes that run's `silo_hdf5/`.

## Precious vs regenerable (inside a run dir)

- **Precious — never delete:** `restart_data/` (sole source of truth), the `case_*.py` copy,
  the `*.inp` inputs.
- **Regenerable bloat — safe to prune for disk:** `silo_hdf5/` (~78M), `viz/`, `run_time.inf`
  (~65M), `MFC.out`/`run.log` (~59M), `D/`. All rebuildable from `restart_data/` + inputs.

## Keep the root clean

If a case ever gets run in-place here, remove the leftovers afterward (all gitignored, none
tracked): `restart_data/`, `D/`, `silo_hdf5/`, `viz/`, `out/`, `*.inp`, `*.dat`, `MFC.out`,
`MFC.sh`, `run_time.inf`, `indices.dat`, `__pycache__/`. Don't commit them; don't leave them.
