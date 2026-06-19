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

Runs are nested **one subdir per axis**, in this order:

```
runs/<case>/<Ma>/<grid>/<smoothing>/
            tc1   ma0p1  w064   sc050
            tc2   ma0p01 w128   sc099
            tc3   ma0p001 w256
```

- `<case>` — Samareh case: `tc1` (Ma→0 rise), `tc2` (low-Ma), `tc3` (3D large-Ma; in the 3D dir)
- `<Ma>` — Marangoni number: `ma0p1`, `ma0p01`, `ma0p001` (`p` = decimal point)
- `<grid>` — cells per drop diameter: `w064`, `w128`, `w256`
- `<smoothing>` — interface `smooth_coeff`: `sc050` (=0.5, the **default/baseline**), `sc099` (=0.99, ~2× sharper)

Current tree:
```
runs/tc1/
  frozen/w064/                 frozen-T density proxy (Ma=0; no Ma/smoothing axis — the one asymmetry)
  ma0p1/   w064/{sc050,sc099}   w128/sc050
  ma0p01/  w064/{sc050,sc099}
  ma0p001/ w064/{sc050,sc099}
```

A leaf dir IS a single run (contains `restart_data/`, `*.inp`, the `case_*.py` copy). A run
only feeds a figure if its path matches what the plot scripts expect — `plot.py` (fig5/fig7)
and `plot_sc099.py` hardcode these leaf paths; `plot.py samareh`'s field loop walks the tree
for any dir with `simulation.inp`. Grep before adding a run (`grep -r tc1/ plot*.py run*.py`).

`smooth_coeff=0.5` is the default, so a baseline run with no explicit smoothing **is** the
`sc050` leaf — don't create a separate unsuffixed sibling.

`runs/2D_thermocapillary_migration-<timestamp>` is an `--archive` snapshot — leave it alone.

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
