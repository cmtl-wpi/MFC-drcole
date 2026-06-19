---
name: fig-tc-case-mapping
description: How run.py figure targets map to run-dir axis names and source case files
metadata:
  type: reference
---

`run.py` targets are named by **figure**, but run dirs are named by **Samareh case** (`tc1/`,
`tc2/`, `tc3/`). The mapping is not stated in the CLAUDE.md run-dir section, so `run.py fig5`
landing in `runs/tc1/` is non-obvious. Verified from `run.py` (2026-06-19):

| `run.py` target | run-dir axis | source case file | physics |
|---|---|---|---|
| `fig5` | `tc1` | `case_Ma_0p001.py` | 2D Ma→0 rise, grid convergence (v/v_YGB → ~0.80) |
| `fig7` | `tc2` | `case_Ma_20.py` | 2D low-Ma migration (U*/U_r, peak ~0.13) |
| `tc3` | `tc3` | `../3D_thermocapillary_migration/case_Ma_1723.py` | 3D large-Ma + mu(T) |

Gotcha: `run.py`'s header docstring says `fig5  case_Ma_0.py`, but the config dict actually
uses `case_Ma_0p001.py`. So `case_Ma_0.py` looks orphaned (present on disk, not the live fig5
case). The docstring is stale — trust the dict.
