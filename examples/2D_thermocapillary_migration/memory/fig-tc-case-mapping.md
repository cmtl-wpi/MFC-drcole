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

Gotcha (fixed 2026-06-19): the `run.py` and `measure.py` header docstrings used to say
`fig5 (case_Ma_0.py)`, but the live fig5 target is `case_Ma_0p001.py` (Ma=0.001, conduction
TC1). Docstrings now corrected.

`case_Ma_0.py` is NOT orphaned — it is a distinct, README-documented case: the literal
**exactly Ma=0 frozen-T** reference (temperature frozen at the imposed linear profile, no
conduction; fast). It is run manually (`./mfc.sh run case_Ma_0.py`), not via run.py's automated
TARGETS, and feeds `runs/tc1/frozen/`. README case-file table + frozen-T grid-sweep section
cover it. Do not delete.
