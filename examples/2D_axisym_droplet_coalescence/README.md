# Film-zoom validation probe (axisymmetric, head-on vs symmetry plane)

**Purpose.** Produce an **in-house DNS reference for the dimple / early-drainage
phase** of the gas film in a head-on droplet collision, to validate the sub-grid
film model of `docs/gas_film_drainage_bouncing/DESIGN.md` (the capillary-compliance
closure) over the thickness window where MFC's continuum solution is trustworthy.
This is **not** a rupture-resolving run: continuum breaks below h ~ λ = 66 nm
(Kn > 1), and the grid numerically merges the interfaces near h ~ 2–3·dx anyway.
The validation window is **h ≈ 2 µm → ~0.3–0.45 µm** (pilot) — exactly the
dimple-formation phase (Chen & Yang 2020, Figs. 13/15).

**Physics.** Chen & Yang Fig. 6 head-on bounce point: We=8.6, Re=105.9, B=0,
D=306 µm, tetradecane/N₂. Same artificial-sound-speed EOS family as the
`2D_droplet_coalescence` case-s runs (c_l = 100 m/s; see `case.py` header for the
artificial-Mach caveat table and the 30 m/s film-velocity abort criterion, GATE-3).
Single drop + reflective plane at x=0 replaces the mirror drop (exact for head-on).
**Full-gap convention: H = 2·h_plane** — compare H, not h_plane, to references.

**Supersedes** the prior Phase-1 sketch in this directory (see `HANDOFF.md`),
which had three fatal flaws: 10 nm cells × real sound speed c_l≈1316 m/s → ~10⁸
acoustic steps (infeasible); **no viscosity** (film drainage *is* gas viscosity);
no symmetry-plane half-domain. Its validation target is retained: **Klaseboer
(2000) h_min(t) within ±15%**, alongside Chen & Yang Figs. 13/15 shapes and the
DESIGN.md compliance-closure prediction over the overlap window.

## Grid (verified by `grid_check.py` — run it after any parameter change)

| | pilot | production |
|---|---|---|
| dx at plane | 150 nm | 75 nm |
| interior dy (rim zone, r ≤ 90 µm) | 0.50 µm | 0.50 µm |
| cells (m+1)×(n+1) | 278k | 556k |
| steps to 0.4 ms (CFL 0.1) | 2.7M | 5.4M |
| est. wall @ 0.01 s/step, 64 ranks | ~7.5 h | ~15 h |

Stretch anatomy: MFC's tanh map has slope exactly 1 at x=0, so the plane cell
equals the pre-stretch spacing; two x-loops expand the far field to ~1.5 mm while
the drop stays ≤ 3 µm-resolved. Radial: axisym half-cell-offset grid, uniform
0.5 µm through the film-disk/rim zone (r_h grows to ~0.2D ≈ 61 µm).

## Gates (run in order; each gates the next)

1. **GATE-1 — static test** (`--mode static`, ~30 min): drop at rest, surface
   1 µm from the plane. Pass: (i) **Laplace jump = 2σ/R = 347.7 Pa ± 10%**
   measured across the interface — this decides whether MFC's surface tension
   carries the axisymmetric geometric terms (**no shipped example validates
   cyl_coord + surface_tension**; a ~σ/R result means the out-of-plane curvature
   is missing → axisym plan falls back to 2D planar); (ii) max spurious |u| <
   0.025 m/s (5 % of Ud) and not growing.
2. **GATE-2 — grid materialization**: run `pre_process`, check the printed
   min/max grid against `grid_check.py` (and `python3 case.py --check`).
3. **GATE-3 — film Mach monitor** (during the probe run): max gas |u| in the
   film band < 30 m/s (artificial Ma 0.3). Above that, the artificial-sound-speed
   trick is polluting the film — stop and reassess (raising c_l costs ∝ c_l steps).

## Run

```sh
cd ~/repos/MFC
./mfc.sh run examples/2D_axisym_droplet_coalescence/case.py -t pre_process simulation \
    --case-optimization -n 64 -- --mode static          # GATE-1
./mfc.sh run ... -- --mode probe                        # pilot (dx 150 nm)
./mfc.sh run ... -- --mode probe --dx-fine 75e-9        # production
# extend a finished run: raise --t-stop and pass --n_start <last step>
```

## Measurements (`analyze_film.py`)

From each frame: interface position x_i(r) (α=0.5 crossing nearest the plane) →
h_plane(r) → **H(r,t) = 2·h_plane**; H_center = H(r≈0); **H_min** and its radius
r_h (the dimple rim); gauge pressure p(0,0) (probe 1, high-rate) and p(r) on the
plane. Deliverables: `film_metrics.csv` + 4-panel figure (H_center, H_min, r_h,
p(0,0) vs t/(D/U_r)) for direct comparison with Chen & Yang Figs. 13/15 and the
compliance closure. Record the numerical-merge time (H_min → 0) — it marks the
end of the validation window, not a physical event.
