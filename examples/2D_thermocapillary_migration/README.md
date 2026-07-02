# Thermocapillary droplet migration — σ(T) validation against Samareh et al. (2014)

Validation of MFC's **temperature-dependent surface tension** closure
(`sigma_model = 1`, linear σ(T)) against the thermocapillary-migration benchmarks of

> B. Samareh, J. Mostaghimi, C. Moreau, **"Thermocapillary migration of a deformable droplet,"**
> *Int. J. Heat Mass Transfer* **73** (2014) 616–626.

A neutrally-buoyant drop in an imposed linear temperature field develops a surface-tension gradient
along its interface (tension falls as temperature rises, σ_T = dσ/dT < 0). The resulting tangential
**Marangoni stress** drags interfacial fluid hot→cold and, by reaction, the drop migrates toward the
**hot** wall.

![Marangoni migration mechanism](figures/mechanism_schematic.png)

*Tension falls with temperature, so σ varies around the interface (high on the cold side, low on the
hot side). That tangential σ-gradient is the Marangoni stress ∇ₛσ; it drives interfacial fluid
hot→cold and, by reaction, the drop migrates toward the hot wall. Source:
`figures/mechanism_schematic.tex`.*

This example covers the two **σ(T)-only** Samareh cases (§4.1.1 and §4.1.2). The large-Marangoni LMS
case (§4.2) additionally needs temperature-dependent viscosity μ(T) and is not included here.

| | Samareh scenario | Regime | Figure | Result |
|---|---|---|---|---|
| **TC1** | Zero Marangoni number (§4.1.1) | Ma = 0 | Fig 5 (2D cylinder) | plateau v_t/v_YGB ≈ **0.80** |
| **TC2** | Low Marangoni, Nas & Tryggvason (§4.1.2) | Re = 5, Ma = 20, Ca = 0.0167 | Fig 7 | peak U* ≈ **0.13** |

## The σ(T) feature

Selected per case by `sigma_model = 1` (`sigma_model = 0` reproduces the constant-σ path bit-for-bit):

| Parameter | Meaning |
|---|---|
| `sigma_model` | 0 = constant σ, 1 = linear σ(T) |
| `sigma_T_ref`  | reference temperature T_ref |
| `sigma_dTdT`   | slope σ_T = dσ/dT (signed) |

The closure fills a cell-centered field `σ(T) = σ + sigma_dTdT·(T − sigma_T_ref)`; the capillary
continuum-surface-force flux then uses a **face-local** σ (averaged from the two adjacent cells), so
the tangential variation of σ along the interface *is* the Marangoni stress — added to the momentum
and CSF-energy source terms. Temperature is recovered from the stiffened-gas mixture EOS
(`f_compute_mixture_temperature`); there is no separate temperature field.

**Imposing the gradient (density proxy).** MFC is compressible, so T = (p+p∞)/((γ−1)ρcv) is an EOS
function of the state. The linear T(y) is imposed through the density IC, ρ(y) = ρ_coeff/T(y), painted
with a single full-box analytic patch so T stays continuous across the drop interface (a two-patch
seam injects a spurious force exactly where the Marangoni force lives —
[`case_Ma_20_2patch.py`](case_Ma_20_2patch.py) is the control that measures that error).

**Dependency.** σ(T) needs a sustained temperature gradient. A bare density proxy is advected away by
the flow (good only for t/t_r ≲ 2), so these cases run with **bulk thermal conduction**
(`thermal_conduction = T`) + isothermal walls to hold the field — i.e. this example builds on the
thermal-conduction feature.

## Cases

| Case | Samareh | What it is |
|---|---|---|
| [`case_Ma_0.py`](case_Ma_0.py)            | §4.1.1 / Fig 5 | 2D Ma=0 frozen-T rise, slip-wall box. Literal Ma=0 reference (fast); drifts above 0.80 at late t/t_r. |
| [`case_Ma_0p001.py`](case_Ma_0p001.py)    | §4.1.1 / Fig 5 | Conduction companion (deep Ma=0.001 limit) + isothermal walls actively holding the imposed gradient. |
| [`case_Ma_0p1.py`](case_Ma_0p1.py)        | §4.1.1 / Fig 5 | Conduction companion at Ma=0.1. |
| [`case_Ma_20.py`](case_Ma_20.py)          | §4.1.2 / Fig 7 | Finite-Ma migration (Nas & Tryggvason), two-fluid drop + conduction. |
| [`case_Ma_20_2patch.py`](case_Ma_20_2patch.py) | §4.1.2 | Control: TC2 with a conventional two-patch drop IC, to isolate the single-patch requirement. |

## Results

**TC1 / Fig 5 — 2D cylinder, Ma = 0 → 0.80.** The committed sweep (`results/summary.json`) runs the
bulk-conduction case at 12.8, 25.6, and 51.2 cells/D, each to t/t_r = 2. The terminal-velocity
plateaus are v_t/v_YGB = 0.81, 0.89, and 0.87 — a non-monotone spread above Samareh's 0.80 within
this window; only the coarsest grid lands on the reference value, and coincidentally so. The
unbounded-cylinder analytic limit is 15/16 = 0.938; the finite slip-wall box brings that down to
≈ 0.80, the value Samareh's own four 2D methods agree on (it is a cylinder-in-a-box, not a sphere).

**TC2 / Fig 7 — finite Ma.** [`case_Ma_20.py`](case_Ma_20.py) runs the Nas & Tryggvason
configuration (Re = 5, Ma = 20, Ca = 0.0167): U* = U/U_r ramps from rest, overshoots, and relaxes
toward a terminal value; the published transient peaks at U* ≈ 0.13.
[`case_Ma_20_2patch.py`](case_Ma_20_2patch.py) repeats it with a conventional two-patch drop IC,
isolating the spurious interfacial force a patch seam injects. `run.py fig7` runs the grid variants
and `measure.py` reduces each to a U*(t*) curve overlaid on the digitised transient
(`figures/case2_fig7.png`); no TC2 results are committed here.

## Caveats

- **Compressible acoustics.** The closed slip-wall box rings (the unbalanced t=0 Laplace jump excites
  the vertical acoustic standing wave). It is a resolved oscillation, not aliasing; read the migration
  as the mean of the velocity cloud. Samareh's incompressible solver has no acoustics.
- **Frozen-T drift.** Without conduction the density-proxy gradient is advected, so `case_Ma_0` ramps,
  plateaus, then drifts — compare the post-overshoot plateau, not the endpoint. The conduction runs
  (the headline figures) hold the imposed gradient instead, though their tails still move on the
  finer grids (`slope_per_tr` in `results/summary.json`).
- **Distinct-fluid finite-Ma is qualitative.** Each fluid's absolute density stratifies ~1/T (a
  compressibility artifact absent in the incompressible reference), so TC2's exact terminal value is
  fidelity-limited; TC1's equal-density case avoids this artifact.

## Reproducing

```bash
# Single TC1 run + measure the migration velocity:
./mfc.sh run examples/2D_thermocapillary_migration/case_Ma_0.py -n 6
python3 examples/2D_thermocapillary_migration/measure.py examples/2D_thermocapillary_migration

# Grid sweeps (launch MPI jobs into runs/; fig5 aggregates results/summary.json,
# fig7 results/fig7_summary.json):
python3 examples/2D_thermocapillary_migration/run.py fig5   # TC1 Fig 5 grids
python3 examples/2D_thermocapillary_migration/run.py fig7   # TC2 Fig 7
# run.py <fig5|fig7> [run|remeasure]; `remeasure` re-reads existing runs/ without simulating.

# Field plots (EOS-recovered T, σ around the interface, recirculation):
python3 examples/2D_thermocapillary_migration/plot.py fields <case_dir> temperature
```

`runs/` is gitignored simulation output.
