# 3D thermocapillary droplet migration — Samareh et al. (2014), Fig 6

The **3D sibling** of [`../2D_thermocapillary_migration`](../2D_thermocapillary_migration). That example
reproduces Samareh's **Fig 5** (a planar 2D *cylinder* → $v_t/v_\mathrm{YGB} \approx 0.80$); this one reproduces their
**Fig 6** — a fully-3D **sphere** in the same imposed linear temperature gradient, converging to
$v_t/v_\mathrm{YGB} \approx 0.95$.

> B. Samareh, J. Mostaghimi, C. Moreau, "Thermocapillary migration of a deformable droplet",
> *Int. J. Heat Mass Transfer* **73** (2014) 616–626.

## Why this needs bulk conduction (and the 2D case does not)

In the no-conduction (frozen-$T$) limit the **2D** rise reaches a quasi-steady plateau, but the **3D**
rise does **not**: the drop's toroidal internal circulation continuously steepens the frozen interfacial
temperature gradient, so the velocity drifts past $v_\mathrm{YGB}$ without saturating (finer grid → faster). With
no plateau there is no validatable 3D number.

**Bulk Fourier conduction** (the `thermal_conduction` feature) diffuses the temperature, tames that
runaway, and restores a steady plateau that can be compared to Samareh's $0.95$. So this example runs
with `thermal_conduction = T` and isothermal Dirichlet `y`-walls **by construction** — `SAMAREH3D_MA > 0`
is enforced in `case.py`.

## Case setup (Samareh §4.1.1, 3D)

| Quantity | Value |
|---|---|
| Drop | $D = 1$ **sphere** (`geometry = 8`), 1.5D above the cold floor at $(0, -2.25, 0)$ |
| Box | $5D \times 5D \times 7.5D$ (`geometry = 9` cuboid background), gradient along **+y** |
| BCs | slip (`−2`) on all six faces; isothermal $T$-walls on $y$ (cold floor / hot ceiling) |
| Fluids | $\rho_d = \rho_b = 0.2$, $\mu_d = \mu_b = 0.1$, $\sigma_0 = 0.1$, $\sigma_T = -0.1$, $\|\nabla T\| = 2/15$ |
| $v_\mathrm{YGB}$ | $8.889 \times 10^{-3}$; $\tau = 0.5$, $t_r = 7.5$ |

Two patches share the same analytic density $\rho(y) = \rho_\mathrm{coeff}/(T_0+\nabla T \cdot y)$ (no real density jump — the
$\mu^*=k^*=1$ YGB limit); only the color function differs. $T$ is EOS-derived from density, shifted up by
$T_0 = 10$ so the proxy stays positive (the gradient and slope $\sigma_T$ are exact).

| Env var | Meaning | Default |
|---|---|---|
| `SAMAREH3D_NX` | cells per box width (Samareh used 64, 128) | 64 |
| `SAMAREH3D_MA` | thermal Marangoni number; **must be > 0** (conduction required) | 1.0 |
| `SAMAREH3D_TR` | run length in capillary-thermal times $t_r$ | 2 |

## Quick start

```bash
# 3D is compute-heavy — use a decomposition with >= ~25 cells per rank-block per split dim.
./mfc.sh run examples/3D_thermocapillary_migration/case.py -n 8
python3 examples/3D_thermocapillary_migration/measure.py examples/3D_thermocapillary_migration
```

`measure.py` is dimension-agnostic (it reshapes to `(nz, ny, nx)` and auto-detects `dim`); it prints the
lab-frame color-weighted rise velocity and the quasi-steady plateau against $v_t/v_\mathrm{YGB} \approx 0.95$.

## Scope

A *converged* Fig-6 number is a heavy run: a 3D plateau needs ~1–2 $t_r$ of wall-clock at a real grid
(64/128 per width), plus a grid-convergence pair — multi-day, not interactive. `case.py` validates and
smoke-runs immediately; the headline $0.95$ comparison is a follow-up production run.

# `case_ygb.py` — recovering $u_\mathrm{YGB}$ as a validation of variable surface tension

`case.py` (above) reproduces Samareh's *confined* $0.95$ using the **density-proxy** temperature.
`case_ygb.py` is a stronger, cleaner test: it recovers the **analytic** Young–Goldstein–Block
terminal velocity $u_\mathrm{YGB}$ ($v_t/u_\mathrm{YGB} \to 1.0$) as a proper validation of MFC's $\sigma(T)$ Marangoni stress.

## Why it differs from `case.py`

Both cases use the same temperature setup — the **density proxy** ($\rho(y)=\rho_\mathrm{coeff}/T(y)$, recovered
from the EOS) with bulk conduction + isothermal walls holding the gradient (identical fluids,
$\mu^*=k^*=1$, so the only driver is the $\sigma(T)$ gradient). The difference is the **validation strategy**:

- `case.py` is a single *confined* run in an offset $5D \times 7.5D$ box that reproduces Samareh's Fig 6
  anchor ($v_t/v_\mathrm{YGB} \approx 0.95$).
- `case_ygb.py` is a **convergence harness**: a cube geometry with the drop centered for clean
  symmetric clearance, plus `YGB_W`/`YGB_MA`/`YGB_NX` knobs to extrapolate the confinement, finite-Ma,
  and grid deficits → recover the analytic $u_\mathrm{YGB}$ ($v_t/u_\mathrm{YGB} \to 1.0$).

## Recovering $u_\mathrm{YGB}$ is a convergence claim, not one number

$u_\mathrm{YGB}$ is the *unbounded, zero-Ma, Stokes* sphere result. MFC falls below $1.0$ by three
*vanishable* deficits; the validation is showing each → 0 drives the ratio → $1.0$:

| deficit | knob | sweep | reduction |
|---|---|---|---|
| confinement | `YGB_W` (cube box width) | 6, 8, 10, 12 | fit ratio vs $1/W$, extrapolate $W \to \infty$ (**headline**) |
| finite Ma | `YGB_MA` | 1.0, 0.5, 0.25, 0.1 | extrapolate $\mathrm{Ma} \to 0$ (perfectly invariant $T$) |
| grid | `YGB_NX` | 64, 96, 128 | Richardson in $dx$ |

$\mathrm{Re}_M = \rho v_\mathrm{YGB} D/\mu \approx 0.018$ is already deep Stokes — no Reynolds sweep. Geometry modes via
`YGB_GEOM`: `cube` (centered drop, the sweep) and `samareh` (offset $5D \times 7.5D$ box, the $\approx 0.95$ anchor).

## Pipeline

```bash
python3 run_ygb.py smoke               # tiny machinery check (serial; minutes)
python3 run_ygb.py anchor              # reproduce Samareh ~0.95 in the confined box
python3 run_ygb.py confinement         # the headline cube sweep (multi-hour each; run under nohup)
python3 run_ygb.py grid ma             # grid + Ma refinement of the converged corner
python3 validate_ygb.py all            # convergence fits -> figures/ygb_vs_{confinement,dx,Ma}.png
python3 fields_ygb.py <run_dir>        # EOS temperature / color midplane sanity field
```

`run_ygb.py` runs into `runs/ygb/<geom>/<W>/<grid>/<Ma>/` and **skips** populated leaves (pass
`--force` to rerun) — a partial sweep is safe to re-invoke. The heavy selectors are a multi-day
background/batch campaign.
