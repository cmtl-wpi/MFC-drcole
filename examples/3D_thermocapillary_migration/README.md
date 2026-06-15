# 3D thermocapillary droplet migration — Samareh et al. (2014), Fig 6

The **3D sibling** of [`../2D_thermocapillary_migration`](../2D_thermocapillary_migration). That example
reproduces Samareh's **Fig 5** (a planar 2D *cylinder* → `v_t/v_YGB ≈ 0.80`); this one reproduces their
**Fig 6** — a fully-3D **sphere** in the same imposed linear temperature gradient, converging to
`v_t/v_YGB ≈ 0.95`.

> B. Samareh, J. Mostaghimi, C. Moreau, "Thermocapillary migration of a deformable droplet",
> *Int. J. Heat Mass Transfer* **73** (2014) 616–626.

## Why this needs bulk conduction (and the 2D case does not)

In the no-conduction (frozen-`T`) limit the **2D** rise reaches a quasi-steady plateau, but the **3D**
rise does **not**: the drop's toroidal internal circulation continuously steepens the frozen interfacial
temperature gradient, so the velocity drifts past `v_YGB` without saturating (finer grid → faster). With
no plateau there is no validatable 3D number.

**Bulk Fourier conduction** (the `thermal_conduction` feature) diffuses the temperature, tames that
runaway, and restores a steady plateau that can be compared to Samareh's `0.95`. So this example runs
with `thermal_conduction = T` and isothermal Dirichlet `y`-walls **by construction** — `SAMAREH3D_MA > 0`
is enforced in `case.py`.

## Case setup (Samareh §4.1.1, 3D)

| Quantity | Value |
|---|---|
| Drop | `D = 1` **sphere** (`geometry = 8`), 1.5D above the cold floor at `(0, −2.25, 0)` |
| Box | `5D × 5D × 7.5D` (`geometry = 9` cuboid background), gradient along **+y** |
| BCs | slip (`−2`) on all six faces; isothermal `T`-walls on `y` (cold floor / hot ceiling) |
| Fluids | `ρ_d = ρ_b = 0.2`, `μ_d = μ_b = 0.1`, `σ_0 = 0.1`, `σ_T = −0.1`, `\|∇T\| = 2/15` |
| `v_YGB` | `8.889×10⁻³`; `τ = 0.5`, `t_r = 7.5` |

Two patches share the same analytic density `ρ(y) = ρ_coeff/(T₀+∇T·y)` (no real density jump — the
`μ*=k*=1` YGB limit); only the color function differs. `T` is EOS-derived from density, shifted up by
`T₀ = 10` so the proxy stays positive (the gradient and slope `σ_T` are exact).

| Env var | Meaning | Default |
|---|---|---|
| `SAMAREH3D_NX` | cells per box width (Samareh used 64, 128) | 64 |
| `SAMAREH3D_MA` | thermal Marangoni number; **must be > 0** (conduction required) | 1.0 |
| `SAMAREH3D_TR` | run length in capillary-thermal times `t_r` | 2 |

## Quick start

```bash
# 3D is compute-heavy — use a decomposition with >= ~25 cells per rank-block per split dim.
./mfc.sh run examples/3D_thermocapillary_migration/case.py -n 8
python3 examples/3D_thermocapillary_migration/measure.py examples/3D_thermocapillary_migration
```

`measure.py` is dimension-agnostic (it reshapes to `(nz, ny, nx)` and auto-detects `dim`); it prints the
lab-frame color-weighted rise velocity and the quasi-steady plateau against `v_t/v_YGB ≈ 0.95`.

## Scope

A *converged* Fig-6 number is a heavy run: a 3D plateau needs ~1–2 `t_r` of wall-clock at a real grid
(64/128 per width), plus a grid-convergence pair — multi-day, not interactive. `case.py` validates and
smoke-runs immediately; the headline `0.95` comparison is a follow-up production run.
