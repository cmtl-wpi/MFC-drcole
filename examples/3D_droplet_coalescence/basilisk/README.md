# Basilisk — 3D binary droplet collision (case p)

A [Basilisk](http://basilisk.fr) port of `examples/3D_droplet_coalescence`
**case p** (Qian & Law tetradecane-in-nitrogen collision). Same physics, solved
the way Basilisk is built for it: **incompressible two-phase Navier–Stokes** with
a geometric VOF interface, height-function surface tension, and **octree AMR**.

## Case p

Near-grazing, high-Weber collision — the stretching-separation regime.

| Quantity | Value |
|---|---|
| Weber `We = ρ_l Ur² D / σ` | 64.9 |
| Reynolds `Re = ρ_l Ur D / μ_l` | 312.8 |
| Impact parameter `B` (offset/D) | 0.71 |
| Ohnesorge `Oh = μ_l/√(ρ_l σ D)` | 0.0258 |
| Drop radius `R` | 177 µm (D = 354 µm) |
| ρ_l, ρ_g | 763, 1.146 kg/m³ (ratio 666) |
| μ_l, μ_g | 2.183e-3, 1.834e-5 Pa·s (ratio 119) |
| σ | 0.0266 N/m |
| Relative velocity `Ur` | 2.528 m/s (each drop ±Ur/2) |
| t* = D/Ur | 0.140 ms; run to 3 ms ≈ 21 t* |

Two drops of radius R, centres offset ±0.68D along x and ±0.355D (= B·D/2)
across it, approaching at ±Ur/2. Units are SI, so the time axis (0–3 ms) and
lengths match the MFC output directly.

## Why incompressible

MFC uses an artificial-compressible 5-equation model with the liquid sound speed
knocked down to `c_l = 100 m/s` to keep the acoustic timestep affordable. The
real Mach number here is ~0.002 — the flow is incompressible — so Basilisk's
incompressible VOF is a faithful (and, on the acoustic side, cleaner) model of
the same physics. `navier-stokes/conserving.h` (momentum-conserving advection)
and `#define FILTERED` are used because the density ratio (666) is large.

## Build & run

Needs Basilisk on `PATH` (`export BASILISK=$HOME/basilisk/src`). `run.sh` sets
this up:

```bash
./run.sh build 8       # serial smoke test, dx_min ~ D/43
./run.sh run 32 10     # 32-rank MPI production run, dx_min ~ D/170
```

`MAXLEVEL` (finest AMR level) sets resolution: `dx_min = 6D / 2^level`.

| level | dx_min | ~ vs MFC |
|---|---|---|
| 8  | D/43  | coarse smoke test |
| 10 | D/170 | ≈ MFC production (D/147–155) |
| 11 | D/341 | finer than MFC production |

The MPI build passes `-D_GNU_SOURCE` — Intel's `mpicc -std=c99` otherwise hides
`madvise()`, which Basilisk's octree memory index needs. For the long production
run on the EPYC node, pin ranks to stay NUMA-local, e.g.
`I_MPI_PIN_DOMAIN=core mpirun -np 32 ...`.

**Cost note.** This is a full-3D AMR run to 21 t* — hours on many cores. The
configuration is mirror-symmetric about the z=0 plane (both drops on z=0, no
out-of-plane velocity), so simulating the half-domain z∈[0,3D] with a symmetry
wall at z=0 halves the cost exactly; do it by changing `origin()`/`size()` in z
and setting `u.n[back]=0`, at the cost of suppressing any z-asymmetric breakup.

## Outputs

- `stats.dat` — per-step `t dt cells ke area volume xmin xmax ymin ymax`
  (volume is the mass-conservation monitor; x/y bounds in units of D track how
  far the ligament stretches).
- `facets-NNNN.dat` — interface polygons every 10 µs (300 frames), for the
  Blender / isosurface pipeline in `../blender`.
- `dump-NNNN` — full-field snapshots every 0.1 ms (restart via
  `restore("dump-...")`; the run auto-restarts from a file named `restart`).
- `f.mp4`, `umag.mp4` — z=0 slice movies of the VOF field and speed.

## Notes

- qcc prints a couple of `dimensional constraints` messages at build time. These
  are Basilisk's optional units checker reacting to the SI constants; they are
  advisory and the build/results are unaffected.
- Validated (2026-07-11): serial and 4-rank MPI agree on mass to 6 digits
  (V/V0 = 0.9919 held constant), interface area and KE match analytic values,
  and the initial z=0 slice shows the two offset drops in the case-p geometry.
