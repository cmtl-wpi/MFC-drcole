# Blender rendering of MFC 3D droplet coalescence

Renders the liquid interface (α = 0.5 isosurface) of an MFC 3D silo run as a
Blender image. Two stages, decoupled so Blender's bundled Python never has to
read MFC data:

1. **`extract_iso.py`** (system Python: numpy, h5py, scipy, scikit-image) reads
   the per-rank Silo-HDF5 files directly, stitches the domain, marching-cubes the
   α = 0.5 interface, and writes a binary PLY.
2. **`render_blender.py`** (Blender's Python) imports the PLY and renders it with
   a studio 3-sun rig, opaque water-blue material, Cycles + denoise.

## Usage

```bash
# 1. extract a timestep's interface to PLY  (ts, out.ply, field, gaussian-sigma)
MFC_SILO=/path/to/run/silo_hdf5 python3 extract_iso.py 30 iso_30.ply alpha1 2.0

# 2. render it  (in.ply, out.png, samples, view[3q|front|side|top], style)
blender -b -P render_blender.py -- iso_30.ply out.png 200 3q opaque
```

Blender isn't a dependency of MFC — grab the portable tarball from
`download.blender.org` (no root needed) and point at `./blender/blender`.
`style` = `opaque` (reliable, recommended), `glass`, `frosted`, or `matte`.

## Non-obvious gotchas (these cost real debugging time)

- **VTK/ParaView here can't read Silo.** The system ParaView build has no Silo
  reader, and VTK/pyvista never do. Reading the Silo-HDF5 with `h5py` directly is
  the zero-dependency path. The per-rank datasets are `#000001..3` = x/y/z node
  coords, `#000004..9` = the 6 wrt fields. Identify fields by physical range
  (α, cf ∈ [0,1]; pres ~1e3+; vel signed O(1)) rather than trusting write order.
- **Silo stores Fortran (column-major) order; h5py reads it C-major** → the two
  trailing axes come out scrambled (looks like horizontal striping / a droplet
  smeared across the whole domain). Fix: `buf.ravel('C').reshape(shape, order='F')`.
- **Sharp interface → marching-cubes terracing.** MFC's interface is ~1 cell
  wide by design; contouring it raw gives a stairstepped surface. A small
  `gaussian_filter(sigma≈2)` on the volume before contouring cleans it up.
- **Isosurface normals must be set explicitly and outward.** skimage's normals
  and Blender's face-averaged normals are both noisy (moiré) and can point
  inward (renders pure black). Compute smooth outward normals from `-∇G` sampled
  at each vertex — this fixes both the moiré and the black-object-in-bright-scene.
- The run on disk is a 2×2×2 rank decomposition with a ~4-cell ghost overlap at
  each split plane; `extract_iso.py` merges it by unique cell-center. Rank count
  is auto-detected from the `pK/` dirs.

## Frames

`frames/` holds a 4-step story of the head-on collision:
`00_spheres` (approach) → `15_merge` (coalescence bridge) → `30_pancake`
(biconcave impact disc) → `90_capsule` (retracted, oscillating).
