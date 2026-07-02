#!/usr/bin/env python3
"""Extract the alpha=0.5 interface from an MFC 3D silo run -> binary PLY.

Reads the 8 per-rank Silo-HDF5 files directly (VTK has no Silo reader),
stitches the 2x2x2 decomposition into one global array via a unique
cell-center merge (handles the ghost-cell overlap at the split planes),
then marching-cubes and maps index-space verts through the stretched
cell-center coordinate arrays to physical space.
"""
import os, sys, struct, numpy as np, h5py
from skimage.measure import marching_cubes
from scipy.ndimage import gaussian_filter, map_coordinates

# Silo output dir (per-rank pK/<ts>.silo). Override with MFC_SILO=/path/to/silo_hdf5
SILO = os.environ.get("MFC_SILO",
                      "/home/daveygravy/repos/MFC/examples/3D_droplet_coalescence/silo_hdf5")
# dataset order confirmed by physical value-range: #4 vel1 #5 vel2 #6 vel3 #7 pres #8 alpha1 #9 cf
FIELD = {"alpha1": "#000008", "cf": "#000009"}

def read_rank(pk, ts, dset):
    g = h5py.File(f"{SILO}/{pk}/{ts}.silo", "r")["/.silo"]
    ks = sorted([k for k in g if k.startswith("#")], key=lambda s: int(s[1:]))
    arr = {k: g[k][()] for k in ks}
    coords = [arr[k] for k in ks if arr[k].ndim == 1]   # x,y,z nodes
    x, y, z = coords
    raw = arr[dset]                                      # Silo stores Fortran-order;
    f = raw.ravel(order="C").reshape(raw.shape, order="F")  # h5py gives it C-order -> reinterpret
    cc = lambda n: 0.5 * (n[:-1] + n[1:])
    return cc(x), cc(y), cc(z), f

def build_axis(center_arrays):
    """Union of cell-centers across ranks -> sorted global axis + key map."""
    allc = np.concatenate(center_arrays)
    keys = np.round(allc / 1e-10).astype(np.int64)      # 0.1 nm dedup tol
    ukeys, idx = np.unique(keys, return_index=True)     # sorted
    return allc[idx], {int(k): i for i, k in enumerate(ukeys)}

def key(a):
    return np.round(a / 1e-10).astype(np.int64)

def write_ply(path, V, N, F):
    with open(path, "wb") as fp:
        hdr = (f"ply\nformat binary_little_endian 1.0\n"
               f"element vertex {len(V)}\n"
               f"property float x\nproperty float y\nproperty float z\n"
               f"property float nx\nproperty float ny\nproperty float nz\n"
               f"element face {len(F)}\n"
               f"property list uchar int vertex_indices\n"
               f"end_header\n")
        fp.write(hdr.encode())
        vn = np.hstack([V, N]).astype("<f4")
        fp.write(vn.tobytes())
        cnt = np.full((len(F), 1), 3, np.uint8)
        for c, tri in zip(cnt, F.astype("<i4")):
            fp.write(c.tobytes()); fp.write(tri.tobytes())

def main():
    ts = int(sys.argv[1]); out = sys.argv[2]
    fld = FIELD[sys.argv[3] if len(sys.argv) > 3 else "alpha1"]
    sigma = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
    nranks = len([d for d in os.listdir(SILO) if d.startswith("p") and d[1:].isdigit()])
    ranks = [read_rank(f"p{i}", ts, fld) for i in range(nranks)]
    gx, mx = build_axis([r[0] for r in ranks])
    gy, my = build_axis([r[1] for r in ranks])
    gz, mz = build_axis([r[2] for r in ranks])
    G = np.zeros((len(gx), len(gy), len(gz)), np.float32)   # 0 = gas background
    for xc, yc, zc, f in ranks:
        i0 = mx[int(key(xc[0]))]; j0 = my[int(key(yc[0]))]; k0 = mz[int(key(zc[0]))]
        G[i0:i0+f.shape[0], j0:j0+f.shape[1], k0:k0+f.shape[2]] = f
    if sigma > 0:
        G = gaussian_filter(G, sigma=sigma)   # round the sharp interface for a clean render surface
    verts, faces, _, _ = marching_cubes(G, level=0.5)
    px = np.interp(verts[:, 0], np.arange(len(gx)), gx)
    py = np.interp(verts[:, 1], np.arange(len(gy)), gy)
    pz = np.interp(verts[:, 2], np.arange(len(gz)), gz)
    V = np.column_stack([px, py, pz]) * 1e4                 # meters -> ~cm-scale Blender units
    # outward, smooth normals from the volume gradient: -grad(G) points liquid->gas.
    # (skimage/geometry normals are unreliable in sign and noisy; this is bulletproof.)
    gi, gj, gk = np.gradient(G)
    smp = verts.T
    n = -np.stack([map_coordinates(gi, smp, order=1),
                   map_coordinates(gj, smp, order=1),
                   map_coordinates(gk, smp, order=1)], axis=1)
    n /= (np.linalg.norm(n, axis=1, keepdims=True) + 1e-30)
    # make face winding right-handed wrt the outward normal (consistent front faces)
    tv = V[faces]
    fn = np.cross(tv[:, 1] - tv[:, 0], tv[:, 2] - tv[:, 0])
    flip = np.einsum("ij,ij->i", fn, n[faces[:, 0]]) < 0
    faces[flip] = faces[flip][:, ::-1]
    normals = n
    write_ply(out, V, normals, faces)
    ext = V.max(0) - V.min(0)
    print(f"ts={ts} field={sys.argv[3] if len(sys.argv)>3 else 'alpha1'} "
          f"grid={G.shape} verts={len(V)} faces={len(faces)} "
          f"bbox(units)={ext.round(3)} center={(V.max(0)+V.min(0)).round(3)/2}")

if __name__ == "__main__":
    main()
