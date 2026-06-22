---
name: tc1-droop-is-numerical-interface-diffusion
description: The TC1/Fig5 conduction late-time velocity droop is numerical color-interface smearing, not conduction strength, grid magnitude, or deformation
metadata:
  type: decision
---

The late-time decay ("droop") of the TC1/Fig 5 conduction rise velocity (v/v_YGB peaks ~1·t_r
then sags toward t/t_r=10) is a **numerical interface-diffusion artifact**, confirmed by
diagnostics on the w064 sharp (sc=0.99) Ma=0.1 run over the full 10·t_r window:

- drop **volume** conserved to +0.08% → NOT mass loss/leakage
- interface **length**/perimeter −0.8%, stays circular (L/L_circle ≈ 0.98 throughout) → NOT deformation
- interface **thickness** grows 2.24·dx → 3.11·dx (×1.39), monotonic
- **corr(v/v_YGB, band width) = −0.98** → the velocity decay IS the interface spreading

Mechanism: the color function is passively advected; WENO numerically diffuses it and MFC
(diffuse-interface CSF, no THINC/anti-diffusion/reconstruction) cannot re-sharpen it, so the
band widens ~indefinitely and the Marangoni CSF force smears over more cells and weakens.

Hypotheses RULED OUT by the run matrix (all w064 unless noted):
- conduction strength (Ma): Ma=0.1 and Ma=0.01 droop **identically** through t/t_r=6
  (droop 0.215 vs 0.171; curves overlap) → not a Marangoni-number effect. [[tc1-conduction-ma-sweep-result]]
- grid magnitude: w064→w128 droop barely changes (sharp 0.215→0.190, thick 0.277→0.272) even
  though finer grid lifts the whole plateau ~+0.06. [[thermocapillary-marangoni-sweep]]

Interface SHARPNESS is the lever that matters: droop sc0.99 < sc0.5 at both grids
(w064 0.215 vs 0.277; w128 0.190 vs 0.272). This is the smooth_coeff=0.99 study
(half-width w=dx/smooth_coeff; sc0.99→w≈dx is ~2× sharper than the sc0.5 default).

Why Samareh stays flat at ~0.83: VOF reconstructs the interface every step → stays sharp →
true steady plateau. The problem is physically steady (drop moves only ~0.5D in a wall-pinned
linear field), so any decay is numerical and vanishes in the sharp-interface/converged limit.
Refines [[frozen-t-proxy-advects-not-frozen]] (that was the frozen-T density proxy; this is the
conduction case and the mechanism is color-band smearing, isolated quantitatively).

Diagnostic recipe: read restart_data color field (c_idx = nvars-1, the last conserved variable),
per snapshot compute volume=sum(c)·dxdy, length=sum|∇c|·dxdy (co-area), thickness=band_area/length,
band_area=count(0.1<c<0.9)·dxdy. Run dirs: runs/tc1/ma0p1/w064/{sc050,sc099}, .../w128/..., full 10·t_r.
