---
name: confinement-to-one
description: "Recovering u_YGB → 1.0 is a GRID-convergence claim, not confinement: the force-free drop makes confinement negligible (W=10 grid-converges to 1.0 = unbounded). 0.95 is NOT a limit — it's Samareh's finite-resolution value, reproduced only as a matched-grid point."
metadata:
  type: project
---

`case_ygb.py` validates MFC's σ(T) Marangoni physics by recovering the Young-Goldstein-Block
analytic terminal velocity `v_t/u_YGB → 1.0`. This is a **convergence claim**, but the real axis
is **GRID**, not confinement. (Corrects the prior framing of this note.)

**The data (`results/ygb_summary.json`, 2026-06):** at fixed box W=10, the cube grid sweep
Nx 64→128 (6.4→12.8 cells/D, Ma=0.5, 3 t_r) gives plateau ratio 0.837→0.863→0.888→0.926; the
dx→0 (Richardson-in-dx) extrapolation is **1.009 ≈ 1.0**. That IS the validation.

- **Confinement is negligible here — do NOT sweep it.** The drop is neutrally buoyant
  (ρ_d=ρ_b=0.2) and equal-conductivity (k_d=k_b), so it is **force-free**: no long-range Stokeslet,
  no thermal perturbation → wall corrections scale like (a/L)³, not (a/L). Proof from the data
  itself: a *finite* W=10 box already recovers the *unbounded* YGB value (1.009). If walls mattered
  at this size you could not hit 1.0 in a finite box. The W∈{6,8,10,12} confinement sweep in
  `run_ygb.py` was **never actually run** and isn't needed.
- **Ma is nearly flat.** At Nx=80: ratio 0.858/0.863/0.869 for Ma 0.25/0.5/1.0 — ~1% spread, so
  Ma→0 extrapolation is moot. Pick a *stable* Ma: Ma=0.1 is UNSTABLE (blew up to −25 at Nx=80).
  Use Ma 0.5–1.0.
- **Grid is the only lever**, and it converges to 1.0.

**0.95 is NOT a convergence limit.** Grid refinement climbs toward 1.0 (=YGB) and sails *past* 0.95
(it crosses ≈0.95 only well beyond Nx=128, ~Nx≈190 by the linear-in-dx fit). Samareh's Fig 6 0.95
is his *finite-resolution* DNS value — ≈5% below the true YGB answer because his run wasn't fully
grid-converged. MFC's converged answer (1.0) is *more* accurate than Samareh's; "converging to
0.95" would mean deliberately under-resolving. You can only hit 0.95 as a **point** at
Samareh-matched resolution.

**The `samareh` anchor** (offset 5D×5D×7.5D box, drop 1.5D above the cold floor) reproduces
Samareh's confined Fig 6 number as that point. The `anchor` selector runs Nx=64 (= 12.8 cells/D,
the SAME resolution as cube Nx=128 → expect ~0.92–0.93 vs Samareh 0.95), Ma=1.0, 2 t_r. Launched
2026-06-23 on EPYC (~5.5 h, 8 ranks; grid is rank-limited by the weno5 ≥25-cell/block floor).
Chosen validation framing = "Both": grid-convergence → 1.0 as the headline + this single anchor as
the literature point.

**Cost:** 3D σ(T)+conduction+WENO5 is acoustic-dt-limited and multi-hour per run; small grids are
rank-limited (`ranks_for` cubes the per-dim split keeping blocks ≥~26 cells). See
[[3d-ygb-validation-harness]], [[mfc-sigmaT-3d-drift-no-conduction]], and the EPYC build gotcha
[[mfc-intel-oneapi-build-pollution]].
