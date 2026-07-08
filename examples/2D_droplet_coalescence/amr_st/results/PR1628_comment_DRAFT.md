# DRAFT — PR 1628 feedback comment (DO NOT POST without Davis's approval)

Status: draft. Numbers are from a reduced np=1 config against pinned SHA `ace2285a`; the
containment result is preliminary (AMR run mid-flight). Finalize the AMR numbers and get
approval before posting.

---

I lifted the `surface_tension`-under-AMR gate at `ace2285a` to characterize the
**contained-interface** case (interface kept well inside the static block, away from the
2:1 seam) — the case §7.1's three attempted fixes didn't cover (all three had the
interface crossing/near the seam). Two things came out of it.

**1. Lifting the gate exposes a surface-tension-specific bug in the parallel fine-halo
exchange — the gate is hiding a broken path, not just a caution.**

On a static-drop Laplace case (5-eq, `mpp_lim`, WENO5, RK3, one drop + static 2:1 block):

| config | surface_tension | result |
|---|---|---|
| AMR, np=1 | on | runs, stable |
| AMR, np=2 | on | **hangs at step 0** (fine-halo deadlock) |
| AMR, np=16 | on | **`MPI_Sendrecv` message-truncated abort** at step 0 |
| AMR, np=16 | **off** | runs to completion (control) |

The identical case runs fine at np=16 **without** surface_tension, so this isn't a general
decomposition bug — enabling the color function (which adds a variable to `sys_size`)
breaks `s_mpi_sendrecv_amr_fine_halo`. The AMR core is unaffected: prolong / restrict /
reflux conservation all report ~1e-16 at every config. The `MPI_SENDRECV`
(`m_mpi_proxy.fpp:181`) uses one `cnt` for both send and recv, assuming the neighbor's
matching fine face is the same size and that the neighbor participates — under the mirror
decomposition adjacent ranks can own different-sized block pieces, and a rank owning no
block cells returns early (`:92`), so a poster can lack a partner (deadlock) or a size can
mismatch (truncation). This needs fixing before AMR+ST is usable at scale regardless of
the gate.

**2. With that worked around (np=1), a contained interface introduces no seam current.**

Interface at radius R, static block edge at 1.5R (25 coarse cells of clearance ≫ the WENO
stencil + interface width), so the color function is exactly {0,1} across the seam
(measured: cf within 2.7e-11 of {0,1} in the seam band). The seam-band max|u| tracks the
uniform-coarse parasitic flow **in the same region** — ratio ≈ 0.7× (the fine block
actually resolves the interface better, so *less* parasitic flow reaches the seam),
non-growing over the run. That's the parasitic baseline, not the 27–540× growing seam
current of §7.1 — consistent with the hypothesis that the inconsistent-normal mechanism
needs `|∇c| > capillary_cutoff` within a stencil of the seam, which containment removes by
construction.

**Suggestion:** once the fine-halo bug (1) is fixed, consider relaxing the hard ST gate to
a **containment guard** — permit `surface_tension` under AMR when the tagged/blocked region
keeps the diffuse interface ≥ (stencil + margin) cells inside every block edge, and abort
otherwise. The PR already ships this contain-and-guard pattern for EL bubble clouds and
moving IBs, so it'd be consistent. Happy to share the test case and harness.

*(Reduced config for feasibility since np=1 is the only working option today: c_l lowered,
domain ±3.2R, 2τ — none of which affect the seam mechanism.)*
