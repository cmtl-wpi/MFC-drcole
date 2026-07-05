# case_s_nf3 — status notes (2026-07-05)

- **Superseded.** This `num_fluids=3` anti-merge prototype (drop A = fluid 1,
  gas = fluid 2, drop B = fluid 3, `int_comp=0`) predates the design doc
  `docs/gas_film_drainage_bouncing/DESIGN.md`. That doc's Phase 0 replaces it for
  head-on cases: symmetry-plane half-domain, `num_fluids=2`, sub-grid film model,
  no markers. Multi-marker work is deferred to the oblique S3 stage.
- **Run killed mid-flight** (2026-07-05, ~60%): with no film physics it cannot
  produce a bounce, so completing it answered nothing after the Chen & Yang (2020)
  reframing (bouncing at We=8.6 is real; all-merge results are a film-resolution
  artifact).
- **Fluid-ordering caveat:** this case puts gas as fluid 2. The design doc's
  N-fluid mTHINC closure convention requires the continuous phase LAST
  (drop A = 1, drop B = 2, gas = 3). Reorder before reusing this case for S3.
