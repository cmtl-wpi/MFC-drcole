---
name: working-mode-understanding-over-results
description: On this thermocapillary validation, the user is exploring to understand the physics/numerics, not to ship a figure; explain mechanisms, treat runs as experiments.
metadata:
  type: feedback
---

The user said (2026-06-23), mid-way through the Fig 5 grid-convergence work: *"i am trying to
understand things here, not necessarily just achieve a result."*

**Why:** I had slipped into ship-the-deliverable mode — optimizing the figure, recommending we
cancel "messy" runs (frozen w256) for a cleaner result, and pushing the user for keep/cancel
decisions. That optimizes for an artifact; the user is here to understand the physics and numerics.

**How to apply:**
- Lead with *mechanism*, not output. When the data looks "messy" (e.g. the frozen-T curve ringing
  harder at finer grid), explain why — don't propose hiding it. The messy part is often the
  instructive part.
- Treat runs as *experiments that test a hypothesis*, not inputs to a figure. A run that "muddies
  the figure" can still be worth doing if it answers an understanding question (frozen w256 tests
  whether the droop is grid-independent and whether acoustic ringing grows with resolution).
- Don't rush the user toward decisions or deliverables. Offer explanations and let them direct where
  to dig.

Concrete worked example from this session — the Fig 5 droop has two separable causes, which is *why*
frozen and conduction grid-converge oppositely: (1) numerical interface smearing (vanishes as dx→0,
so conduction climbs toward Samareh on refinement); (2) physical proxy-advection of the imposed
gradient (grid-independent, dominates frozen-T); plus (3) acoustic ringing that *grows* with grid
because scheme dissipation ∝ dx. See [[tc1-droop-is-numerical-interface-diffusion]] and
[[temperature-via-density-proxy]].
