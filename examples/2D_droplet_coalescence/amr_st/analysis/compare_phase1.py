#!/usr/bin/env python3
"""
Aggregate the Phase-1 runs into a verdict + overlay figure.

The rigorous seam-current metric compares the AMR **seam-band** max|u| to the
uniform-coarse max|u| in the SAME cell-index band (control) at matched physical time --
this isolates any current the coarse/fine seam ADDS from the ambient parasitic flow that
reaches that region anyway. The containment hypothesis passes iff, over the run, that
ratio stays O(1) (<~2x) and non-growing. The PR-1628 failure is 27-540x AND growing.

Inputs (written by seam_analysis.py):
  coarse.json       uniform-coarse, domain metrics
  coarse_band.json  uniform-coarse with --band-block = AMR block  (the control)
  fine.json         uniform-fine, domain metrics
  amr.json          AMR run, domain + seam-band metrics

Usage: compare_phase1.py
"""
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXP = "/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st"
RES = os.path.join(EXP, "results")


def load(lab):
    fp = os.path.join(RES, f"{lab}.json")
    return json.load(open(fp)) if os.path.exists(fp) else None


d = {lab: load(lab) for lab in ("coarse", "coarse_band", "fine", "amr")}
present = {k: v for k, v in d.items() if v}
if not present:
    raise SystemExit("no run results found")

tau = next(iter(present.values()))["summary"]["tau"]
Us = next(iter(present.values()))["summary"]["U_sigma"]


def series(lab, key):
    s = present[lab]["series"]
    return np.array(s["t"]) / tau, np.array(s[key])


# ---- overlay figure -----------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
if d["coarse"]:
    t, v = series("coarse", "max_speed"); ax.plot(t, v / Us, "-", color="#1f77b4", label="coarse: domain max|u|")
if d["fine"]:
    t, v = series("fine", "max_speed"); ax.plot(t, v / Us, "-", color="#2ca02c", label="fine: domain max|u|")
if d["amr"]:
    t, v = series("amr", "max_speed"); ax.plot(t, v / Us, "-", color="#d62728", label="amr: domain max|u|")
    t, v = series("amr", "max_speed_seam"); ax.plot(t, v / Us, "--o", ms=3, color="black", label="amr: SEAM-band max|u|")
if d["coarse_band"]:
    t, v = series("coarse_band", "max_speed_seam"); ax.plot(t, v / Us, ":s", ms=3, color="gray", label="coarse: same-band max|u| (control)")
ax.set_xlabel("t / tau"); ax.set_ylabel("max|u| / U_sigma")
ax.set_title("Phase 1: contained-interface seam current vs uniform parasitic baseline")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fp = os.path.join(RES, "figures", "phase1_overlay.png")
os.makedirs(os.path.dirname(fp), exist_ok=True)
fig.savefig(fp, dpi=200, bbox_inches="tight"); plt.close(fig)

# ---- verdict: AMR seam vs coarse same-band, at the latest common time ----------
summary = {"runs": {lab: present[lab]["summary"] for lab in present}, "overlay": fp}
if d["amr"] and d["coarse_band"]:
    ta, sa = series("amr", "max_speed_seam")
    tc, sc = series("coarse_band", "max_speed_seam")
    t_common = min(ta.max(), tc.max())
    ia = int(np.argmin(np.abs(ta - t_common)))
    ic = int(np.argmin(np.abs(tc - t_common)))
    amr_seam = float(sa[ia]); coarse_band = float(sc[ic])
    ratio = amr_seam / coarse_band if coarse_band else float("inf")
    seam_growth = present["amr"]["summary"].get("seam_growth_ratio")
    non_growing = (seam_growth is None) or (seam_growth < 1.5)
    contained = bool(present["amr"]["summary"].get("containment_pass"))
    verdict = {
        "t_over_tau_compared": float(t_common),
        "amr_seam_maxu": amr_seam,
        "coarse_same_band_maxu": coarse_band,
        "amr_over_coarse_band": ratio,
        "within_2x": ratio < 2.0,
        "amr_seam_growth_ratio": seam_growth,
        "non_growing": non_growing,
        "containment_audit_pass": contained,
        "hypothesis_supported": bool(ratio < 2.0 and non_growing and contained),
        "amr_complete": present["amr"]["summary"]["t_final_over_tau"] >= 1.9,
        "note": ("AMR seam-band max|u| vs uniform-coarse in the same cell band; "
                 "PR-1628 failure is 27-540x AND growing"),
    }
    summary["verdict"] = verdict
    print("VERDICT " + json.dumps(verdict, indent=2))

with open(os.path.join(RES, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(f"[compare] wrote {os.path.join(RES, 'summary.json')} and {fp}")
