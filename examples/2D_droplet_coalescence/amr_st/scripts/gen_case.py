#!/usr/bin/env python3
"""Write a standalone case.py into a run dir with the case-args baked into sys.argv.

Sidesteps mfc.sh's `-- <args>` forwarding entirely (mfc passes only --mfc, which the
wrapper ignores). The restart index sits on its own `_NSTART` line for the resume loop.

  gen_case.py <run_id> <run_dir> -- <case args for case_laplace.py ...>
"""
import os
import sys

CASE = "/home/daveygravy/repos/MFC/examples/2D_droplet_coalescence/amr_st/cases/case_laplace.py"

run_id, run_dir = sys.argv[1], sys.argv[2]
args = sys.argv[3:]
if args and args[0] == "--":
    args = args[1:]

os.makedirs(run_dir, exist_ok=True)
argv_list = ", ".join(f'"{a}"' for a in args)
with open(os.path.join(run_dir, "case.py"), "w") as f:
    f.write(
        f"import sys as _s   # AUTO-GENERATED standalone case ({run_id}); mfc args ignored\n"
        f"_NSTART = 0\n"
        f'_s.argv = [_s.argv[0], {argv_list}, "--n-start", str(_NSTART), "--mfc", "{{}}"]\n'
        f'exec(open("{CASE}").read())\n'
    )
print(f"wrote {run_dir}/case.py  (baked: {' '.join(args)})")
