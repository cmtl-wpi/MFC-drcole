---
name: gpu-build-on-nighthawk
description: How to build/run MFC with --gpu acc on the "nighthawk" workstation (single V100, NVHPC 24.11, no modules)
metadata:
  type: reference
---

`nighthawk` is a bare workstation (not an HPC cluster — `LMOD_SYSHOST` empty, no `module`),
one **Tesla V100-PCIE-16GB**, NVHPC SDK at `/opt/nvidia/hpc_sdk/Linux_x86_64/24.11`. There is
no `-c <slug>` in `toolchain/modules` for it, so `source ./mfc.sh load` does not apply; set the
toolchain by hand. Gotchas hit in order while getting `./mfc.sh build --gpu acc` to work:

1. **Broken venv**: `build/venv` had no pip and system `python3` lacks `ensurepip`
   (`python3-venv` not installed, no sudo). Fix: `rm -rf build/venv &&
   python3 -m venv --without-pip build/venv && build/venv/bin/python3 build/get-pip.py`
   (get-pip.py from https://bootstrap.pypa.io/pip/get-pip.py). mfc.sh then reuses it.
2. **CMake picks gfortran** → "GPU not compatible with GNU". Must `export CC=nvc CXX=nvc++ FC=nvfortran`.
3. **No MPI for nvfortran**: system OpenMPI is gfortran-built. Use NVHPC's bundled MPI
   (`comm_libs/mpi/bin` → hpcx). Put it first on PATH; set `OMPI_CC/CXX/FC`.
4. **Runtime missing libs** (`libnvToolsExt.so.1`, `libmpi.so.40`): add CUDA + hpcx-ompi libs
   to `LD_LIBRARY_PATH`.

Working env script (sourced before build AND every run):
```
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11
export PATH=$NVHPC/comm_libs/mpi/bin:$NVHPC/compilers/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC/comm_libs/12.6/hpcx/hpcx-2.20/ompi/lib:$NVHPC/compilers/lib:$NVHPC/math_libs/lib64:$NVHPC/cuda/12.6/targets/x86_64-linux/lib:$NVHPC/cuda/12.6/lib64:$NVHPC/math_libs/12.6/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export CC=nvc CXX=nvc++ FC=nvfortran OMPI_CC=nvc OMPI_CXX=nvc++ OMPI_FC=nvfortran
```
Run with `OMPI_MCA_hwloc_base_binding_policy=none` and (for n>1 on the single GPU)
`OMPI_MCA_rmaps_base_oversubscribe=1`. **Do not** use `./mfc.sh run --no-build`: cases with
analytic ICs (the thermocapillary cases) compile their IC into the binary, so each case needs
its own build; `--no-build` looks up a stale install hash and fails. Only one V100, so all MPI
ranks share device 0 (`acc_set_device_num(mod(rank,1))`).

5. **Always pass `--no-debug` for GPU runs.** Unspecified config flags inherit the last build's
   `build/lock.yaml` (via `argparse_gen.py` `set_defaults`), and the lock here is often
   `debug:true` from a prior CPU debug build. `./mfc.sh run … --gpu acc` (no `--no-debug`) then
   becomes a Debug-GPU build that maps to a half-configured staging dir and dies with
   `gmake: Makefile: No such file or directory` / "Failed to build the syscheck target".
   `--no-debug` selects the good Release-GPU config. The trailing "Terminated"/exit-143 line on
   such a failure is just mfc's outer wrapper — read higher in the log for the real error.
6. **Run the example `validate.py` scripts with `build/venv/bin/python3`**, not system `python3`
   (system python has no numpy/matplotlib; the MFC venv has numpy 2.2.6 + matplotlib 3.10.9).
7. **Run `./mfc.sh test` (and `--generate`) with `--gpu acc --no-debug`**, the same cached config
   as builds. A fresh CPU config can't build here: `--no-gpu` (MPI) fails CMake with "Could NOT
   find MPI (missing MPI_Fortran)" (no nvfortran MPI wrapper on PATH in a plain shell), and
   `--no-gpu --no-mpi` configures with system gfortran but then fails to *link* — the prebuilt
   `build/install/fftw/libfftw3.a` is non-PIE (`relocation R_X86_64_32S ... can not be used when
   making a PIE object`). The cached gpu-acc config builds via absolute compiler paths in the
   CMake cache, so even a bare `./mfc.sh build` (no env sourced) recompiles it. Goldens generated
   on the V100 pass the default 1e-12 tolerance (these short test cases are cross-platform
   reproducible — the whole CI suite compares one golden across CPU and GPU).

Pushing from this box: HTTPS remote, no system creds. `gh` 2.95.0 is installed at `~/.local/bin/gh`
(logged in as `drcole17`, token in `~/.config/gh/hosts.yml`, `gh auth setup-git` done), so `git push`
works once `~/.local/bin` is on PATH.

Related: [[thermocap-nan-pressure-relaxation]] (the case's NaN was a 6-eq pressure-relaxation
0/0, now fixed — not a GPU or build issue).
