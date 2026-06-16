#!/usr/bin/env bash
##
## milanx.mako — MFC run template for a single-node dual AMD EPYC Milan-X workstation.
##
## Tested on: 2× EPYC 7773X (128 physical cores, 16 CCDs/L3 domains, 2 NUMA nodes,
##            1 TiB DDR4, system OpenMPI 5, gfortran). No scheduler.
##
## Use:        ./mfc.sh run case.py -c milanx -n <ranks>
##
## Pinning strategy:
##   - `--map-by l3cache` round-robins ranks across the 16 L3 (V-Cache) domains, so
##     ranks that are MPI neighbors share a 96 MB V-Cache slice — halo exchanges hit
##     L3 instead of DRAM. `--bind-to core` then pins each rank to a single physical
##     core within its assigned L3 domain.
##   - For pure-MPI MFC, ranks = 128 (one per physical core, SMT siblings unused) is
##     the usual sweet spot. 256 ranks doubles SMT contention per core; 64 ranks
##     leaves V-Cache underused unless you raise OMP_NUM_THREADS to fill the CCD.
##
<%namespace name="helpers" file="helpers.mako"/>

<%
if engine == 'batch':
    raise Exception("milanx.mako is a single-node interactive template; batch engine is not supported.")
%>

${helpers.template_prologue()}

% if mpi:
    # Resolve mpirun (OpenMPI). Allow override via --binary.
    for binary in ${binary or ''} mpirun mpiexec; do
        if command -v "$binary" > /dev/null; then break; fi
    done
    if ! command -v "$binary" > /dev/null; then
        error ":( Could not find mpirun/mpiexec on PATH."
        exit 1
    fi
    ok ":) Using MPI launcher $MAGENTA$binary$COLOR_RESET."

    # Pure-MPI runtime defaults. MFC simulation calls OpenMP only via GPU offload
    # macros, so on CPU we want one thread per rank.
    export OMP_NUM_THREADS=${'${OMP_NUM_THREADS:-1}'}
    export OMP_PROC_BIND=${'${OMP_PROC_BIND:-close}'}
    export OMP_PLACES=${'${OMP_PLACES:-cores}'}

    # System OpenMPI (Ubuntu) disables UCX and uses the shared-memory BTL "sm"
    # with the ob1 PML, which is optimal for a single-node box. No UCX env needed.

    # Warn (don't fail) if the user oversubscribes physical cores.
    phys_cores=$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)
    total_ranks=$(( ${nodes} * ${tasks_per_node} ))
    if [ "$total_ranks" -gt "$phys_cores" ]; then
        warn ":! $total_ranks ranks > $phys_cores physical cores; SMT siblings will share L1/L2."
    fi
% endif

% for target in targets:
    ${helpers.run_prologue(target)}

    % if not mpi:
        (set -x; ${profiler} "${target.get_install_binpath(case)}")
    % else:
        # L3-cache-aware mapping is the single biggest run-time tuning knob on Milan-X.
        # If a case decomposes into chunks larger than ~80 MB per rank, this still helps;
        # if smaller, the V-Cache hit rate climbs sharply.
        (set -x; ${profiler}                                    \
            "$binary" -np ${nodes*tasks_per_node}               \
                      --map-by l3cache                          \
                      --bind-to core                            \
                      --mca pml ob1                             \
                      --mca btl self,sm                         \
                      ${'${MFC_MPIRUN_EXTRA_FLAGS:-}'}          \
                      "${target.get_install_binpath(case)}")
    % endif

    ${helpers.run_epilogue(target)}

    echo
% endfor

${helpers.template_epilogue()}
