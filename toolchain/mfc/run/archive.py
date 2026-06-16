import dataclasses
import datetime
import os
import shutil
import sys
import tarfile

from ..common import MFCException, does_command_exist, file_dump_yaml, generate_git_tagline, system
from ..printer import cons
from ..state import ARG, CFG, gARG
from . import input


@dataclasses.dataclass
class ArchivePlan:
    run_dir: str  # directory the simulation executes in (inputs + outputs land here)
    dest: str  # final deliverable: == run_dir for 'dir', tarball path otherwise
    archive_format: str
    stem: str


def plan_archive(case):
    """Validate --archive settings and reserve a unique run directory.

    Runs before the simulation executes so bad paths, bad formats, or
    unwritable roots fail fast. Returns None if --archive is unset.

    The simulation executes *inside* `<archive_root>/<stem>/` (the run
    directory), so its input namelists and outputs are written there
    directly rather than next to case.py. The stem is
    `<case_dir_name>-<timestamp>` so archives from different cases dropped
    in the same archive root are self-identifying (e.g.
    1D_sodshocktube-20260424-123045) rather than all sharing the generic
    --name default.

    For `--archive-format dir` the run directory is the final artifact; for
    `tar`/`tar.zst` it is packed into a sibling tarball afterward.

    If the computed path already exists, appends "-2", "-3", ... to the
    stem until a free name is found, so two runs starting in the same
    second never collide.
    """
    archive_root = ARG("archive")
    if archive_root is None:
        return None

    archive_format = ARG("archive_format") or "dir"
    suffix_map = {"dir": "", "tar": ".tar", "tar.zst": ".tar.zst"}
    if archive_format not in suffix_map:
        raise MFCException(f"Archive: unsupported format '{archive_format}'. Must be one of: {', '.join(suffix_map)}.")
    suffix = suffix_map[archive_format]

    # Derive the stem from the case's parent directory name so archives
    # from different cases are distinguishable. Fall back to "case" if
    # the case somehow lives at the filesystem root.
    case_dir_name = os.path.basename(os.path.abspath(case.dirpath).rstrip(os.sep)) or "case"
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    base_stem = f"{case_dir_name}-{timestamp}"

    archive_root = os.path.abspath(os.path.expanduser(archive_root))
    try:
        os.makedirs(archive_root, exist_ok=True)
    except OSError as e:
        raise MFCException(f"Archive: cannot create archive root {archive_root}: {e}") from e

    def paths_for(stem: str):
        run_dir = os.path.join(archive_root, stem)
        dest = run_dir if archive_format == "dir" else os.path.join(archive_root, stem + suffix)
        return run_dir, dest

    stem = base_stem
    run_dir, dest = paths_for(stem)
    counter = 2
    while os.path.exists(run_dir) or os.path.exists(dest):
        stem = f"{base_stem}-{counter}"
        run_dir, dest = paths_for(stem)
        counter += 1

    if stem != base_stem:
        cons.print(f"[yellow]Archive: destination existed, using {stem} to avoid collision.[/yellow]")

    return ArchivePlan(run_dir=run_dir, dest=dest, archive_format=archive_format, stem=stem)


def prepare_run_dir(plan: "ArchivePlan", case: input.MFCInputFile) -> None:
    """Create the run directory and redirect the run into it.

    The simulation's input namelists, job script, working directory, and
    outputs are all keyed off `case.dirpath` / `ARG("input")`. Repointing
    both at the run directory makes the whole run execute there, so results
    are written directly into the archive rather than next to case.py.

    The case file is copied in for provenance only; it is NOT re-executed
    (its parameters were already read at load time), so analytic-IC and
    relative-path behavior in case.py is unaffected. Call this after the
    binary has been built.
    """
    run_dir = plan.run_dir
    os.makedirs(run_dir, exist_ok=True)

    case_basename = os.path.basename(case.filename)
    case_copy = os.path.join(run_dir, case_basename)
    if os.path.isfile(case.filename) and os.path.abspath(case.filename) != os.path.abspath(case_copy):
        shutil.copy2(case.filename, case_copy)

    # A chemistry case may reference its mechanism file relative to the case
    # directory; copy it alongside so it still resolves after the redirect.
    if case.params.get("chemistry", "F") == "T":
        cantera_file = case.params.get("cantera_file")
        if cantera_file:
            local_mech = os.path.join(case.dirpath, cantera_file)
            if os.path.isfile(local_mech):
                shutil.copy2(local_mech, os.path.join(run_dir, os.path.basename(cantera_file)))

    case.dirpath = run_dir
    case.filename = case_copy
    gARG["input"] = case_copy


def __build_manifest(case: input.MFCInputFile, targets, plan: "ArchivePlan") -> dict:
    run_dir = plan.run_dir
    contents = sorted(name for name in os.listdir(run_dir) if name != "manifest.yaml")

    return {
        "timestamp": datetime.datetime.now().astimezone().isoformat(),
        "source_case": os.path.abspath(case.filename),
        "source_dir": os.path.abspath(run_dir),
        "invocation": sys.argv[1:],
        "git": generate_git_tagline(),
        "targets": [t.name for t in targets],
        "archive_format": plan.archive_format,
        "archive_path": plan.dest,
        "build_lock": dataclasses.asdict(CFG()),
        "contents": contents,
    }


def __pack_tar(run_dir: str, dest: str, compressed: bool) -> None:
    """Pack the whole run directory into a tarball rooted at its basename."""
    parent = os.path.dirname(run_dir)
    arcroot = os.path.basename(run_dir)

    if compressed:
        if not does_command_exist("tar"):
            raise MFCException("Archive: 'tar' binary not found; required for --archive-format tar.zst.")

        result = system(["tar", "--zstd", "-cf", dest, "-C", parent, arcroot], print_cmd=False)
        if result.returncode != 0:
            raise MFCException(f"Archive: 'tar --zstd' failed with exit code {result.returncode}. Ensure GNU tar >= 1.31.")
        return

    with tarfile.open(dest, "w") as tf:
        tf.add(run_dir, arcname=arcroot)


def finalize_archive(plan: "ArchivePlan", case: input.MFCInputFile, targets) -> None:
    """Stamp the manifest and, for tar formats, pack the run directory.

    Caller must have run the simulation after prepare_run_dir(), so the run
    directory already holds the case file, input namelists, and outputs.
    """
    run_dir = plan.run_dir
    if not os.path.isdir(run_dir):
        cons.print(f"[yellow]Archive: run directory {run_dir} is missing; skipping.[/yellow]")
        return

    manifest = __build_manifest(case, targets, plan)
    file_dump_yaml(os.path.join(run_dir, "manifest.yaml"), manifest)

    cons.print()
    if plan.archive_format == "dir":
        cons.print(f"[bold]Archived[/bold] run at [magenta]{run_dir}[/magenta] (+ manifest.yaml).")
        return

    cons.print(f"[bold]Archiving[/bold] to [magenta]{plan.dest}[/magenta] ({plan.archive_format})")
    cons.indent()
    try:
        __pack_tar(run_dir, plan.dest, compressed=(plan.archive_format == "tar.zst"))
        shutil.rmtree(run_dir)
        cons.print("Packed the run directory into the tarball and removed the working copy.")
    finally:
        cons.unindent()
