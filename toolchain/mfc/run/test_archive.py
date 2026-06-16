"""
Tests for the archive module.

Covers run-directory planning (stem naming + collision tally), the
redirect that makes the run execute inside the run directory, the manifest
stamp, the tar/tar.zst pack-and-remove path, and the no-op path when
--archive is unset.
"""

import datetime
import os
import shutil
import subprocess
import tarfile
import tempfile
import types
import unittest
from contextlib import contextmanager
from unittest.mock import patch

from .. import state


def _make_fake_case(dirpath: str, params: dict = None):
    from .input import MFCInputFile

    case_py = os.path.join(dirpath, "case.py")
    with open(case_py, "w") as f:
        f.write("# fake case\n")

    return MFCInputFile(case_py, dirpath, params if params is not None else {})


def _fake_targets():
    return [types.SimpleNamespace(name="simulation")]


def _simulate_run_outputs(run_dir: str):
    """Drop the kind of files a real run would write into the run dir."""
    with open(os.path.join(run_dir, "simulation.inp"), "w") as f:
        f.write("&user_inputs /\n")
    with open(os.path.join(run_dir, "MFC.out"), "w") as f:
        f.write("log\n")
    os.makedirs(os.path.join(run_dir, "restart_data"))
    with open(os.path.join(run_dir, "restart_data", "lustre_0.dat"), "w") as f:
        f.write("data\n")


class _StateSandbox(unittest.TestCase):
    def setUp(self):
        self._saved_gARG = dict(state.gARG)
        state.gARG.update({"name": "MFC", "output_summary": None, "input": None})

    def tearDown(self):
        state.gARG.clear()
        state.gARG.update(self._saved_gARG)


@contextmanager
def _run_archive(fmt: str):
    from . import archive as archive_mod

    src = tempfile.mkdtemp()
    dest_root = tempfile.mkdtemp()
    try:
        state.gARG["archive"] = dest_root
        state.gARG["archive_format"] = fmt
        case = _make_fake_case(src)

        plan = archive_mod.plan_archive(case)
        archive_mod.prepare_run_dir(plan, case)
        _simulate_run_outputs(case.dirpath)
        archive_mod.finalize_archive(plan, case, _fake_targets())

        entries = sorted(os.listdir(dest_root))
        assert len(entries) == 1, f"expected one archive entry, got {entries}"
        yield plan, os.path.join(dest_root, entries[0])
    finally:
        shutil.rmtree(src, ignore_errors=True)
        shutil.rmtree(dest_root, ignore_errors=True)


class TestPrepareRunDir(_StateSandbox):
    def test_redirect_into_run_dir(self):
        from . import archive as archive_mod

        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as dest_root:
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "dir"
            case = _make_fake_case(src)

            plan = archive_mod.plan_archive(case)
            archive_mod.prepare_run_dir(plan, case)

            self.assertTrue(os.path.isdir(plan.run_dir))
            self.assertEqual(case.dirpath, plan.run_dir)
            self.assertEqual(case.filename, os.path.join(plan.run_dir, "case.py"))
            self.assertEqual(state.gARG["input"], case.filename)
            self.assertTrue(os.path.isfile(os.path.join(plan.run_dir, "case.py")))


class TestArchiveFormats(_StateSandbox):
    def test_format_dir(self):
        with _run_archive("dir") as (plan, path):
            self.assertTrue(os.path.isdir(path))
            self.assertEqual(path, plan.run_dir)
            self.assertTrue(os.path.isfile(os.path.join(path, "manifest.yaml")))
            self.assertTrue(os.path.isfile(os.path.join(path, "case.py")))
            self.assertTrue(os.path.isfile(os.path.join(path, "simulation.inp")))
            self.assertTrue(os.path.isdir(os.path.join(path, "restart_data")))
            self.assertTrue(os.path.isfile(os.path.join(path, "restart_data", "lustre_0.dat")))

    def test_format_tar(self):
        with _run_archive("tar") as (plan, path):
            self.assertTrue(path.endswith(".tar"))
            self.assertTrue(tarfile.is_tarfile(path))
            # The working run directory is removed once packed.
            self.assertFalse(os.path.isdir(plan.run_dir))
            with tarfile.open(path) as tf:
                names = tf.getnames()
            base = os.path.basename(path)[: -len(".tar")]
            self.assertIn(f"{base}/manifest.yaml", names)
            self.assertIn(f"{base}/case.py", names)
            self.assertIn(f"{base}/restart_data/lustre_0.dat", names)

    def test_format_tar_zst(self):
        try:
            r = subprocess.run(["tar", "--zstd", "--version"], capture_output=True, check=False, timeout=5)
            if r.returncode != 0:
                self.skipTest("tar --zstd not available")
        except (FileNotFoundError, subprocess.TimeoutExpired):
            self.skipTest("tar --zstd not available")

        with _run_archive("tar.zst") as (plan, path):
            self.assertTrue(path.endswith(".tar.zst"))
            self.assertGreater(os.path.getsize(path), 0)
            self.assertFalse(os.path.isdir(plan.run_dir))

            listing = subprocess.run(["tar", "--zstd", "-tf", path], capture_output=True, text=True, check=True)
            base = os.path.basename(path)[: -len(".tar.zst")]
            self.assertIn(f"{base}/manifest.yaml", listing.stdout)
            self.assertIn(f"{base}/case.py", listing.stdout)


class TestArchiveBehavior(_StateSandbox):
    def test_plan_returns_none_when_archive_unset(self):
        from . import archive as archive_mod

        state.gARG["archive"] = None
        state.gARG["archive_format"] = "dir"
        # No case needed since plan_archive returns before touching it.
        self.assertIsNone(archive_mod.plan_archive(case=None))

    def test_plan_bad_format_raises(self):
        from ..common import MFCException
        from . import archive as archive_mod

        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as dest_root:
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "bogus"
            case = _make_fake_case(src)
            with self.assertRaises(MFCException):
                archive_mod.plan_archive(case)

    def test_plan_uses_case_dir_name_as_stem(self):
        from . import archive as archive_mod

        fixed = datetime.datetime(2026, 1, 1, 12, 0, 0)

        with tempfile.TemporaryDirectory() as parent, tempfile.TemporaryDirectory() as dest_root, patch("mfc.run.archive.datetime.datetime") as MockDT:
            MockDT.now.return_value = fixed
            case_dir = os.path.join(parent, "my_cool_case")
            os.makedirs(case_dir)
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "dir"
            case = _make_fake_case(case_dir)

            plan = archive_mod.plan_archive(case)

        self.assertEqual(plan.stem, "my_cool_case-20260101-120000")
        self.assertTrue(plan.run_dir.endswith("my_cool_case-20260101-120000"))
        self.assertEqual(plan.dest, plan.run_dir)

    def test_plan_collision_gets_tally_suffix(self):
        from . import archive as archive_mod

        fixed = datetime.datetime(2026, 1, 1, 12, 0, 0)

        with tempfile.TemporaryDirectory() as src, tempfile.TemporaryDirectory() as dest_root, patch("mfc.run.archive.datetime.datetime") as MockDT:
            MockDT.now.return_value = fixed
            state.gARG["archive"] = dest_root
            state.gARG["archive_format"] = "dir"
            case = _make_fake_case(src)

            plan1 = archive_mod.plan_archive(case)
            os.makedirs(plan1.run_dir)  # simulate an existing run at that path
            plan2 = archive_mod.plan_archive(case)
            os.makedirs(plan2.run_dir)
            plan3 = archive_mod.plan_archive(case)

        self.assertTrue(plan2.run_dir.endswith("-2"))
        self.assertTrue(plan3.run_dir.endswith("-3"))
        self.assertNotEqual(plan1.run_dir, plan2.run_dir)
        self.assertNotEqual(plan2.run_dir, plan3.run_dir)


if __name__ == "__main__":
    unittest.main()
