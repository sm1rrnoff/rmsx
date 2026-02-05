import os
import unittest
from pathlib import Path


class TestNotebookSafety(unittest.TestCase):
    def test_import_rmsx_without_core_deps(self) -> None:
        # Importing the top-level package should not eagerly import heavy deps
        # (e.g., pandas/MDAnalysis) so notebook users can at least access flipbook tools.
        import rmsx  # noqa: F401

    def test_run_flipbook_exposed_and_lightweight(self) -> None:
        from rmsx import run_flipbook

        self.assertEqual(run_flipbook.__module__, "rmsx.flipbook")

    def test_run_flipbook_raises_not_sysexit(self) -> None:
        from rmsx import run_flipbook

        with self.assertRaises(NotADirectoryError):
            run_flipbook("/this/does/not/exist", viewer="chimerax")

    def test_vmd_scripts_packaged_and_loader_is_cwd_independent(self) -> None:
        import rmsx.vmd_scripts as vmd_scripts

        scripts_dir = Path(vmd_scripts.__file__).resolve().parent
        wait_to_load = scripts_dir / "wait_to_load.tcl"
        grid_script = scripts_dir / "grid_color_scale_centered_xaxis_hotkeys.tcl"

        self.assertTrue(wait_to_load.is_file(), "wait_to_load.tcl should be present next to vmd_scripts package")
        self.assertTrue(grid_script.is_file(), "grid TCL script should be present next to vmd_scripts package")

        text = wait_to_load.read_text(encoding="utf-8")
        self.assertIn("[info script]", text)
        self.assertIn("after idle", text)
        self.assertIn("VMDMODULATENEWTUBE", text)

    def test_flipbook_can_find_loader_script_on_disk(self) -> None:
        # When installed normally (wheel unpacked), this path should exist.
        # In this repo checkout, it should also exist.
        import rmsx.flipbook as flipbook

        loader = os.path.join(os.path.dirname(flipbook.__file__), "vmd_scripts", "wait_to_load.tcl")
        self.assertTrue(os.path.exists(loader))


if __name__ == "__main__":
    unittest.main()
