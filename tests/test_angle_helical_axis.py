from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from bnp_na_lib.angle_helical_axisV2_2 import (
    Atom,
    BILD_SAVE_FILETYPES,
    PDB_OPEN_FILETYPES,
    RunConfig,
    _file_dialog_options,
    build_arg_parser,
    fit_axis_pca,
    parse_axis_ranges,
    run,
)


class AngleHelicalAxisGuiTests(unittest.TestCase):
    def test_pdb_dialog_patterns_are_separate_tcl_list_items(self) -> None:
        label, patterns = PDB_OPEN_FILETYPES[0]

        self.assertEqual(label, "PDB files")
        self.assertIsInstance(patterns, tuple)
        self.assertEqual(patterns, ("*.pdb", "*.ent", "*.pdb1", "*.pdb.txt", "*.txt"))
        self.assertTrue(all(" " not in pattern for pattern in patterns))

    def test_macos_dialog_omits_native_filetype_filter(self) -> None:
        options = _file_dialog_options(
            "Select PDB file",
            PDB_OPEN_FILETYPES,
            platform="darwin",
        )

        self.assertEqual(options, {"title": "Select PDB file"})

    def test_macos_save_dialog_keeps_default_extension(self) -> None:
        options = _file_dialog_options(
            "Save BILD file",
            BILD_SAVE_FILETYPES,
            defaultextension=".bild",
            platform="darwin",
        )

        self.assertEqual(
            options,
            {"title": "Save BILD file", "defaultextension": ".bild"},
        )

    def test_other_platforms_keep_filetype_filter(self) -> None:
        options = _file_dialog_options(
            "Select PDB file",
            PDB_OPEN_FILETYPES,
            platform="linux",
        )

        self.assertIs(options["filetypes"], PDB_OPEN_FILETYPES)


class AxisRangeTests(unittest.TestCase):
    @staticmethod
    def _axis_atoms() -> list[Atom]:
        atoms = []
        serial = 1
        for chain, residues in (("A", (1, 2, 3)), ("B", (6, 5, 4))):
            for index, resseq in enumerate(residues):
                atoms.append(
                    Atom(serial, "C1'", "DA", chain, resseq, 0.0, 0.0, float(index))
                )
                serial += 1
        atoms.append(Atom(serial, "C1'", "DA", "A", 99, 100.0, 100.0, 100.0))
        return atoms

    @staticmethod
    def _pdb_atom_line(atom: Atom) -> str:
        chars = list(" " * 80)
        chars[0:6] = "ATOM  "
        chars[6:11] = f"{atom.serial:5d}"
        chars[12:16] = f"{atom.name:^4s}"
        chars[17:20] = f"{atom.resName:>3s}"
        chars[21] = atom.chainID
        chars[22:26] = f"{atom.resSeq:4d}"
        chars[30:38] = f"{atom.x:8.3f}"
        chars[38:46] = f"{atom.y:8.3f}"
        chars[46:54] = f"{atom.z:8.3f}"
        return "".join(chars) + "\n"

    def test_parses_multiple_ascending_and_descending_ranges(self) -> None:
        ranges = parse_axis_ranges("A1-A3,B60-B58")

        self.assertEqual(
            [[(ref.chainID, ref.resSeq) for ref in residue_range] for residue_range in ranges],
            [[("A", 1), ("A", 2), ("A", 3)], [("B", 60), ("B", 59), ("B", 58)]],
        )

    def test_first_range_order_sets_positive_axis_direction(self) -> None:
        atoms = self._axis_atoms()
        forward = fit_axis_pca(
            atoms,
            ["C1'"],
            residue_ranges=parse_axis_ranges("A1-A3,B6-B4"),
        )
        reverse = fit_axis_pca(
            atoms,
            ["C1'"],
            residue_ranges=parse_axis_ranges("A3-A1,B4-B6"),
        )

        self.assertGreater(float(np.dot(forward.direction, np.array([0.0, 0.0, 1.0]))), 0.99)
        self.assertLess(float(np.dot(reverse.direction, np.array([0.0, 0.0, 1.0]))), -0.99)
        self.assertTrue(np.allclose(forward.point, np.array([0.0, 0.0, 1.0])))
        self.assertTrue(np.allclose(reverse.point, np.array([0.0, 0.0, 1.0])))

    def test_axis_range_excludes_atoms_outside_selected_residues(self) -> None:
        axis = fit_axis_pca(
            self._axis_atoms(),
            ["C1'"],
            residue_ranges=parse_axis_ranges("A1-A3,B6-B4"),
        )

        self.assertTrue(np.allclose(axis.point, np.array([0.0, 0.0, 1.0])))

    def test_cli_accepts_underscore_and_hyphen_spellings(self) -> None:
        parser = build_arg_parser()

        underscore = parser.parse_args(["--axis_range", "A1-A3,B6-B4"])
        hyphen = parser.parse_args(["--axis-range", "A3-A1"])

        self.assertEqual(underscore.axis_range, "A1-A3,B6-B4")
        self.assertEqual(hyphen.axis_range, "A3-A1")

    def test_run_reports_axis_range_and_writes_bild(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / "axis.pdb"
            bild = root / "axis.bild"
            pdb.write_text(
                "".join(self._pdb_atom_line(atom) for atom in self._axis_atoms()),
                encoding="utf-8",
            )

            report = run(
                RunConfig(
                    pdb=pdb,
                    axis_atoms=["C1'"],
                    axis_range="A1-A3,B6-B4",
                    point1="1 0 0",
                    point2="0 1 0",
                    out_bild=bild,
                )
            )

            self.assertIn("Axis residue ranges: A1-A3,B6-B4", report)
            self.assertIn("start-to-end order of the first axis range", report)
            self.assertIn("signed around axis (point1->point2): 90.000 deg", report)
            self.assertTrue(bild.is_file())
            self.assertIn(".arrow", bild.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
