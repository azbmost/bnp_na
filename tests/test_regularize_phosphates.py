from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from bnp_na_lib.regularize_phosphates import (
    RegularizePhosphatesError,
    default_regularized_output_path,
    regularize_phosphates,
)


def atom_line(serial: int, atom: str, chain: str, resseq: int, coord: np.ndarray) -> str:
    chars = list(" " * 80)
    chars[0:6] = "ATOM  "
    chars[6:11] = f"{serial:5d}"
    chars[12:16] = f"{atom:>4s}"
    chars[17:20] = " DA"
    chars[21] = chain
    chars[22:26] = f"{resseq:4d}"
    chars[30:38] = f"{coord[0]:8.3f}"
    chars[38:46] = f"{coord[1]:8.3f}"
    chars[46:54] = f"{coord[2]:8.3f}"
    chars[54:60] = f"{1.0:6.2f}"
    chars[60:66] = f"{0.0:6.2f}"
    chars[76:78] = " P" if atom == "P" else " O" if atom.startswith("O") else " C"
    return "".join(chars) + "\n"


def step_coord(coord: np.ndarray, index: int) -> np.ndarray:
    result = np.asarray(coord, dtype=float)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    for _ in range(index):
        result = rotation @ result + np.array([0.0, 0.0, 2.0])
    return result


def read_atom_coords(path: Path) -> dict[tuple[int, str], np.ndarray]:
    coords = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        coords[(int(line[22:26]), line[12:16].strip())] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    return coords


class RegularizePhosphatesTests(unittest.TestCase):
    def _write_helix(self, path: Path, *, include_three_prime: bool = False, count: int = 6) -> None:
        templates = {
            "P": np.array([7.0, 1.0, -0.5]),
            "OP1": np.array([8.0, 1.0, -0.5]),
            "OP2": np.array([7.0, 2.0, -0.5]),
            "O5'": np.array([6.5, 0.5, -0.5]),
            "C5'": np.array([6.0, 0.0, -0.25]),
            "O3'": np.array([5.5, -0.5, 0.25]),
        }
        internal_noise = [
            np.array([0.10, 0.00, 0.00]),
            np.array([-0.10, 0.00, 0.00]),
            np.array([0.00, 0.10, 0.00]),
            np.array([0.00, -0.10, 0.00]),
            np.array([0.00, 0.00, 0.00]),
        ]
        lines = ["REMARK ORIGINAL\n"]
        serial = 1
        for index in range(count):
            for atom_name, template in templates.items():
                if atom_name == "O3'":
                    local_noise = (
                        np.array([-3.0, 2.0, -1.0])
                        if index == count - 1
                        else internal_noise[index]
                    )
                else:
                    local_noise = (
                        np.array([3.0, -2.0, 1.0])
                        if index == 0
                        else internal_noise[index - 1]
                    )
                lines.append(atom_line(serial, atom_name, "A", index + 1, step_coord(template + local_noise, index)))
                serial += 1
            lines.append(atom_line(serial, "C1'", "A", index + 1, step_coord(np.array([5.0, 0.0, 0.0]), index)))
            serial += 1
            lines.append(atom_line(serial, "C2'", "A", index + 1, step_coord(np.array([4.5, 0.5, 0.0]), index)))
            serial += 1
        if include_three_prime:
            index = count
            for atom_name in ("P", "OP1", "OP2"):
                template = templates[atom_name]
                terminal_noise = np.array([-4.0, 3.0, -1.5])
                lines.append(atom_line(serial, atom_name, "A", index + 1, step_coord(template + terminal_noise, index)))
                serial += 1
        lines.append("END\n")
        path.write_text("".join(lines), encoding="utf-8")

    def test_terminal_group_is_excluded_then_regularized(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.pdb"
            output_path = Path(directory) / "output.pdb"
            self._write_helix(input_path)
            original = read_atom_coords(input_path)

            result = regularize_phosphates(input_path, output_path)
            coords = read_atom_coords(output_path)

            templates = {
                "P": np.array([7.0, 1.0, -0.5]),
                "OP1": np.array([8.0, 1.0, -0.5]),
                "OP2": np.array([7.0, 2.0, -0.5]),
                "O5'": np.array([6.5, 0.5, -0.5]),
                "C5'": np.array([6.0, 0.0, -0.25]),
                "O3'": np.array([5.5, -0.5, 0.25]),
            }
            for index in range(6):
                for atom_name, template in templates.items():
                    np.testing.assert_allclose(coords[(index + 1, atom_name)], step_coord(template, index), atol=0.001)
            np.testing.assert_array_equal(coords[(1, "C2'")], original[(1, "C2'")])
            self.assertEqual(result.chains[0].terminal_group_count, 1)
            self.assertEqual(result.chains[0].regularized_atom_count, 36)
            self.assertIn("TERMINAL ATOMS EXCLUDED FROM CONSENSUS", output_path.read_text(encoding="utf-8"))

            linkage_distances = []
            for resseq in range(2, 7):
                linkage_distances.append(float(np.linalg.norm(coords[(resseq, "P")] - coords[(resseq - 1, "O3'")])))
            np.testing.assert_allclose(linkage_distances, [linkage_distances[0]] * 5, atol=0.002)
            for resseq in range(1, 7):
                self.assertAlmostEqual(
                    float(np.linalg.norm(coords[(resseq, "P")] - coords[(resseq, "O5'")])),
                    float(np.linalg.norm(templates["P"] - templates["O5'"])),
                    delta=0.002,
                )

    def test_three_prime_phosphate_only_residue_is_regularized_after_fit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.pdb"
            output_path = Path(directory) / "output.pdb"
            self._write_helix(input_path, include_three_prime=True)

            result = regularize_phosphates(input_path, output_path)
            coords = read_atom_coords(output_path)

            self.assertEqual(result.chains[0].terminal_group_count, 2)
            np.testing.assert_allclose(
                coords[(7, "P")],
                step_coord(np.array([7.0, 1.0, -0.5]), 6),
                atol=0.001,
            )

    def test_requires_four_c1_atoms(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "short.pdb"
            self._write_helix(input_path, count=3)
            with self.assertRaisesRegex(RegularizePhosphatesError, "at least 4"):
                regularize_phosphates(input_path)

    def test_default_output_name(self) -> None:
        self.assertEqual(
            default_regularized_output_path(Path("model.pdb")),
            Path("model_regularized_phosphates.pdb"),
        )


if __name__ == "__main__":
    unittest.main()
