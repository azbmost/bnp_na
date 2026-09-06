from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from bnp_na_lib.opposing_phosphate_xdisp import (
    OpposingPhosphateCancelled,
    OpposingPhosphateError,
    find_opposing_phosphate_xdisp,
    format_result_report,
    parse_p_atoms,
    phosphate_pair_rows,
    quantize_x,
    summarize_rows,
    validate_fixed_overrides,
    x_label,
)


# Reference answers for the A31 helix at the bnp_na B-DNA defaults, from
# test/oppo_phos_azbmost/NOTE.md and reproduced by the CLI in this module.
REFERENCE_DSSR_X = 3.0872
REFERENCE_PHENIX_X = 2.7992


def atom_line(serial: int, chain: str, resseq: int, coord: Tuple[float, float, float]) -> str:
    chars = list(" " * 80)
    chars[0:6] = "ATOM  "
    chars[6:11] = f"{serial:5d}"
    chars[12:16] = f"{'P':>4s}"
    chars[17:20] = " DA"
    chars[21] = chain
    chars[22:26] = f"{resseq:4d}"
    chars[30:38] = f"{coord[0]:8.3f}"
    chars[38:46] = f"{coord[1]:8.3f}"
    chars[46:54] = f"{coord[2]:8.3f}"
    chars[54:60] = f"{1.0:6.2f}"
    chars[60:66] = f"{0.0:6.2f}"
    chars[76:78] = " P"
    return "".join(chars) + "\n"


def write_pdb(directory: Path, atoms: List[Tuple[str, int, Tuple[float, float, float]]]) -> Path:
    path = directory / "helix.pdb"
    lines = [atom_line(index + 1, chain, resseq, coord) for index, (chain, resseq, coord) in enumerate(atoms)]
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def opposed_helix(n_bp: int, radius: float = 8.9, rise: float = 3.4) -> List[Tuple[str, int, Tuple[float, float, float]]]:
    """Chain A at angle theta_i, chain B placed exactly across the z-axis."""
    atoms: List[Tuple[str, int, Tuple[float, float, float]]] = []
    for i in range(1, n_bp + 1):
        theta = math.radians(34.2857 * i)
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        atoms.append(("A", i, (x, y, rise * i)))
        atoms.append(("B", n_bp + 2 - i, (-x, -y, rise * i)))
    return atoms


class GeometryTests(unittest.TestCase):
    def test_perfectly_opposed_pairs_score_zero(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = write_pdb(Path(tmp), opposed_helix(6))
            rows = phosphate_pair_rows(path)
            summary = summarize_rows(rows)

        # PDB coordinates are written to 0.001 A, so perfect opposition survives
        # only to about that precision. This is the same floor the real search
        # hits: its best models sit near 1e-4 A, not at zero.
        self.assertEqual(len(rows), 5)  # i = 2..6
        for row in rows:
            self.assertAlmostEqual(float(row["signed_cross_xy"]), 0.0, delta=1e-3)
            self.assertAlmostEqual(float(row["axis_line_distance_A"]), 0.0, delta=1e-4)
            self.assertAlmostEqual(float(row["radial_angle_deg"]), 180.0, delta=1e-3)
            self.assertAlmostEqual(float(row["axis_fraction_on_segment"]), 0.5, delta=1e-4)
            self.assertTrue(row["axis_between_atoms"])
        self.assertAlmostEqual(summary["rms_axis_line_distance_A"], 0.0, delta=1e-4)
        self.assertAlmostEqual(summary["mean_radial_angle_deg"], 180.0, delta=1e-3)

    def test_pairing_rule_is_a_i_with_b_n_plus_2_minus_i(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = write_pdb(Path(tmp), opposed_helix(31))
            rows = phosphate_pair_rows(path)

        pairs = [(int(row["A_res"]), int(row["B_res"])) for row in rows]
        self.assertEqual(pairs[0], (2, 31))  # the A2/B31 case from the reference note
        self.assertEqual(pairs[-1], (31, 2))
        self.assertEqual(len(pairs), 30)
        for a_res, b_res in pairs:
            self.assertEqual(a_res + b_res, 33)

    def test_signed_cross_changes_sign_with_rotation_direction(self) -> None:
        """The search relies on this sign flip to bracket the root."""

        def cross_for(offset_deg: float) -> float:
            atoms: List[Tuple[str, int, Tuple[float, float, float]]] = []
            for i in (1, 2):
                theta = math.radians(90.0 * i)
                atoms.append(("A", i, (8.0 * math.cos(theta), 8.0 * math.sin(theta), 3.4 * i)))
                mirrored = theta + math.pi + math.radians(offset_deg)
                atoms.append(("B", 4 - i, (8.0 * math.cos(mirrored), 8.0 * math.sin(mirrored), 3.4 * i)))
            with tempfile.TemporaryDirectory() as tmp:
                rows = phosphate_pair_rows(write_pdb(Path(tmp), atoms))
            return summarize_rows(rows)["mean_signed_cross_xy"]

        self.assertGreater(cross_for(-3.0), 0.0)
        self.assertLess(cross_for(3.0), 0.0)
        self.assertAlmostEqual(cross_for(0.0), 0.0, places=6)

    def test_axis_line_distance_is_perpendicular_offset(self) -> None:
        # With n_bp = 2 the only pair is A2/B2, so the offset goes on B2.
        # A at (5, 0), B at (-5, 4): the xy line between them misses the axis.
        atoms = [
            ("A", 1, (1.0, 0.0, 0.0)),
            ("B", 1, (-1.0, 0.0, 0.0)),
            ("A", 2, (5.0, 0.0, 3.4)),
            ("B", 2, (-5.0, 4.0, 3.4)),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            rows = phosphate_pair_rows(write_pdb(Path(tmp), atoms))

        row = next(r for r in rows if int(r["A_res"]) == 2)
        # |cross| / |B-A| = |5*4 - 0*(-5)| / sqrt(10^2 + 4^2)
        self.assertAlmostEqual(
            float(row["axis_line_distance_A"]), 20.0 / math.hypot(10.0, 4.0), places=9
        )

    def test_parse_p_atoms_ignores_non_phosphorus(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.pdb"
            line = atom_line(1, "A", 1, (1.0, 2.0, 3.0))
            other = line[:12] + " C1'" + line[16:]
            path.write_text(line + other + "END\n", encoding="utf-8")
            atoms = parse_p_atoms(path)

        self.assertEqual(list(atoms), [("A", 1)])

    def test_missing_partner_is_reported(self) -> None:
        atoms = [("A", 1, (1.0, 0.0, 0.0)), ("A", 2, (0.0, 1.0, 3.4)), ("B", 1, (-1.0, 0.0, 0.0))]
        with tempfile.TemporaryDirectory() as tmp:
            path = write_pdb(Path(tmp), atoms)
            with self.assertRaises(OpposingPhosphateError) as ctx:
                phosphate_pair_rows(path)

        self.assertIn("A2/B2", str(ctx.exception))

    def test_single_chain_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = write_pdb(Path(tmp), [("A", 1, (1.0, 0.0, 0.0)), ("A", 2, (0.0, 1.0, 3.4))])
            with self.assertRaises(OpposingPhosphateError):
                phosphate_pair_rows(path)


class HelperTests(unittest.TestCase):
    def test_quantize_x_rounds_half_up_to_four_decimals(self) -> None:
        self.assertEqual(quantize_x(3.08715), REFERENCE_DSSR_X)
        self.assertEqual(quantize_x(2.79915), REFERENCE_PHENIX_X)
        self.assertEqual(quantize_x(3.0872), REFERENCE_DSSR_X)

    def test_x_label_is_filename_safe(self) -> None:
        self.assertEqual(x_label(REFERENCE_DSSR_X), "3p0872")
        self.assertEqual(x_label(-1.5), "m1p5000")

    def test_validate_fixed_overrides_rejects_x_disp(self) -> None:
        with self.assertRaises(OpposingPhosphateError) as ctx:
            validate_fixed_overrides({"X-disp": 3.0})
        self.assertIn("searched parameter", str(ctx.exception))

    def test_validate_fixed_overrides_rejects_unknown_and_non_numeric(self) -> None:
        with self.assertRaises(OpposingPhosphateError):
            validate_fixed_overrides({"Twist": 30.0})
        with self.assertRaises(OpposingPhosphateError):
            validate_fixed_overrides({"h-Twist": "thirty"})

    def test_validate_fixed_overrides_accepts_known_keys(self) -> None:
        self.assertEqual(validate_fixed_overrides({"h-Twist": "30"}), {"h-Twist": 30.0})
        self.assertEqual(validate_fixed_overrides(None), {})


class FakeModel:
    """Analytic stand-in for a built helix, linear in X-disp around a chosen root."""

    def __init__(
        self,
        *,
        dssr_root: float = REFERENCE_DSSR_X,
        phenix_root: float = REFERENCE_PHENIX_X,
        floor: float = 0.00039,
    ) -> None:
        self.dssr_root = dssr_root
        self.phenix_root = phenix_root
        self.floor = floor
        self.calls: List[Tuple[float, bool]] = []

    def __call__(self, x_disp: float, run_phenix: bool) -> Dict[str, float]:
        self.calls.append((x_disp, run_phenix))
        root = self.phenix_root if run_phenix else self.dssr_root
        slope = -21.5 if run_phenix else -20.15
        cross = slope * (x_disp - root)
        return {
            "mean_signed_cross_xy": cross,
            "mean_abs_cross_xy": abs(cross),
            "max_abs_cross_xy": abs(cross),
            "rms_axis_line_distance_A": self.floor + abs(x_disp - root) * 0.15,
            "mean_axis_line_distance_A": self.floor,
            "max_axis_line_distance_A": self.floor * 2.0,
            "mean_radial_angle_deg": 180.0 - abs(cross) * 0.01,
            "min_radial_angle_deg": 179.98,
            "max_radial_angle_deg": 180.0,
            "min_axis_fraction_on_segment": 0.4,
            "max_axis_fraction_on_segment": 0.6,
        }

    def count(self, run_phenix: bool) -> int:
        return sum(1 for _x, phenix in self.calls if phenix is run_phenix)


class SearchTests(unittest.TestCase):
    def test_raw_dssr_search_finds_reference_root(self) -> None:
        model = FakeModel()
        result = find_opposing_phosphate_xdisp("A31", measure=model)

        self.assertAlmostEqual(float(result["dssr"]["x_disp"]), REFERENCE_DSSR_X, places=4)
        self.assertIsNone(result["phenix"])
        self.assertEqual(result["recommended_mode"], "dssr")
        self.assertAlmostEqual(float(result["recommended_x_disp"]), REFERENCE_DSSR_X, places=4)
        self.assertEqual(model.count(True), 0)

    def test_phenix_refinement_finds_its_own_root(self) -> None:
        """The minimized root sits ~0.29 A away, well outside any window on the raw answer."""
        model = FakeModel()
        result = find_opposing_phosphate_xdisp(
            "A31", measure=model, refine_with_phenix=True, params_file="unused.params"
        )

        self.assertAlmostEqual(float(result["dssr"]["x_disp"]), REFERENCE_DSSR_X, places=4)
        self.assertAlmostEqual(float(result["phenix"]["x_disp"]), REFERENCE_PHENIX_X, places=4)
        self.assertEqual(result["recommended_mode"], "phenix")
        self.assertAlmostEqual(float(result["recommended_x_disp"]), REFERENCE_PHENIX_X, places=4)

    def test_phenix_refinement_is_frugal_with_builds(self) -> None:
        """Each Phenix build costs ~15 s, so the seeded secant must not wander."""
        model = FakeModel()
        find_opposing_phosphate_xdisp(
            "A31", measure=model, refine_with_phenix=True, params_file="unused.params"
        )
        self.assertLessEqual(model.count(True), 12)

    def test_selection_prefers_axis_distance_over_residual_cross(self) -> None:
        """Reproduces the real 2.7991 / 2.7992 tie: 2.7991 has the smaller |cross|,
        but 2.7992 has the smaller axis-line distance and must win."""

        def measure(x_disp: float, run_phenix: bool) -> Dict[str, float]:
            table = {
                2.7991: (0.000377, 0.000416),
                2.7992: (-0.000697, 0.000391),
            }
            cross, distance = table.get(round(x_disp, 4), (-20.15 * (x_disp - 2.79915), 0.5))
            return {
                "mean_signed_cross_xy": cross,
                "mean_abs_cross_xy": abs(cross),
                "max_abs_cross_xy": abs(cross),
                "rms_axis_line_distance_A": distance,
                "mean_axis_line_distance_A": distance,
                "max_axis_line_distance_A": distance,
                "mean_radial_angle_deg": 179.99,
                "min_radial_angle_deg": 179.98,
                "max_radial_angle_deg": 180.0,
                "min_axis_fraction_on_segment": 0.4,
                "max_axis_fraction_on_segment": 0.6,
            }

        result = find_opposing_phosphate_xdisp("A31", measure=measure)
        self.assertAlmostEqual(float(result["dssr"]["x_disp"]), 2.7992, places=4)

    def test_fixed_overrides_reach_the_report(self) -> None:
        model = FakeModel()
        result = find_opposing_phosphate_xdisp(
            "A31", measure=model, fixed_overrides={"h-Twist": 30.0}
        )
        self.assertEqual(result["fixed_overrides"], {"h-Twist": 30.0})
        self.assertIn("h-Twist=30.0000", format_result_report(result))

    def test_report_names_both_answers(self) -> None:
        model = FakeModel()
        result = find_opposing_phosphate_xdisp(
            "A31", measure=model, refine_with_phenix=True, params_file="unused.params"
        )
        report = format_result_report(result)
        self.assertIn(f"{REFERENCE_DSSR_X:.4f}", report)
        self.assertIn(f"{REFERENCE_PHENIX_X:.4f}", report)
        self.assertIn("Recommended X-disp", report)

    def test_scan_that_never_crosses_is_reported(self) -> None:
        model = FakeModel(dssr_root=99.0)
        with self.assertRaises(OpposingPhosphateError) as ctx:
            find_opposing_phosphate_xdisp("A31", measure=model)
        self.assertIn("never crossed", str(ctx.exception))

    def test_blank_sequence_is_rejected(self) -> None:
        with self.assertRaises(OpposingPhosphateError):
            find_opposing_phosphate_xdisp("   ", measure=FakeModel())

    def test_phenix_refinement_requires_a_params_file(self) -> None:
        with self.assertRaises(OpposingPhosphateError) as ctx:
            find_opposing_phosphate_xdisp("A31", refine_with_phenix=True)
        self.assertIn("params file", str(ctx.exception))

    def test_cancellation_stops_the_search(self) -> None:
        model = FakeModel()
        state = {"builds": 0}

        def should_cancel() -> bool:
            state["builds"] += 1
            return state["builds"] > 3

        with self.assertRaises(OpposingPhosphateCancelled):
            find_opposing_phosphate_xdisp("A31", measure=model, should_cancel=should_cancel)
        self.assertLessEqual(model.count(False), 4)

    def test_progress_messages_are_emitted(self) -> None:
        messages: List[str] = []
        find_opposing_phosphate_xdisp("A31", measure=FakeModel(), progress=messages.append)

        self.assertTrue(any("Sequence: A31" in line for line in messages))
        self.assertTrue(any("Raw DSSR X-disp" in line for line in messages))

    def test_repeated_x_disp_values_are_only_built_once(self) -> None:
        model = FakeModel()
        result = find_opposing_phosphate_xdisp("A31", measure=model)

        evaluated = [x for x, phenix in model.calls if not phenix]
        self.assertEqual(len(evaluated), len(set(evaluated)))
        self.assertEqual(model.count(False), int(result["build_counts"]["dssr"]))

    def test_max_phenix_builds_is_enforced(self) -> None:
        model = FakeModel(phenix_root=6.5)  # far from the raw seed, so the secant needs room
        with self.assertRaises(OpposingPhosphateError) as ctx:
            find_opposing_phosphate_xdisp(
                "A31",
                measure=model,
                refine_with_phenix=True,
                params_file="unused.params",
                max_phenix_builds=2,
            )
        self.assertIn("without converging", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
