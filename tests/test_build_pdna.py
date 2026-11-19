from __future__ import annotations

import itertools
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from bnp_na_lib import core_PDNA
from bnp_na_lib.build_common import PipelineError
from bnp_na_lib.build_pdna import (
    PURINE_CHI_OFFSET_DEG,
    build_pdna,
    build_pdna_model,
    closest_nonbonded_contact,
    level_pairs,
    model_pdb_lines,
    resolve_pdna_sequences,
    superpose_onto_reference,
)
from bnp_na_lib.standard_bases import (
    PURINES,
    base_frame,
    chi_deg,
    graft_base,
    kabsch,
    standard_base_atoms,
)

BACKBONE = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")


def radius(xyz: np.ndarray) -> float:
    return float(math.hypot(xyz[0], xyz[1]))


def atoms_of(residue):
    return {name: xyz for name, xyz in residue[1]}


class ResolveSequencesTest(unittest.TestCase):
    def test_strand_two_defaults_to_reverse_complement(self) -> None:
        seq1, seq2, defaulted = resolve_pdna_sequences("ACGTTG")
        self.assertEqual(seq1, "ACGTTG")
        self.assertEqual(seq2, "CAACGT")
        self.assertTrue(defaulted)

    def test_blank_strand_two_also_defaults(self) -> None:
        self.assertEqual(resolve_pdna_sequences("AAC", "  ")[1], "GTT")

    def test_explicit_strand_two_is_kept(self) -> None:
        seq1, seq2, defaulted = resolve_pdna_sequences("A3C2", "GGGGG")
        self.assertEqual((seq1, seq2, defaulted), ("AAACC", "GGGGG", False))

    def test_unequal_lengths_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "same length"):
            resolve_pdna_sequences("ACGT", "ACG")

    def test_rna_letters_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            resolve_pdna_sequences("ACGU")
        with self.assertRaises(ValueError):
            resolve_pdna_sequences("")

    def test_level_pairs_follow_the_register(self) -> None:
        self.assertEqual(level_pairs("ACG", "TTC")[0], (1, "A", 3, "C"))
        self.assertEqual(level_pairs("ACG", "TTC")[-1], (3, "G", 1, "T"))


class TemplateTest(unittest.TestCase):
    def test_template_level_is_the_dyad_image_of_itself(self) -> None:
        level = core_PDNA.level_atoms()
        strand_i = dict(level["I"])
        strand_ii = dict(level["II"])
        for name, (x, y, z) in strand_i.items():
            self.assertEqual(strand_ii[name], (x, -y, -z))

    def test_helical_parameters(self) -> None:
        params = core_PDNA.helical_parameters()
        self.assertAlmostEqual(params["h-Rise (Å/nt)"], 6.10)
        self.assertAlmostEqual(params["h-Twist (°/nt)"], 131.84)
        self.assertAlmostEqual(params["nt per turn"], 2.7306, places=4)
        self.assertAlmostEqual(params["Pitch (Å)"], 16.6566, places=3)


class ModelGeometryTest(unittest.TestCase):
    def setUp(self) -> None:
        self.seq1, self.seq2, _ = resolve_pdna_sequences("ACGTGCATTG", "GGATCCAATC")
        self.model = build_pdna_model(self.seq1, self.seq2)

    def test_chain_lengths_and_residue_names(self) -> None:
        self.assertEqual([r[0] for r in self.model["A"]], ["D" + b for b in self.seq1])
        self.assertEqual([r[0] for r in self.model["B"]], ["D" + b for b in self.seq2])

    def test_backbone_follows_the_template_screw(self) -> None:
        twist = math.radians(core_PDNA.TWIST_DEG)
        rot = np.array([[math.cos(twist), -math.sin(twist), 0.0], [math.sin(twist), math.cos(twist), 0.0], [0, 0, 1]])
        shift = np.array([0.0, 0.0, core_PDNA.RISE])
        for chain, step in (("A", 1), ("B", -1)):
            residues = self.model[chain]
            pairs = zip(residues, residues[1:]) if step == 1 else zip(residues[1:], residues)
            for lower, upper in pairs:
                low, up = atoms_of(lower), atoms_of(upper)
                for name in BACKBONE:
                    np.testing.assert_allclose(rot @ low[name] + shift, up[name], atol=1e-9)

    def test_axis_radii_match_the_template(self) -> None:
        for chain in ("A", "B"):
            for residue in self.model[chain]:
                atoms = atoms_of(residue)
                self.assertAlmostEqual(radius(atoms["P"]), core_PDNA.P_RADIUS, delta=0.01)
                self.assertAlmostEqual(radius(atoms["C1'"]), core_PDNA.C1_RADIUS, delta=0.01)

    def test_strand_one_runs_up_and_level_one_is_at_the_origin(self) -> None:
        first = atoms_of(self.model["A"][0])
        last = atoms_of(self.model["A"][-1])
        self.assertGreater(last["C1'"][2], first["C1'"][2])
        level_one = (first["C1'"][2] + atoms_of(self.model["B"][-1])["C1'"][2]) / 2.0
        self.assertAlmostEqual(level_one, 0.0, places=6)

    def test_levels_hold_one_nucleotide_of_each_strand(self) -> None:
        n = len(self.seq1)
        for k in range(1, n + 1):
            a = atoms_of(self.model["A"][k - 1])["C1'"]
            b = atoms_of(self.model["B"][n - k])["C1'"]
            self.assertLess(abs(a[2] - b[2]), core_PDNA.LEVEL_C1_DZ + 1e-6)
            self.assertAlmostEqual(a[2], (k - 1) * core_PDNA.RISE, delta=0.2)

    def test_linkages_and_glycosidic_bonds(self) -> None:
        for chain in ("A", "B"):
            residues = self.model[chain]
            for lower, upper in zip(residues, residues[1:]):
                d = float(np.linalg.norm(atoms_of(lower)["O3'"] - atoms_of(upper)["P"]))
                self.assertAlmostEqual(d, 1.641, delta=0.002)
            for res_name, atoms in residues:
                names = dict(atoms)
                n_atom = "N9" if res_name[-1] in PURINES else "N1"
                self.assertAlmostEqual(float(np.linalg.norm(names["C1'"] - names[n_atom])), 1.47, delta=0.012)

    def test_bases_keep_standard_geometry(self) -> None:
        for chain in ("A", "B"):
            for res_name, atoms in self.model[chain]:
                names = dict(atoms)
                standard = standard_base_atoms(res_name[-1])
                moving = np.array([xyz for _n, xyz in standard])
                target = np.array([names[n] for n, _xyz in standard])
                self.assertLess(kabsch(moving, target)[2], 1e-6)

    def test_chi_follows_the_template(self) -> None:
        template = {name: np.array(xyz) for name, xyz in core_PDNA.level_atoms()["I"]}
        template_chi = chi_deg(dict(graft_base(template, "C", "C")[0]), "C")
        for chain in ("A", "B"):
            for res_name, atoms in self.model[chain]:
                base = res_name[-1]
                expected = template_chi + (PURINE_CHI_OFFSET_DEG if base in PURINES else 0.0)
                self.assertAlmostEqual(chi_deg(dict(atoms), base), expected, delta=0.5)

    def test_no_clashes_for_any_base_combination(self) -> None:
        for b1, b2 in itertools.product("ACGT", repeat=2):
            model = build_pdna_model(b1 * 4, b2 * 4)
            residues = [(chain, i, atoms) for chain in ("A", "B") for i, (_n, atoms) in enumerate(model[chain])]
            contact, _where = closest_nonbonded_contact(residues)
            self.assertGreater(contact, 3.4, msg=f"{b1}/{b2}")

    def test_pdb_lines_are_well_formed(self) -> None:
        lines = model_pdb_lines(self.model)
        atom_lines = [line for line in lines if line.startswith("ATOM")]
        self.assertEqual(atom_lines[0][12:16], " P  ")
        self.assertEqual(atom_lines[0][17:20], " DA")
        self.assertEqual(atom_lines[0][21], "A")
        self.assertTrue(all(len(line.rstrip("\n")) == 78 for line in atom_lines))
        self.assertEqual(sum(line.startswith("TER") for line in lines), 2)
        self.assertEqual(lines[-1], "END\n")


class GraftBaseTest(unittest.TestCase):
    def setUp(self) -> None:
        self.template = {name: np.array(xyz) for name, xyz in core_PDNA.level_atoms()["I"]}

    def test_round_trip_through_another_base(self) -> None:
        direct = dict(graft_base(self.template, "C", "C")[0])
        via_g = dict(graft_base(dict(graft_base(self.template, "C", "G")[0]), "G", "C")[0])
        for name, xyz in direct.items():
            np.testing.assert_allclose(via_g[name], xyz, atol=1e-9)

    def test_backbone_is_untouched(self) -> None:
        grafted = dict(graft_base(self.template, "C", "A")[0])
        for name in BACKBONE:
            np.testing.assert_allclose(grafted[name], self.template[name])
        self.assertNotIn("N4", grafted)
        self.assertIn("N6", grafted)

    def test_pivot_fit_keeps_the_glycosidic_bond(self) -> None:
        rot, trans, _rmsd = base_frame(self.template, "C")
        self.assertAlmostEqual(float(np.linalg.det(rot)), 1.0, places=9)
        grafted = dict(graft_base(self.template, "C", "T")[0])
        standard = dict(standard_base_atoms("T", include_c1=True))
        expected = float(np.linalg.norm(standard["N1"] - standard["C1'"]))
        self.assertAlmostEqual(float(np.linalg.norm(grafted["N1"] - grafted["C1'"])), expected, places=6)


class BuildPipelineTest(unittest.TestCase):
    def test_build_without_phenix_writes_an_aligned_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = build_pdna("ACGTACGT", "P8test", tmp, run_phenix=False, run_regularize_phosphates=False)
            aligned = Path(result["pdb_aligned"])
            self.assertTrue(aligned.exists())
            self.assertEqual(aligned.name, "P8test-1PFI_out_aligned2Z.pdb")
            self.assertEqual(result["strand2"], "ACGTACGT")
            self.assertTrue(result["strand2_default"])
            self.assertEqual(result["length"], 8)
            self.assertIn("Fixed helical parameters", result["log_text"])
            atom_lines = [l for l in aligned.read_text().splitlines() if l.startswith("ATOM")]
            self.assertEqual({l[21] for l in atom_lines}, {"A", "B"})

    def test_minimization_needs_a_params_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(PipelineError):
                build_pdna("ACGT", None, tmp, run_phenix=True, params_file=None)

    def test_bad_sequence_raises_pipeline_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(PipelineError):
                build_pdna("ACGT", None, tmp, seq_strand2="AC")

    def test_superpose_restores_the_frame(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = build_pdna("GATTACA", "P7", tmp, run_phenix=False, run_regularize_phosphates=False)
            reference = Path(result["pdb_normalized"])
            angle = math.radians(37.0)
            rot = np.array([[1, 0, 0], [0, math.cos(angle), -math.sin(angle)], [0, math.sin(angle), math.cos(angle)]])
            moved_lines = []
            for line in reference.read_text().splitlines(keepends=True):
                if line.startswith("ATOM"):
                    xyz = rot @ np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) + [5.0, -3.0, 9.0]
                    line = line[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + line[54:]
                moved_lines.append(line)
            moved = Path(tmp) / "moved.pdb"
            moved.write_text("".join(moved_lines))
            fit = superpose_onto_reference(moved, reference, Path(tmp) / "back.pdb")
            self.assertLess(fit["rmsd"], 0.002)


if __name__ == "__main__":
    unittest.main()
