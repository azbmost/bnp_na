from __future__ import annotations

import argparse
import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np

from bnp_na_lib import core_BP
from bnp_na_lib.build_bp import build_bp_structure
from bnp_na_lib.build_common import PipelineError
from bnp_na_lib.build_pdna import build_pdna_model, model_pdb_lines
from bnp_na_lib.make_BP import (
    detect_helix_type,
    load_helix,
    new_core,
    parse_types,
    run_pipeline,
    screw_axis_line,
)
from bnp_na_lib.standard_bases import (
    DNA_COMPLEMENT,
    DNA_RESIDUE_NAMES,
    base_code,
    graft_base,
    is_base_atom,
    kabsch,
)

ORDER = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'")


def core_residue_atoms(chain: str, res_seq: int):
    atoms = {}
    res_name = ""
    for line in core_BP.PDB_LINES:
        if line.startswith("ATOM") and line[21] == chain and int(line[22:26]) == res_seq:
            atoms[line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            res_name = line[17:20].strip()
    return res_name, atoms


def b_dna_from_core(seq: str):
    """Ideal B-DNA rebuilt from the core's two B pairs and the screw between them."""
    (e_out, f_out), (e_in, f_in) = core_BP.B_END_BASEPAIRS[1], core_BP.B_END_BASEPAIRS[0]
    outer = [core_residue_atoms(*e_out), core_residue_atoms(*f_out)]
    inner = [core_residue_atoms(*e_in), core_residue_atoms(*f_in)]
    src, dst = [], []
    for (_n0, a0), (_n1, a1) in zip(outer, inner):
        for name in ORDER:
            src.append(a0[name])
            dst.append(a1[name])
    rot, trans, rmsd = kabsch(np.array(src), np.array(dst))
    assert rmsd < 0.01, rmsd
    strand_a, strand_b = [], []
    pair = outer
    for base in seq:
        (name_e, atoms_e), (name_f, atoms_f) = pair
        grafted_e = dict(graft_base(atoms_e, base_code(name_e), base)[0])
        grafted_f = dict(graft_base(atoms_f, base_code(name_f), DNA_COMPLEMENT[base])[0])
        strand_a.append((DNA_RESIDUE_NAMES[base], list(grafted_e.items())))
        strand_b.insert(0, (DNA_RESIDUE_NAMES[DNA_COMPLEMENT[base]], list(grafted_f.items())))
        pair = [(n, {k: rot @ v + trans for k, v in a.items()}) for n, a in pair]
    return {"A": strand_a, "B": strand_b}


def write_model(model, path: Path) -> Path:
    path.write_text("".join(model_pdb_lines(model)), encoding="utf-8")
    return path


def read_chains(path: Path):
    chains = {}
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            chain, res_seq = line[21], int(line[22:26])
            res = chains.setdefault(chain, {}).setdefault(res_seq, {"name": line[17:20].strip(), "atoms": {}})
            res["atoms"][line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return {c: [chains[c][k] for k in sorted(chains[c])] for c in chains}


def run_quiet(files, out: Path, types: str = ""):
    args = argparse.Namespace(dna_files=[str(f) for f in files], out=str(out), types=types,
                              axis_mode="codirectional", axis_source="screw")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        result = run_pipeline(args)
    return result, buf.getvalue()


class BPBuilderTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.b_seq = "GATTCAGC"
        self.p_seq1, self.p_seq2 = "ACGTTGCA", "TTTGGGCC"
        self.b_pdb = write_model(b_dna_from_core(self.b_seq), self.dir / "B8.pdb")
        self.p_pdb = write_model(build_pdna_model(self.p_seq1, self.p_seq2), self.dir / "P8.pdb")

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def assert_ligated(self, path: Path, expected_sequences, max_link: float = 1.75) -> None:
        chains = read_chains(path)
        self.assertEqual(sorted(chains), ["A", "B"])
        seqs = ["".join(base_code(r["name"]) for r in chains[c]) for c in ("A", "B")]
        self.assertEqual(seqs, list(expected_sequences))
        for residues in chains.values():
            for lower, upper in zip(residues, residues[1:]):
                d = float(np.linalg.norm(lower["atoms"]["O3'"] - upper["atoms"]["P"]))
                self.assertLess(d, max_link, msg=f"O3'-P {d:.3f} A")
        # No heavy-atom clashes between residues; covalent neighbors count only
        # for contacts that involve a base atom.
        flat = []
        for c, residues in chains.items():
            for i, r in enumerate(residues):
                names = list(r["atoms"])
                flat.append((c, i, np.array([r["atoms"][n] for n in names]),
                             np.array([is_base_atom(n) for n in names])))
        worst = 9.0
        for a in range(len(flat)):
            ca, ia, xa, base_a = flat[a]
            for b in range(a + 1, len(flat)):
                cb, ib, xb, base_b = flat[b]
                d = np.linalg.norm(xa[:, None, :] - xb[None, :, :], axis=2)
                if ca == cb and abs(ia - ib) == 1:
                    d = np.where(base_a[:, None] | base_b[None, :], d, np.inf)
                worst = min(worst, float(d.min()))
        self.assertGreater(worst, 2.6)

    def test_helix_types_are_detected(self) -> None:
        self.assertEqual(detect_helix_type(load_helix(self.b_pdb, "s1")), "B")
        self.assertEqual(detect_helix_type(load_helix(self.p_pdb, "s2")), "P")

    def test_p_dna_screw_axis(self) -> None:
        line = screw_axis_line(load_helix(self.p_pdb, "s2"))
        self.assertGreater(abs(line.direction[2]), 0.9999)
        self.assertLess(float(np.hypot(line.point[0], line.point[1])), 1e-3)
        self.assertIn("rise=6.100", line.detail)

    def test_b_then_p_keeps_every_base(self) -> None:
        out = self.dir / "BP.pdb"
        result, log = run_quiet([self.b_pdb, self.p_pdb], out)
        self.assertEqual(result["types"], ["B", "P"])
        rev_b = "".join(DNA_COMPLEMENT[b] for b in reversed(self.b_seq))
        self.assert_ligated(out, [self.b_seq + self.p_seq1, self.p_seq2 + rev_b])
        self.assertTrue(Path(str(result["pdb_raw"])).exists())
        self.assertIn("Ligation graph: components=2", log)

    def test_p_then_b(self) -> None:
        out = self.dir / "PB.pdb"
        result, _log = run_quiet([self.p_pdb, self.b_pdb], out)
        self.assertEqual(result["types"], ["P", "B"])
        rev_b = "".join(DNA_COMPLEMENT[b] for b in reversed(self.b_seq))
        # Chain A of the first input (P strand I) continues into chain A of the B-DNA.
        self.assert_ligated(out, [self.p_seq1 + self.b_seq, rev_b + self.p_seq2])

    def test_b_p_b(self) -> None:
        b2_seq = "CCGATAGC"
        b2 = write_model(b_dna_from_core(b2_seq), self.dir / "B8b.pdb")
        out = self.dir / "BPB.pdb"
        run_quiet([self.b_pdb, self.p_pdb, b2], out)
        rev = lambda s: "".join(DNA_COMPLEMENT[b] for b in reversed(s))  # noqa: E731
        self.assert_ligated(out, [self.b_seq + self.p_seq1 + b2_seq, rev(b2_seq) + self.p_seq2 + rev(self.b_seq)])

    def test_types_must_alternate(self) -> None:
        with self.assertRaises(SystemExit):
            run_quiet([self.b_pdb, self.b_pdb], self.dir / "BB.pdb")
        with self.assertRaises(SystemExit):
            run_quiet([self.b_pdb, self.p_pdb], self.dir / "x.pdb", types="BB")

    def test_short_p_helix_is_rejected(self) -> None:
        short = write_model(build_pdna_model("ACG", "CGT"), self.dir / "P3.pdb")
        with self.assertRaises(SystemExit):
            run_quiet([self.b_pdb, short], self.dir / "short.pdb")

    def test_p_between_two_junctions_needs_seven_levels(self) -> None:
        b2 = write_model(b_dna_from_core("CCGATAGC"), self.dir / "B8b.pdb")
        six = write_model(build_pdna_model("ACGTTG", "CAACGT"), self.dir / "P6.pdb")
        with self.assertRaises(SystemExit):
            run_quiet([self.b_pdb, six, b2], self.dir / "bpb6.pdb")
        seven = write_model(build_pdna_model("ACGTTGA", "TCAACGT"), self.dir / "P7.pdb")
        out = self.dir / "bpb7.pdb"
        run_quiet([self.b_pdb, seven, b2], out)
        rev = lambda s: "".join(DNA_COMPLEMENT[b] for b in reversed(s))  # noqa: E731
        self.assert_ligated(out, [self.b_seq + "ACGTTGA" + "CCGATAGC", rev("CCGATAGC") + "TCAACGT" + rev(self.b_seq)])

    def test_p_b_p_with_a_five_bp_b_dna(self) -> None:
        b5 = write_model(b_dna_from_core("GATCC"), self.dir / "B5.pdb")
        p2 = write_model(build_pdna_model("TTGCA", "GGCAT"), self.dir / "P5.pdb")
        out = self.dir / "pbp.pdb"
        run_quiet([self.p_pdb, b5, p2], out)
        self.assert_ligated(out, [self.p_seq1 + "GATCC" + "TTGCA", "GGCAT" + "GGATC" + self.p_seq2])

    def test_wrapper_reports_errors_as_pipeline_errors(self) -> None:
        with self.assertRaises(PipelineError):
            build_bp_structure([self.b_pdb], self.dir / "one.pdb")
        with self.assertRaises(PipelineError):
            build_bp_structure([self.b_pdb, self.b_pdb], self.dir / "bb.pdb", axis_source="screw")

    def test_parse_types(self) -> None:
        self.assertEqual(parse_types("b p,b", 3), ["B", "P", "B"])
        with self.assertRaises(ValueError):
            parse_types("BZ", 2)
        with self.assertRaises(ValueError):
            parse_types("BPB", 2)


class CoreTest(unittest.TestCase):
    def test_core_layout(self) -> None:
        core = new_core("c")
        self.assertEqual(sorted({key[1] for key in core.residues}), ["E", "F"])
        self.assertEqual(len(core.residues), 12)
        self.assertEqual(len(core.basepairs), 6)

    def test_core_linkages_are_bonded(self) -> None:
        chains = {}
        for line in core_BP.PDB_LINES:
            if line.startswith("ATOM"):
                res = chains.setdefault(line[21], {}).setdefault(int(line[22:26]), {})
                res[line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        for chain in ("E", "F"):
            residues = [chains[chain][k] for k in sorted(chains[chain])]
            for lower, upper in zip(residues, residues[1:]):
                self.assertLess(float(np.linalg.norm(lower["O3'"] - upper["P"])), 1.7)


if __name__ == "__main__":
    unittest.main()
