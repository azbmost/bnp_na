from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from bnp_na_lib.combine_pdb import (
    CombinePDBError,
    combine_pdb_files,
    list_pdb_chains,
    parse_chain_selection,
)


def atom_line(serial: int, chain: str, resseq: int, x: float = 0.0) -> str:
    chars = list(" " * 80)
    chars[0:6] = "ATOM  "
    chars[6:11] = f"{serial:5d}"
    chars[12:16] = " P  "
    chars[17:20] = " DA"
    chars[21] = chain
    chars[22:26] = f"{resseq:4d}"
    chars[30:38] = f"{x:8.3f}"
    chars[38:46] = f"{0.0:8.3f}"
    chars[46:54] = f"{0.0:8.3f}"
    chars[54:60] = f"{1.0:6.2f}"
    chars[60:66] = f"{0.0:6.2f}"
    chars[76:78] = " P"
    return "".join(chars) + "\n"


def anisou_line(serial: int, chain: str, resseq: int) -> str:
    chars = list(" " * 80)
    chars[0:6] = "ANISOU"
    chars[6:11] = f"{serial:5d}"
    chars[12:16] = " P  "
    chars[17:20] = " DA"
    chars[21] = chain
    chars[22:26] = f"{resseq:4d}"
    return "".join(chars) + "\n"


def ter_line(serial: int, chain: str, resseq: int) -> str:
    chars = list(" " * 80)
    chars[0:6] = "TER   "
    chars[6:11] = f"{serial:5d}"
    chars[17:20] = " DA"
    chars[21] = chain
    chars[22:26] = f"{resseq:4d}"
    return "".join(chars) + "\n"


def link_line(chain_1: str, resseq_1: int, chain_2: str, resseq_2: int) -> str:
    chars = list(" " * 80)
    chars[0:6] = "LINK  "
    chars[12:16] = " P  "
    chars[17:20] = " DA"
    chars[21] = chain_1
    chars[22:26] = f"{resseq_1:4d}"
    chars[42:46] = " O3'"
    chars[47:50] = " DT"
    chars[51] = chain_2
    chars[52:56] = f"{resseq_2:4d}"
    chars[73:78] = f"{1.61:5.2f}"
    return "".join(chars) + "\n"


def het_line(chain: str, resseq: int) -> str:
    chars = list(" " * 30)
    chars[0:6] = "HET   "
    chars[7:10] = "X33"
    chars[12] = chain
    chars[13:17] = f"{resseq:4d}"
    chars[20:25] = f"{3:5d}"
    return "".join(chars) + "\n"


class CombinePDBTests(unittest.TestCase):
    def test_assigns_chains_and_remaps_serial_records(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            output = root / "combined.pdb"
            first.write_text(
                "".join(
                    [
                        atom_line(10, "X", 1, 1.25),
                        anisou_line(10, "X", 1),
                        ter_line(11, "X", 1),
                        atom_line(20, "Y", 1, 2.5),
                        ter_line(21, "Y", 1),
                        "CONECT   10   20\n",
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            second.write_text(
                "".join(
                    [
                        "REMARK 950 RE_SCRIPT SOFTWARE name=re_helix version=V3.17\n",
                        "REMARK 950 RE_SCRIPT COMMAND text=re_helix.py input.pdb 13A A:1:DA\n",
                        "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=A start=A:1:DA end=A:2:DT count=2\n",
                        "REMARK 950 RE_SCRIPT CHAIN_RESIDUES chain=A part=1/1 residues=A:1:DA,A:2:DT\n",
                        "REMARK 950 RE_SCRIPT JUNCTION residues=A:1:DA,A:2:DT core=A:1:DA "
                        "original_prev1=A:9:DA\n",
                        "REMARK 950 RE_SCRIPT SPECIAL event=standalone_x33 chain=A source=A:9 "
                        "residue=A:2:X33\n",
                        "REMARK BNP_NA_L_RESIDUE KIND L-DNA CHAIN A RESSEQ 1 ICODE . RESNAME DA\n",
                        "REMARK    A.DA1 [A-T] A.DT2\n",
                        het_line("A", 2),
                        "HETNAM     X33 3'-3' PHOSPHODIESTER LINKER PHOSPHATE\n",
                        link_line("A", 1, "A", 2),
                        atom_line(5, "A", 1, 3.75),
                        atom_line(6, " ", 2, 4.0),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )

            result = combine_pdb_files([first, second], output)
            lines = output.read_text(encoding="utf-8").splitlines()
            atoms = [line for line in lines if line.startswith("ATOM  ")]
            anisou = next(line for line in lines if line.startswith("ANISOU"))
            conect = next(line for line in lines if line.startswith("CONECT"))
            link = next(line for line in lines if line.startswith("LINK  "))
            het = next(line for line in lines if line.startswith("HET   "))

            self.assertEqual([line[21] for line in atoms], ["A", "B", "C", "D"])
            self.assertEqual([int(line[6:11]) for line in atoms], [1, 3, 5, 6])
            self.assertEqual(int(anisou[6:11]), 1)
            self.assertEqual(anisou[21], "A")
            self.assertEqual([int(conect[index : index + 5]) for index in range(6, len(conect), 5)], [1, 3])
            self.assertEqual([float(line[30:38]) for line in atoms], [1.25, 2.5, 3.75, 4.0])
            self.assertEqual((link[21], link[51]), ("C", "C"))
            self.assertEqual(het[12], "C")
            self.assertIn(
                "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=C start=C:1:DA end=C:2:DT count=2", lines
            )
            self.assertIn(
                "REMARK 950 RE_SCRIPT CHAIN_RESIDUES chain=C part=1/1 residues=C:1:DA,C:2:DT", lines
            )
            self.assertIn(
                "REMARK 950 RE_SCRIPT JUNCTION residues=C:1:DA,C:2:DT core=C:1:DA "
                "original_prev1=A:9:DA",
                lines,
            )
            self.assertIn(
                "REMARK 950 RE_SCRIPT SPECIAL event=standalone_x33 chain=C source=A:9 "
                "residue=C:2:X33",
                lines,
            )
            self.assertIn("REMARK 950 RE_SCRIPT COMMAND text=re_helix.py input.pdb 13A A:1:DA", lines)
            self.assertIn(
                "REMARK BNP_NA_L_RESIDUE KIND L-DNA CHAIN C RESSEQ 1 ICODE . RESNAME DA", lines
            )
            self.assertIn("REMARK    C.DA1 [A-T] C.DT2", lines)
            self.assertEqual(lines[-1], "END")
            self.assertEqual(sum(line == "END" for line in lines), 1)
            self.assertEqual(result.chain_count, 4)
            self.assertEqual(result.atom_count, 4)
            self.assertEqual(result.remark_count, 8)
            self.assertEqual(result.link_count, 1)

    def test_allows_the_same_input_file_to_be_repeated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.pdb"
            output = root / "combined.pdb"
            source.write_text(
                "".join(
                    [
                        "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=X start=X:1:DA end=X:1:DA count=1\n",
                        "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=Y start=Y:1:DT end=Y:1:DT count=1\n",
                        link_line("X", 1, "Y", 1),
                        atom_line(10, "X", 1, 1.0),
                        atom_line(20, "Y", 1, 2.0),
                        "CONECT   10   20\n",
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )

            result = combine_pdb_files([source, source], output)
            lines = output.read_text(encoding="utf-8").splitlines()
            atoms = [line for line in lines if line.startswith("ATOM  ")]
            links = [line for line in lines if line.startswith("LINK  ")]
            conect = [line for line in lines if line.startswith("CONECT")]

            self.assertEqual([line[21] for line in atoms], ["A", "B", "C", "D"])
            self.assertEqual([(line[21], line[51]) for line in links], [("A", "B"), ("C", "D")])
            self.assertEqual(
                [[int(line[index : index + 5]) for index in range(6, len(line), 5)] for line in conect],
                [[1, 2], [3, 4]],
            )
            self.assertIn(
                "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=A start=A:1:DA end=A:1:DA count=1", lines
            )
            self.assertIn(
                "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=C start=C:1:DA end=C:1:DA count=1", lines
            )
            self.assertEqual(result.input_pdbs, [source.resolve(), source.resolve()])
            self.assertEqual(result.chain_count, 4)
            self.assertEqual(result.atom_count, 4)
            self.assertEqual(result.remark_count, 4)
            self.assertEqual(result.link_count, 2)

    def test_chain_selection_keeps_only_the_requested_chains(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            output = root / "combined.pdb"
            first.write_text(
                "".join(
                    [
                        atom_line(10, "X", 1, 1.25),
                        anisou_line(10, "X", 1),
                        ter_line(11, "X", 1),
                        atom_line(20, "Y", 1, 2.5),
                        ter_line(21, "Y", 1),
                        "CONECT   10   20\n",
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            second.write_text(
                "".join(
                    [
                        "REMARK 950 RE_SCRIPT SOFTWARE name=re_helix version=V3.17\n",
                        "REMARK 950 RE_SCRIPT COMMAND text=re_helix.py input.pdb 13A A:1:DA\n",
                        "REMARK 950 RE_SCRIPT CHAIN_RANGE chain=A start=A:1:DA end=A:2:DT count=2\n",
                        "REMARK 950 RE_SCRIPT SPECIAL event=standalone_x33 chain=A source=A:9 "
                        "residue=A:2:X33\n",
                        "REMARK BNP_NA_L_RESIDUE KIND L-DNA CHAIN A RESSEQ 1 ICODE . RESNAME DA\n",
                        "REMARK    A.DA1 [A-T] A.DT2\n",
                        het_line("A", 2),
                        "HETNAM     X33 3'-3' PHOSPHODIESTER LINKER PHOSPHATE\n",
                        link_line("A", 1, "A", 2),
                        atom_line(5, "A", 1, 3.75),
                        atom_line(6, " ", 2, 4.0),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )

            result = combine_pdb_files([first, second], output, ["Y", "_"])
            lines = output.read_text(encoding="utf-8").splitlines()
            atoms = [line for line in lines if line.startswith("ATOM  ")]

            self.assertEqual([line[21] for line in atoms], ["A", "B"])
            self.assertEqual([int(line[6:11]) for line in atoms], [1, 3])
            self.assertEqual([float(line[30:38]) for line in atoms], [2.5, 4.0])
            self.assertEqual([line[21] for line in lines if line.startswith("TER   ")], ["A"])

            # The excluded chains take their coordinates and metadata with them.
            self.assertEqual([line for line in lines if line.startswith("ANISOU")], [])
            self.assertEqual([line for line in lines if line.startswith("CONECT")], [])
            self.assertEqual([line for line in lines if line.startswith("LINK  ")], [])
            self.assertEqual([line for line in lines if line.startswith("HET   ")], [])
            self.assertEqual([line for line in lines if line.startswith("HETNAM")], [])
            self.assertNotIn("CHAIN_RANGE", output.read_text(encoding="utf-8"))
            self.assertNotIn("BNP_NA_L_RESIDUE", output.read_text(encoding="utf-8"))
            self.assertNotIn("[A-T]", output.read_text(encoding="utf-8"))

            # Provenance that names no current chain is still preserved verbatim.
            self.assertIn("REMARK 950 RE_SCRIPT SOFTWARE name=re_helix version=V3.17", lines)
            self.assertIn("REMARK 950 RE_SCRIPT COMMAND text=re_helix.py input.pdb 13A A:1:DA", lines)

            self.assertIn(
                "REMARK BNP_NA_COMBINE_PDB INPUT 1 first.pdb CHAINS Y->A SKIPPED X", lines
            )
            self.assertIn(
                "REMARK BNP_NA_COMBINE_PDB INPUT 2 second.pdb CHAINS (blank)->B SKIPPED A", lines
            )
            self.assertEqual(result.chain_count, 2)
            self.assertEqual(result.atom_count, 2)
            self.assertEqual(result.remark_count, 2)
            self.assertEqual(result.link_count, 0)
            self.assertEqual(result.mappings[0].excluded_chains, ("X",))
            self.assertEqual(result.mappings[1].excluded_chains, ("A",))

    def test_chain_selection_prunes_dropped_conect_partners(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            output = root / "combined.pdb"
            first.write_text(
                "".join(
                    [
                        atom_line(1, "P", 1, 1.0),
                        atom_line(2, "P", 2, 2.0),
                        atom_line(3, "Q", 1, 3.0),
                        "CONECT    1    2    3\n",
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            second.write_text(atom_line(1, "Z", 1, 4.0), encoding="utf-8")

            result = combine_pdb_files([first, second], output, ["P", None])
            lines = output.read_text(encoding="utf-8").splitlines()
            conect = [line for line in lines if line.startswith("CONECT")]

            self.assertEqual(conect, ["CONECT    1    2"])
            self.assertEqual(
                [line[21] for line in lines if line.startswith("ATOM  ")], ["A", "A", "B"]
            )
            self.assertEqual(result.chain_count, 2)
            self.assertEqual(result.atom_count, 3)

    def test_chain_selection_keeps_ter_records_with_a_blank_chain_column(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            output = root / "combined.pdb"
            first.write_text(
                "".join(
                    [
                        atom_line(1, "X", 1, 1.0),
                        ter_line(2, " ", 1),
                        atom_line(3, "Y", 1, 2.0),
                        ter_line(4, " ", 1),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            second.write_text(atom_line(1, "Z", 1, 3.0), encoding="utf-8")

            combine_pdb_files([first, second], output, ["Y", "all"])
            lines = output.read_text(encoding="utf-8").splitlines()
            serial_records = [
                line for line in lines if line.startswith(("ATOM  ", "TER   "))
            ]

            self.assertEqual(
                [(line[:6], int(line[6:11]), line[21]) for line in serial_records],
                [("ATOM  ", 1, "A"), ("TER   ", 2, "A"), ("ATOM  ", 3, "B")],
            )

    def test_rejects_a_chain_that_is_not_in_the_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            first.write_text(atom_line(1, "X", 1), encoding="utf-8")
            second.write_text(atom_line(1, "Y", 1), encoding="utf-8")
            with self.assertRaisesRegex(CombinePDBError, r"has no chain Z; its chains are X\."):
                combine_pdb_files([first, second], root / "combined.pdb", ["Z", None])

    def test_rejects_a_chain_selection_count_that_misses_an_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            first.write_text(atom_line(1, "X", 1), encoding="utf-8")
            second.write_text(atom_line(1, "Y", 1), encoding="utf-8")
            with self.assertRaisesRegex(CombinePDBError, "given for 1 inputs, but 2 input PDB files"):
                combine_pdb_files([first, second], root / "combined.pdb", ["X"])

    def test_parse_chain_selection_accepts_the_documented_forms(self) -> None:
        self.assertIsNone(parse_chain_selection(""))
        self.assertIsNone(parse_chain_selection("   "))
        self.assertIsNone(parse_chain_selection("all"))
        self.assertIsNone(parse_chain_selection("*"))
        self.assertEqual(parse_chain_selection("A B"), ["A", "B"])
        self.assertEqual(parse_chain_selection("A,B,A"), ["A", "B"])
        self.assertEqual(parse_chain_selection("_"), [" "])
        self.assertEqual(parse_chain_selection("blank C"), [" ", "C"])
        with self.assertRaisesRegex(CombinePDBError, "not a single-character PDB chain ID"):
            parse_chain_selection("AB")

    def test_list_pdb_chains_reports_first_appearance_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.pdb"
            source.write_text(
                "".join(
                    [
                        atom_line(1, "Y", 1),
                        atom_line(2, "X", 1),
                        atom_line(3, "Y", 2),
                        atom_line(4, " ", 1),
                    ]
                ),
                encoding="utf-8",
            )
            self.assertEqual(list_pdb_chains(source), ["Y", "X", " "])

    def test_rejects_more_than_26_chains(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            first.write_text(
                "".join(atom_line(index, chr(ord("A") + index - 1), 1) for index in range(1, 27)),
                encoding="utf-8",
            )
            second.write_text(atom_line(1, "a", 1), encoding="utf-8")
            with self.assertRaisesRegex(CombinePDBError, "26 uppercase alphabetic chain IDs"):
                combine_pdb_files([first, second], root / "combined.pdb")


if __name__ == "__main__":
    unittest.main()
