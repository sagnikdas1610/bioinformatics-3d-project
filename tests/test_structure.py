"""Tests for PDB/mmCIF structure parsing, model-0 discipline, provenance, and missing residues."""

from __future__ import annotations

from pathlib import Path

from bio3d.config import DATA_DIR, UBIQUITIN_SEQUENCE
from bio3d.structure import parse_structure


def test_parse_1ubq_benchmark():
    text = (DATA_DIR / "1ubq_ubiquitin.pdb").read_text(encoding="utf-8")
    parsed = parse_structure(text, filename="1ubq_ubiquitin.pdb")
    assert parsed.file_format == "pdb"
    assert parsed.pdb_id == "1UBQ"
    assert parsed.is_predicted is False
    assert parsed.chains == ["A"]
    assert parsed.polymer_residue_count == 76
    assert parsed.chain_sequences["A"] == UBIQUITIN_SEQUENCE
    assert parsed.resolution_angstrom == 1.8
    assert parsed.water_count > 0
    assert parsed.bfactor_stats is not None
    assert parsed.plddt_stats is None


def test_parse_1crn_benchmark():
    text = (DATA_DIR / "1crn_crambin.pdb").read_text(encoding="utf-8")
    parsed = parse_structure(text, filename="1crn_crambin.pdb")
    assert parsed.pdb_id == "1CRN"
    assert parsed.chains == ["A"]
    assert parsed.polymer_residue_count == 46
    assert parsed.resolution_angstrom == 1.5
    assert parsed.bfactor_stats is not None
    assert parsed.plddt_stats is None


def test_multimodel_and_altloc_deduplication():
    pdb_sample = """HEADER    NMR TEST                                01-JAN-24   1TST
EXPDTA    SOLUTION NMR
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.70 10.00           N
ATOM      2  N  BALA A   1       0.100   0.000   0.000  0.30 20.00           N
ATOM      3  CA  ALA A   1       1.458   0.000   0.000  1.00 12.00           C
HETATM    4  O   HOH A 101       5.000   5.000   5.000  1.00 30.00           O
ENDMDL
MODEL        2
ATOM      1  N   ALA A   1       9.000   9.000   9.000  1.00 15.00           N
ATOM      2  CA  ALA A   1       9.458   9.000   9.000  1.00 15.00           C
ENDMDL
"""
    parsed = parse_structure(pdb_sample, filename="nmr.pdb")
    assert parsed.model_count == 2
    assert parsed.polymer_residue_count == 1
    assert parsed.polymer_atom_count == 2  # altLoc B deduplicated, Model 2 excluded
    assert parsed.water_count == 1
    assert parsed.hetero_ligand_count == 0
    assert parsed.polymer_residues["A"][0].atoms["N"].occupancy == 0.70


def test_mmcif_missing_residues_and_plddt_separation():
    cif_sample = """data_AF_TEST
_entry.id AF_TEST
_ma_model_list.ordinal_id 1
loop_
_pdbx_unobs_or_zero_occ_residues.PDB_model_num
_pdbx_unobs_or_zero_occ_residues.auth_asym_id
_pdbx_unobs_or_zero_occ_residues.auth_comp_id
_pdbx_unobs_or_zero_occ_residues.auth_seq_id
1 A GLY 1
1 A SER 2
loop_
_atom_site.group_PDB
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.type_symbol
_atom_site.pdbx_PDB_model_num
ATOM CA . ALA A 3 ? 0.0 0.0 0.0 1.0 0.92 C 1
ATOM CA . GLY A 4 ? 3.8 0.0 0.0 1.0 0.78 C 1
"""
    parsed = parse_structure(cif_sample, filename="test.cif")
    assert parsed.file_format == "cif"
    assert parsed.is_predicted is True
    assert parsed.bfactor_stats is None
    assert parsed.plddt_stats is not None
    # 0-1 ESMFold-like scale normalized to 0-100
    assert parsed.plddt_stats["mean"] == 85.0
    assert len(parsed.missing_residues) == 2
    assert parsed.missing_residue_source == "mmCIF _pdbx_unobs_or_zero_occ_residues"


def test_pdb_remark465_missing_residues():
    pdb_sample = """HEADER    TEST                                    01-JAN-24   1MIS
EXPDTA    X-RAY DIFFRACTION
REMARK 465
REMARK 465 MISSING RESIDUES
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)
REMARK 465   M RES C SSSEQI
REMARK 465     MET A     1
REMARK 465     GLN A     2
ATOM      1  CA  ILE A   3       0.000   0.000   0.000  1.00 18.00           C
"""
    parsed = parse_structure(pdb_sample, filename="missing.pdb")
    assert len(parsed.missing_residues) == 2
    assert parsed.missing_residue_source == "REMARK 465"
    assert parsed.missing_residues[0]["resname"] == "MET"
    assert parsed.missing_residues[0]["resseq"] == 1
