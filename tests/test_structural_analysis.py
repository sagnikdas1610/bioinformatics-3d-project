"""Scientific validation tests for Ramachandran, Secondary Structure, and Contact Maps."""

from __future__ import annotations

from bio3d.config import DATA_DIR
from bio3d.structural_analysis import (
    compute_contact_map,
    compute_ramachandran,
    compute_secondary_structure,
)
from bio3d.structure import parse_structure


def test_1ubq_ramachandran_and_secondary_structure():
    text = (DATA_DIR / "1ubq_ubiquitin.pdb").read_text(encoding="utf-8")
    parsed = parse_structure(text, filename="1ubq_ubiquitin.pdb")

    rama = compute_ramachandran(parsed, chain_id="A")
    assert rama.total_evaluated == 74  # 76 residues minus N- and C-termini
    assert "MolProbity Top8000" in rama.reference_source
    assert any(r.rama_class == "Ile/Val" for r in rama.residues)
    assert rama.favored_pct >= 94.0
    assert rama.outlier_pct <= 2.0

    ss = compute_secondary_structure(parsed, chain_id="A")
    assert ss.total_residues == 76
    # Human ubiquitin is an alpha+beta fold (~20% helix, ~30-38% beta sheet)
    assert 15.0 <= ss.helix_pct <= 32.0
    assert 22.0 <= ss.sheet_pct <= 45.0
    assert ss.header_available is True


def test_1crn_ramachandran_and_secondary_structure():
    text = (DATA_DIR / "1crn_crambin.pdb").read_text(encoding="utf-8")
    parsed = parse_structure(text, filename="1crn_crambin.pdb")
    rama = compute_ramachandran(parsed)
    assert rama.total_evaluated == 44
    assert rama.favored_pct >= 90.0

    ss = compute_secondary_structure(parsed)
    assert ss.total_residues == 46
    assert ss.helix_pct >= 30.0  # Crambin has two prominent alpha-helices


def test_contact_map_glycine_cb_fallback_to_ca():
    text = (DATA_DIR / "1ubq_ubiquitin.pdb").read_text(encoding="utf-8")
    parsed = parse_structure(text, filename="1ubq_ubiquitin.pdb")
    gly_count = parsed.chain_sequences["A"].count("G")
    assert gly_count > 0

    cmap_cb = compute_contact_map(parsed, chain_id="A", atom_mode="CB", threshold_angstrom=8.0, min_seq_sep=6)
    assert cmap_cb.residue_count == 76
    assert cmap_cb.glycine_fallback_count == gly_count
    assert cmap_cb.total_contacts > 0
    assert cmap_cb.long_range_contacts > 0

    cmap_ca = compute_contact_map(parsed, chain_id="A", atom_mode="CA", threshold_angstrom=8.0, min_seq_sep=6)
    assert cmap_ca.glycine_fallback_count == 0
    assert cmap_ca.residue_count == 76
