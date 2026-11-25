"""Tests for PDF and plain-text scientific report generation."""

from __future__ import annotations

from bio3d.alignment import run_pairwise_alignment
from bio3d.config import DATA_DIR, HBA_HUMAN_SEQ, HBB_HUMAN_SEQ, UBIQUITIN_SEQUENCE
from bio3d.reporting import generate_pdf_report, generate_text_report
from bio3d.sequence import analyze_sequence, kyte_doolittle_profile
from bio3d.structural_analysis import (
    compute_contact_map,
    compute_ramachandran,
    compute_secondary_structure,
)
from bio3d.structure import parse_structure


def test_generate_full_pdf_and_text_report():
    parsed = parse_structure((DATA_DIR / "1ubq_ubiquitin.pdb").read_text(encoding="utf-8"), "1ubq_ubiquitin.pdb")
    seq_profile = analyze_sequence(UBIQUITIN_SEQUENCE)
    hydro = kyte_doolittle_profile(UBIQUITIN_SEQUENCE, window=9)
    aln = run_pairwise_alignment(HBA_HUMAN_SEQ, HBB_HUMAN_SEQ, mode="global")
    rama = compute_ramachandran(parsed, chain_id="A")
    ss = compute_secondary_structure(parsed, chain_id="A")
    cmap = compute_contact_map(parsed, chain_id="A", atom_mode="CB", threshold_angstrom=8.0, min_seq_sep=6)

    pdf_bytes = generate_pdf_report(parsed, seq_profile, hydro, aln, rama, ss, cmap)
    assert pdf_bytes.startswith(b"%PDF-")
    assert len(pdf_bytes) > 15_000

    txt = generate_text_report(parsed, seq_profile, aln, rama, ss, cmap)
    assert "1UBQ" in txt
    assert "Mean B-factor" in txt


def test_generate_partial_missing_data_pdf():
    pdf_bytes = generate_pdf_report(parsed=None, seq_profile=None)
    assert pdf_bytes.startswith(b"%PDF-")
    assert len(pdf_bytes) > 1_000
