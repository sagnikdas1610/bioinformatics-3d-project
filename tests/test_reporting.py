"""Tests for PDF layout, Unicode typography, metadata accuracy, and plain-text reporting."""

from __future__ import annotations

import zlib

import pytest
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph

from bio3d.alignment import run_pairwise_alignment
from bio3d.config import DATA_DIR, HBA_HUMAN_SEQ, HBB_HUMAN_SEQ, UBIQUITIN_SEQUENCE
from bio3d.reporting import (
    USABLE_WIDTH,
    _build_pdf_styles,
    _styled_table,
    generate_pdf_report,
    generate_text_report,
)
from bio3d.sequence import analyze_sequence, kyte_doolittle_profile
from bio3d.structural_analysis import (
    compute_contact_map,
    compute_ramachandran,
    compute_secondary_structure,
)
from bio3d.structure import parse_structure


def _extract_pdf_streams_text(pdf_bytes: bytes) -> str:
    """Decompress FlateDecode PDF content streams for text assertions without external deps."""
    chunks: list[str] = []
    pos = 0
    while True:
        s_idx = pdf_bytes.find(b"stream", pos)
        if s_idx == -1:
            break
        start = s_idx + 6
        if pdf_bytes[start : start + 2] == b"\r\n":
            start += 2
        elif pdf_bytes[start : start + 1] in (b"\n", b"\r"):
            start += 1
        e_idx = pdf_bytes.find(b"endstream", start)
        if e_idx == -1:
            break
        raw = pdf_bytes[start:e_idx].strip(b"\r\n")
        try:
            chunks.append(zlib.decompress(raw).decode("latin-1", errors="replace"))
        except zlib.error:
            chunks.append(raw.decode("latin-1", errors="replace"))
        pos = e_idx + 9
    return "\n".join(chunks)


def test_styled_table_wraps_all_cells_within_usable_width():
    styles = _build_pdf_styles()
    widths = [1.85 * inch, 3.95 * inch, 1.40 * inch]
    rows = [
        ["Field", "Value", "Provenance"],
        [
            "Missing Residues Check",
            "No formal missing-residue records (REMARK 465 / _pdbx_unobs_or_zero_occ_residues) and no internal backbone numbering breaks detected.",
            "Experimental (Measured)",
        ],
        [
            "pLDDT Confidence (0-100)",
            "Mean 98.1 | Median 98.8 | >=90: 95.8% | <50: 0.7% (raw scale: 0-100)",
            "AI-Predicted",
        ],
    ]
    tbl = _styled_table(rows, widths, styles)
    for row in tbl._cellvalues:
        for cell in row:
            assert isinstance(cell, Paragraph)
    w, h = tbl.wrap(USABLE_WIDTH, 10 * inch)
    assert w <= USABLE_WIDTH + 1e-3
    assert h > 0


def test_generate_full_pdf_and_text_report_metadata():
    parsed = parse_structure((DATA_DIR / "1ubq_ubiquitin.pdb").read_text(encoding="utf-8"), "1ubq_ubiquitin.pdb")
    seq_profile = analyze_sequence(UBIQUITIN_SEQUENCE)
    hydro = kyte_doolittle_profile(UBIQUITIN_SEQUENCE, window=9)
    aln = run_pairwise_alignment(
        HBA_HUMAN_SEQ,
        HBB_HUMAN_SEQ,
        mode="global",
        seq1_name="HBA_HUMAN (Hemoglobin subunit alpha)",
        seq1_accession="P69905",
        seq2_name="HBB_HUMAN (Hemoglobin subunit beta)",
        seq2_accession="P68871",
    )
    rama = compute_ramachandran(parsed, chain_id="A")
    ss = compute_secondary_structure(parsed, chain_id="A")
    cmap = compute_contact_map(parsed, chain_id="A", atom_mode="CB", threshold_angstrom=8.0, min_seq_sep=6)

    pdf_bytes = generate_pdf_report(
        parsed,
        seq_profile,
        hydro,
        aln,
        rama,
        ss,
        cmap,
        sequence_name="1UBQ_chain_A",
    )
    assert pdf_bytes.startswith(b"%PDF-")
    assert len(pdf_bytes) > 15_000
    # Ensure Unicode TrueType subset fonts (DejaVuSans / BioSans) are embedded rather than WinAnsi only
    assert b"DejaVuSans" in pdf_bytes or b"/F2+0" in pdf_bytes
    assert "\ufffd".encode("utf-8") not in pdf_bytes

    txt = generate_text_report(parsed, seq_profile, aln, rama, ss, cmap, sequence_name="1UBQ_chain_A")
    assert "1UBQ" in txt
    assert "B-factor" in txt
    assert "P69905" in txt and "P68871" in txt
    assert "142 aa" in txt and "147 aa" in txt
    assert "M^-1 cm^-1" in txt


def test_predicted_structure_pdf_never_fabricates_resolution():
    cif_pred = """data_AF_TEST
_entry.id AF-TEST-F1
_ma_model_list.ordinal_id 1
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
ATOM N . ALA A 1 ? 0.0 0.0 0.0 1.0 91.5 N 1
ATOM CA . ALA A 1 ? 1.4 0.0 0.0 1.0 91.5 C 1
ATOM C . ALA A 1 ? 2.0 1.4 0.0 1.0 91.5 C 1
"""
    parsed = parse_structure(cif_pred, "AF-TEST-F1.cif", origin_hint="predicted_alphafold")
    assert parsed.is_predicted is True
    assert parsed.resolution_angstrom is None
    pdf_bytes = generate_pdf_report(parsed, None)
    assert pdf_bytes.startswith(b"%PDF-")
    txt = generate_text_report(parsed, None)
    assert "Resolution: Unavailable" in txt
    assert "[AI-PREDICTED] Mean pLDDT" in txt


def test_generate_partial_missing_data_pdf():
    pdf_bytes = generate_pdf_report(parsed=None, seq_profile=None)
    assert pdf_bytes.startswith(b"%PDF-")
    assert len(pdf_bytes) > 1_000
