"""Comprehensive P5 regression tests for state invalidation, extraction, alignment provenance, metadata, and PDF checks."""

from __future__ import annotations

import base64
from pathlib import Path
import re
import zlib

import pytest
from streamlit.testing.v1 import AppTest

from bio3d.alignment import run_pairwise_alignment
from bio3d.config import BUNDLED_STRUCTURES, DATA_DIR, UBIQUITIN_SEQUENCE
from bio3d.reporting import generate_pdf_report, generate_text_report
from bio3d.sequence import analyze_sequence, kyte_doolittle_profile, validate_sequence
from bio3d.structural_analysis import (
    compute_contact_map,
    compute_ramachandran,
    compute_secondary_structure,
)
from bio3d.structure import parse_structure
from ui.state import set_active_structure

APP_PATH = Path(__file__).resolve().parent.parent / "app.py"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
ESMFOLD_76_PATH = FIXTURES_DIR / "esmfold_ubiquitin_76aa.pdb"
CRN_PATH = DATA_DIR / "1crn_crambin.pdb"
LYZ_PATH = DATA_DIR / "1lyz_lysozyme.pdb"


def _decompress_pdf_stream(raw: bytes) -> str:
    stripped = raw.strip(b"\r\n")
    for candidate in (stripped,):
        try:
            return zlib.decompress(candidate).decode("latin-1", errors="replace")
        except zlib.error:
            pass
    if stripped.endswith(b"~>") or b"~" in stripped[-4:]:
        try:
            a85 = stripped if stripped.endswith(b"~>") else stripped + b"~>"
            decoded = base64.a85decode(a85, adobe=True)
            return zlib.decompress(decoded).decode("latin-1", errors="replace")
        except Exception:
            pass
    return stripped.decode("latin-1", errors="replace")


def _decode_pdf_text_unicode(pdf_bytes: bytes) -> str:
    """Pure-Python ToUnicode CMap + PDF content stream extractor for ReportLab PDFs."""
    objs: dict[int, bytes] = {}
    for m in re.finditer(rb"(\d+)\s+0\s+obj\b(.*?)endobj", pdf_bytes, re.DOTALL):
        objs[int(m.group(1))] = m.group(2)

    # Map /FontResourceName -> {cid: unicode} using /ToUnicode references per font dictionary
    cmap_by_obj: dict[int, dict[int, str]] = {}
    for obj_id, body in objs.items():
        s_pos = body.find(b"stream")
        e_pos = body.rfind(b"endstream")
        if s_pos == -1 or e_pos == -1:
            continue
        text = _decompress_pdf_stream(body[s_pos + 6 : e_pos])
        if "beginbfchar" in text or "beginbfrange" in text:
            local_map: dict[int, str] = {}
            for block in re.findall(r"beginbfchar(.*?)endbfchar", text, re.DOTALL):
                for cid_hex, uni_hex in re.findall(r"<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>", block):
                    local_map[int(cid_hex, 16)] = bytes.fromhex(uni_hex).decode("utf-16-be", errors="replace")
            for block in re.findall(r"beginbfrange(.*?)endbfrange", text, re.DOTALL):
                for lo_h, hi_h, dst_h in re.findall(
                    r"<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>", block
                ):
                    lo, hi, dst = int(lo_h, 16), int(hi_h, 16), int(dst_h, 16)
                    for offset, cid in enumerate(range(lo, hi + 1)):
                        local_map[cid] = chr(dst + offset)
            cmap_by_obj[obj_id] = local_map

    # Link font objects (e.g., /ToUnicode 12 0 R) and map content-stream /F* font tags
    obj_to_cmap: dict[int, dict[int, str]] = {}
    for obj_id, body in objs.items():
        m = re.search(rb"/ToUnicode\s+(\d+)\s+0\s+R", body)
        if m and int(m.group(1)) in cmap_by_obj:
            obj_to_cmap[obj_id] = cmap_by_obj[int(m.group(1))]

    tag_to_cmap: dict[str, dict[int, str]] = {}
    for body in objs.values():
        for tag, f_obj in re.findall(rb"/(F\d+(?:\+\d+)?)\s+(\d+)\s+0\s+R", body):
            f_id = int(f_obj)
            if f_id in obj_to_cmap:
                tag_to_cmap[tag.decode("latin-1")] = obj_to_cmap[f_id]

    all_cids: dict[int, str] = {}
    for m in cmap_by_obj.values():
        all_cids.update(m)

    out_chunks: list[str] = []
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
        decoded = _decompress_pdf_stream(pdf_bytes[start:e_idx])

        if "BT" in decoded and ("Tj" in decoded or "TJ" in decoded):
            active_cmap = all_cids
            for token in re.finditer(
                r"/(F\d+(?:\+\d+)?)\s+[\d.]+\s+Tf|\(([^()\\]*(?:\\.[^()\\]*)*)\)", decoded
            ):
                if token.group(1):
                    active_cmap = tag_to_cmap.get(token.group(1), all_cids)
                elif token.group(2) is not None:
                    lit = token.group(2)
                    unesc = re.sub(r"\\([0-7]{1,3})", lambda m: chr(int(m.group(1), 8)), lit)
                    unesc = unesc.replace(r"\(", "(").replace(r"\)", ")").replace(r"\\", "\\")
                    out_chunks.append("".join(active_cmap.get(ord(ch), ch) for ch in unesc))
        pos = e_idx + 9
    return " ".join(out_chunks)


def _count_pdf_pages(pdf_bytes: bytes) -> int:
    # Match either '/Type /Pages ... /Count N' or '/Count N ... /Type /Pages' inside PDF dict objects
    counts: list[int] = []
    for body in re.findall(rb"<<(.*?)>>", pdf_bytes, re.DOTALL):
        if re.search(rb"/Type\s*/Pages\b", body):
            m = re.search(rb"/Count\s+(\d+)", body)
            if m:
                counts.append(int(m.group(1)))
    return max(counts) if counts else 0


def _synthetic_two_chain_pdb() -> str:
    """Build a valid synthetic 2-chain PDB with distinct chain sequences (Chain A: 20 aa, Chain B: 12 aa)."""
    aa1_to_3 = {
        "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
        "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
        "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
        "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    }
    seq_a = "AAAAAAAAAAALLLLLLLLL"  # 20 aa, zero aromatic, MW ~ 1756 Da
    seq_b = "DDDDKKKKEERR"          # 12 aa, charged
    lines = [
        "HEADER    SYNTHETIC TWO-CHAIN TEST              01-JAN-26   2CHN",
        "EXPDTA    X-RAY DIFFRACTION",
        "REMARK   2 RESOLUTION.    1.90 ANGSTROMS.",
    ]
    serial = 1
    for ch, seq, y_off in (("A", seq_a, 0.0), ("B", seq_b, 20.0)):
        for i, aa in enumerate(seq, start=1):
            res3 = aa1_to_3[aa]
            x = i * 3.8
            for atom_name, dx, dy, elem in (("N", -0.5, 0.0, "N"), ("CA", 0.0, 0.0, "C"), ("C", 0.8, 0.2, "C"), ("O", 1.1, 1.3, "O")):
                lines.append(
                    f"ATOM  {serial:5d}  {atom_name:<3s} {res3:>3s} {ch}{i:4d}    "
                    f"{x + dx:8.3f}{y_off + dy:8.3f}{0.0:8.3f}  1.00 12.50           {elem:>2s}  "
                )
                serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


# 1. Generate a report for ESMFold 76-aa structure
def test_01_esmfold_76aa_report():
    parsed = parse_structure(ESMFOLD_76_PATH.read_text(encoding="utf-8"), ESMFOLD_76_PATH.name, origin_hint="predicted_esmfold")
    assert parsed.is_predicted is True
    assert parsed.polymer_residue_count == 76
    assert parsed.polymer_atom_count == 601
    assert parsed.resolution_angstrom is None
    assert parsed.bfactor_stats is None
    assert parsed.plddt_stats is not None

    seq = parsed.chain_sequences["A"]
    assert len(seq) == 76
    prof = analyze_sequence(seq)
    hydro = kyte_doolittle_profile(seq)
    rama = compute_ramachandran(parsed, chain_id="A")
    ss = compute_secondary_structure(parsed, chain_id="A")
    cmap = compute_contact_map(parsed, chain_id="A", atom_mode="CA", threshold_angstrom=8.0, min_seq_sep=6)

    pdf = generate_pdf_report(
        parsed, prof, hydro, alignment=None, rama=rama, ss=ss, contact_map=cmap,
        sequence_name="ESMFold_chain_A", sequence_source="structure", report_chain="A",
    )
    txt = generate_text_report(
        parsed, prof, alignment=None, rama=rama, ss=ss, contact_map=cmap,
        sequence_name="ESMFold_chain_A", sequence_source="structure", report_chain="A",
    )
    assert pdf.startswith(b"%PDF-")
    assert "PAIRWISE ALIGNMENT [UNAVAILABLE]" in txt
    assert "No pairwise sequence alignment was performed for the current structure/report session." in txt
    assert "76 aa" in txt
    assert "[AI-PREDICTED] Mean pLDDT" in txt


# 2. Generate a report for 1CRN 46-aa structure + P3 metadata verification
def test_02_1crn_46aa_experimental_metadata_and_report():
    parsed = parse_structure(CRN_PATH.read_text(encoding="utf-8"), CRN_PATH.name, origin_hint="experimental")
    assert parsed.pdb_id == "1CRN"
    assert parsed.polymer_residue_count == 46
    assert parsed.polymer_atom_count == 327
    assert parsed.water_count == 0
    assert parsed.resolution_angstrom == pytest.approx(1.50, abs=1e-3)
    assert parsed.bfactor_stats is not None
    assert parsed.bfactor_stats["mean"] == pytest.approx(5.81, abs=0.05)
    assert parsed.bfactor_stats["occupancy_min"] == pytest.approx(1.0, abs=1e-3)
    assert parsed.bfactor_stats["occupancy_max"] == pytest.approx(1.0, abs=1e-3)

    seq = parsed.chain_sequences["A"]
    assert len(seq) == 46
    prof = analyze_sequence(seq)
    assert prof.molecular_weight_da == pytest.approx(4736.43, abs=0.1)
    assert "46" in str(prof.length)

    txt = generate_text_report(
        parsed, prof, alignment=None, sequence_name="1CRN_chain_A", sequence_source="structure", report_chain="A"
    )
    assert "46 polymer residues, 327 protein atoms" in txt
    assert "Resolution: 1.50 A" in txt
    assert "4736.43 Da" in txt
    assert "76 aa" not in txt


# 3 & 4. Switch from ESMFold 76 aa to 1CRN 46 aa WITHOUT restarting Streamlit; verify 1CRN never shows 76-aa metrics
def test_03_04_switch_esmfold_to_1crn_without_restart_no_stale_metrics():
    at = AppTest.from_file(str(APP_PATH), default_timeout=25)
    esm_text = ESMFOLD_76_PATH.read_text(encoding="utf-8")
    at.run()
    assert not at.exception

    # Switch active structure to ESMFold 76-aa in running session state
    at.session_state["structure_raw"] = esm_text
    at.session_state["structure_filename"] = ESMFOLD_76_PATH.name
    at.session_state["structure_origin_hint"] = "predicted_esmfold"
    at.run()
    assert not at.exception

    # Generate PDF for ESMFold 76-aa
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Reports & Exports").run()
    next(b for b in at.button if b.key == "btn_gen_pdf").click().run()
    assert not at.exception
    esm_pdf = at.session_state["generated_pdf_bytes"]
    assert len(esm_pdf) > 10_000

    # Now switch in the SAME session (no restart!) to bundled 1CRN via sidebar sample button
    crn_label = next(k for k in BUNDLED_STRUCTURES if "1CRN" in k)
    sel = next(s for s in at.sidebar.selectbox if s.key == "sidebar_sample_select")
    sel.set_value(crn_label).run()
    next(b for b in at.sidebar.button if b.key == "btn_load_sample").click().run()
    assert not at.exception

    # Old PDF bytes must have been invalidated immediately upon switching to 1CRN
    assert "generated_pdf_bytes" not in at.session_state

    # Check Sequence Analytics tab shows 46 aa, NOT 76 aa
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Sequence Analytics").run()
    metric_vals = [m.value for m in at.metric]
    assert "46 aa" in metric_vals
    assert "76 aa" not in metric_vals

    # Generate PDF report for 1CRN and verify zero 76-aa metrics remain
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Reports & Exports").run()
    next(b for b in at.button if b.key == "btn_gen_pdf").click().run()
    assert not at.exception
    crn_pdf = at.session_state["generated_pdf_bytes"]
    crn_text = _decode_pdf_text_unicode(crn_pdf)
    assert "46 amino acids" in crn_text
    assert "4,736.43 Da" in crn_text
    assert "76 amino acids" not in crn_text
    assert "8,564.74" not in crn_text
    assert "-0.489" not in crn_text

    # Switch to Custom FASTA with Ubiquitin while 1CRN is loaded -> mismatch warning must appear
    seq_ctrl = next(w for w in at.sidebar.segmented_control if w.key == "sequence_source_mode")
    seq_ctrl.set_value("Custom FASTA").run()
    fasta_box = next(t for t in at.sidebar.text_area if t.key == "custom_fasta_text")
    fasta_box.set_value(f">Mismatch_Test\n{UBIQUITIN_SEQUENCE}").run()
    assert any("Sequence mismatch" in w.value for w in at.sidebar.warning)


# 5. Generate a report for 1LYZ and verify 129-aa extraction, metadata & modified residues
def test_05_1lyz_129aa_extraction_metadata_and_modified_residues():
    parsed = parse_structure(LYZ_PATH.read_text(encoding="utf-8"), LYZ_PATH.name, origin_hint="experimental")
    assert parsed.pdb_id == "1LYZ"
    assert parsed.polymer_residue_count == 129
    assert parsed.polymer_atom_count == 1001
    assert parsed.water_count == 101
    assert parsed.resolution_angstrom == pytest.approx(2.00, abs=1e-3)
    assert parsed.bfactor_stats is not None
    # Distinguish CA vs all-atom B-factors and occupancy (col 55-60)
    assert parsed.bfactor_stats["mean"] == pytest.approx(6.00, abs=0.05)
    assert parsed.bfactor_stats["all_atom_mean"] == pytest.approx(6.65, abs=0.05)
    assert parsed.bfactor_stats["all_atom_max"] == pytest.approx(16.00, abs=0.05)
    assert parsed.bfactor_stats["occupancy_min"] == pytest.approx(1.00, abs=1e-3)

    seq = parsed.chain_sequences["A"]
    assert len(seq) == 129
    val = validate_sequence(seq)
    assert val.is_valid is True
    prof = analyze_sequence(val.cleaned_sequence)
    txt = generate_text_report(
        parsed, prof, alignment=None, sequence_name="1LYZ_chain_A", sequence_source="structure", report_chain="A"
    )
    assert "Length: 129 aa" in txt
    assert "129 polymer residues, 1001 protein atoms" in txt
    assert "Waters: 101" in txt
    assert "No valid standard amino-acid sequence" not in txt

    # Verify HETATM modified backbone residues (MSE/SEP/TPO/PTR) map cleanly to parent AA
    mod_pdb = """HEADER    MODIFIED RESIDUES TEST                  01-JAN-26   MODT
EXPDTA    X-RAY DIFFRACTION
REMARK   2 RESOLUTION.    1.60 ANGSTROMS.
HETATM    1  N   MSE A   1       0.000   0.000   0.000  1.00 10.00           N
HETATM    2  CA  MSE A   1       1.450   0.000   0.000  1.00 10.00           C
HETATM    3  C   MSE A   1       2.050   1.350   0.000  1.00 10.00           C
HETATM    4  N   SEP A   2       3.400   1.400   0.000  1.00 10.00           N
HETATM    5  CA  SEP A   2       4.100   2.650   0.000  1.00 10.00           C
HETATM    6  C   SEP A   2       5.600   2.600   0.000  1.00 10.00           C
HETATM    7  N   TPO A   3       6.200   3.800   0.000  1.00 10.00           N
HETATM    8  CA  TPO A   3       7.650   3.900   0.000  1.00 10.00           C
HETATM    9  C   TPO A   3       8.200   5.250   0.000  1.00 10.00           C
HETATM   10  N   PTR A   4       9.550   5.400   0.000  1.00 10.00           N
HETATM   11  CA  PTR A   4      10.200   6.700   0.000  1.00 10.00           C
HETATM   12  C   PTR A   4      11.700   6.700   0.000  1.00 10.00           C
END
"""
    p_mod = parse_structure(mod_pdb, "mod_test.pdb", origin_hint="experimental")
    assert p_mod.chain_sequences["A"] == "MSTY"
    assert len(p_mod.modified_residues_mapped["A"]) == 4


# 6. Generate an alignment report, then switch structures and ensure alignment state does not leak
def test_06_alignment_provenance_and_zero_leakage_on_switch():
    at = AppTest.from_file(str(APP_PATH), default_timeout=25)
    at.run()

    # Visit Pairwise Alignment and explicitly run HBA vs HBB alignment
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Pairwise Alignment").run()
    next(b for b in at.button if b.key == "btn_run_alignment").click().run()
    assert not at.exception
    assert at.session_state["alignment_active_in_session"] is True

    # Generate PDF report for 1UBQ with this alignment
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Reports & Exports").run()
    next(b for b in at.button if b.key == "btn_gen_pdf").click().run()
    pdf_with_aln = _decode_pdf_text_unicode(at.session_state["generated_pdf_bytes"])
    assert "P69905" in pdf_with_aln and "P68871" in pdf_with_aln
    assert "142 aa" in pdf_with_aln and "147 aa" in pdf_with_aln
    assert "43.62%" in pdf_with_aln

    # Switch structure to 1CRN -> alignment MUST NOT leak!
    crn_label = next(k for k in BUNDLED_STRUCTURES if "1CRN" in k)
    sel = next(s for s in at.sidebar.selectbox if s.key == "sidebar_sample_select")
    sel.set_value(crn_label).run()
    next(b for b in at.sidebar.button if b.key == "btn_load_sample").click().run()
    assert not at.session_state.get("alignment_active_in_session", False)

    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Reports & Exports").run()
    next(b for b in at.button if b.key == "btn_gen_pdf").click().run()
    pdf_after_switch = _decode_pdf_text_unicode(at.session_state["generated_pdf_bytes"])
    assert "[UNAVAILABLE]" in pdf_after_switch
    assert "43.62%" not in pdf_after_switch
    assert "P69905" not in pdf_after_switch
    assert "P68871" not in pdf_after_switch


# 7. Generate reports after changing chains
def test_07_report_updates_when_changing_chains():
    two_chain_pdb = _synthetic_two_chain_pdb()
    at = AppTest.from_file(str(APP_PATH), default_timeout=25)
    at.run()
    at.session_state["structure_raw"] = two_chain_pdb
    at.session_state["structure_filename"] = "two_chains.pdb"
    at.session_state["structure_origin_hint"] = "experimental"
    at.run()
    assert not at.exception

    # Default Chain A (20 aa)
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Reports & Exports").run()
    next(b for b in at.button if b.key == "btn_gen_pdf").click().run()
    pdf_a = _decode_pdf_text_unicode(at.session_state["generated_pdf_bytes"])
    assert "20 amino acids" in pdf_a
    assert "Chain A" in pdf_a

    # Switch sidebar chain selector to Chain B (12 aa) -> PDF bytes must invalidate & update to 12 aa
    ch_sel = next(s for s in at.sidebar.selectbox if s.key == "sidebar_structure_chain")
    ch_sel.set_value("B").run()
    assert "generated_pdf_bytes" not in at.session_state
    next(b for b in at.button if b.key == "btn_gen_pdf").click().run()
    pdf_b = _decode_pdf_text_unicode(at.session_state["generated_pdf_bytes"])
    assert "12 amino acids" in pdf_b
    assert "Chain B" in pdf_b
    assert "20 amino acids" not in pdf_b


# 8. Test unavailable sequences and missing metadata
def test_08_unavailable_sequences_and_missing_metadata():
    unk_pdb = """ATOM      1  N   UNK A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  UNK A   1       1.450   0.000   0.000  1.00  0.00           C
ATOM      3  C   UNK A   1       2.050   1.350   0.000  1.00  0.00           C
END
"""
    parsed = parse_structure(unk_pdb, "unknown_only.pdb", origin_hint=None)
    assert parsed.pdb_id is None
    assert parsed.experimental_method is None
    assert parsed.resolution_angstrom is None
    assert parsed.unmapped_polymer_residues["A"]
    assert "UNK" in parsed.chain_extraction_notes["A"]
    val = validate_sequence(parsed.chain_sequences["A"])
    assert val.is_valid is False

    reason = parsed.chain_extraction_notes["A"]
    txt = generate_text_report(
        parsed, None, sequence_unavailable_reason=reason, sequence_source="structure", report_chain="A"
    )
    assert "PDB ID: Unavailable" in txt
    assert "Method: Unavailable" in txt
    assert "Resolution: Unavailable" in txt
    assert "SEQUENCE METRICS [UNAVAILABLE]" in txt
    assert "UNK" in txt


# 9. Validate PDF text extraction, page count and Unicode characters across generated reports
def test_09_pdf_text_extraction_page_count_and_unicode():
    for path, hint, expected_res_text in (
        (ESMFOLD_76_PATH, "predicted_esmfold", "76 amino acids"),
        (CRN_PATH, "experimental", "46 amino acids"),
        (LYZ_PATH, "experimental", "129 amino acids"),
    ):
        parsed = parse_structure(path.read_text(encoding="utf-8"), path.name, origin_hint=hint)
        seq = parsed.chain_sequences["A"]
        prof = analyze_sequence(seq)
        hydro = kyte_doolittle_profile(seq)
        rama = compute_ramachandran(parsed, chain_id="A")
        ss = compute_secondary_structure(parsed, chain_id="A")
        cmap = compute_contact_map(parsed, chain_id="A", atom_mode="CA", threshold_angstrom=8.0, min_seq_sep=6)
        pdf = generate_pdf_report(
            parsed, prof, hydro, alignment=None, rama=rama, ss=ss, contact_map=cmap,
            sequence_name=f"{parsed.pdb_id or path.stem}_chain_A",
            sequence_source="structure",
            report_chain="A",
        )
        pages = _count_pdf_pages(pdf)
        assert 3 <= pages <= 6
        assert b"DejaVuSans" in pdf
        decoded = _decode_pdf_text_unicode(pdf)
        assert "\ufffd" not in decoded
        assert expected_res_text in decoded
        assert re.search(r"M\s*-1\s*cm\s*-1", decoded) is not None
        assert "Å" in decoded
        assert "[UNAVAILABLE]" in decoded
