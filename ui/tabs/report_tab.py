"""Scientific PDF & Data Export Center tab."""

from __future__ import annotations

import hashlib

import streamlit as st

from bio3d.alignment import (
    AlignmentResult,
    format_alignment_blocks,
    parse_sequence_identity,
    run_pairwise_alignment,
)
from bio3d.reporting import generate_pdf_report, generate_text_report
from bio3d.sequence import SequenceProfile, kyte_doolittle_profile, parse_fasta_records, validate_sequence
from bio3d.structure import ParsedStructure
from ui.state import ActiveSequenceContext, compute_sequence_fingerprint


def _recompute_active_alignment(current_struct_fp: str) -> AlignmentResult | None:
    """Recompute pairwise alignment only when explicitly run for the current structure session."""
    if not st.session_state.get("alignment_active_in_session", False):
        return None
    if st.session_state.get("alignment_structure_fp") != current_struct_fp:
        return None

    snap = st.session_state.get("alignment_inputs_snapshot") or {}
    raw_s1 = st.session_state.get("aln_seq1")
    raw_s2 = st.session_state.get("aln_seq2")

    if raw_s1 and raw_s2:
        recs1 = parse_fasta_records(raw_s1)
        recs2 = parse_fasta_records(raw_s2)
        if not recs1 or not recs2:
            return None
        name1, acc1 = parse_sequence_identity(recs1[0].header, fallback_name="Sequence 1")
        name2, acc2 = parse_sequence_identity(recs2[0].header, fallback_name="Sequence 2")
        seq1, seq2 = recs1[0].sequence, recs2[0].sequence
        mode = st.session_state.get("aln_mode") or snap.get("mode", "global")
        matrix = st.session_state.get("aln_matrix") or snap.get("matrix_name", "BLOSUM62")
        open_gap = float(st.session_state.get("aln_open_gap", snap.get("open_gap", -10.0)))
        extend_gap = float(st.session_state.get("aln_extend_gap", snap.get("extend_gap", -0.5)))
    elif snap.get("seq1") and snap.get("seq2"):
        seq1, seq2 = snap["seq1"], snap["seq2"]
        name1, acc1 = snap.get("seq1_name", "Sequence 1"), snap.get("seq1_accession")
        name2, acc2 = snap.get("seq2_name", "Sequence 2"), snap.get("seq2_accession")
        mode = snap.get("mode", "global")
        matrix = snap.get("matrix_name", "BLOSUM62")
        open_gap = float(snap.get("open_gap", -10.0))
        extend_gap = float(snap.get("extend_gap", -0.5))
    else:
        return None

    try:
        return run_pairwise_alignment(
            seq1=seq1,
            seq2=seq2,
            mode=mode,
            matrix_name=matrix,
            open_gap_score=open_gap,
            extend_gap_score=extend_gap,
            seq1_name=name1,
            seq1_accession=acc1,
            seq2_name=name2,
            seq2_accession=acc2,
        )
    except ValueError:
        return None


def render_report_tab(
    parsed: ParsedStructure | None,
    seq_ctx: ActiveSequenceContext,
    current_struct_fp: str,
    origin_hint: str | None,
    cached_analyze_fn,
    cached_rama_fn,
    cached_ss_fn,
    cached_cmap_fn,
) -> None:
    st.markdown("Export publication-style PDF reports and reproducible analysis artifacts.")

    if seq_ctx.mismatch_warning:
        st.warning(seq_ctx.mismatch_warning, icon=":material/warning:")

    val = validate_sequence(seq_ctx.sequence)
    if val.is_valid:
        seq_hash = compute_sequence_fingerprint(val.cleaned_sequence)
        seq_profile: SequenceProfile | None = cached_analyze_fn(seq_hash, val.cleaned_sequence)
        hydro = kyte_doolittle_profile(val.cleaned_sequence, window=9)
        seq_unavailable_reason = None
    else:
        seq_profile = None
        hydro = None
        reason_parts = [p for p in (seq_ctx.extraction_note, val.message) if p]
        seq_unavailable_reason = " ".join(reason_parts) or "No valid standard amino-acid sequence was available."

    fresh_alignment = _recompute_active_alignment(current_struct_fp)
    include_alignment = st.checkbox(
        "Include active pairwise alignment in report",
        value=fresh_alignment is not None,
        disabled=fresh_alignment is None,
        help="Includes alignment ONLY when explicitly executed for the active structure session.",
        key="chk_include_alignment",
    )
    alignment = fresh_alignment if include_alignment else None

    rama = None
    ss = None
    cmap = None
    report_chain: str | None = None
    if parsed is not None:
        if seq_ctx.chain_id and seq_ctx.chain_id in parsed.chains:
            report_chain = seq_ctx.chain_id
        elif parsed.chains:
            report_chain = parsed.chains[0]
        else:
            report_chain = "All"

        cm_cfg = st.session_state.get("active_cm_params", {})
        use_cm_cfg = cm_cfg.get("structure_fp") == current_struct_fp
        cm_chain = report_chain if report_chain != "All" else parsed.chains[0]
        cm_mode = cm_cfg.get("atom_mode", "CA") if use_cm_cfg else "CA"
        cm_cutoff = float(cm_cfg.get("cutoff", 8.0)) if use_cm_cfg else 8.0
        cm_sep = int(cm_cfg.get("min_seq_sep", 6)) if use_cm_cfg else 6

        rama_scope = report_chain if (report_chain and len(parsed.chains) > 1) else "All"
        rama = cached_rama_fn(current_struct_fp, parsed.raw_text, parsed.filename, origin_hint, rama_scope)
        ss = cached_ss_fn(current_struct_fp, parsed.raw_text, parsed.filename, origin_hint, rama_scope)
        cmap = cached_cmap_fn(
            current_struct_fp,
            parsed.raw_text,
            parsed.filename,
            origin_hint,
            cm_chain,
            cm_mode,
            cm_cutoff,
            cm_sep,
        )

    # Checklist showing exactly what is available vs unavailable
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("##### Report contents status")
        st.write(
            "Structure & Quality:",
            f"Included (`{parsed.filename}`, Chain `{report_chain}`)" if parsed else "Unavailable",
        )
        if seq_profile:
            src_txt = "structure-derived" if seq_ctx.source == "structure" else "custom FASTA"
            st.write("Sequence Profile:", f"Included (`{seq_ctx.label}`, {seq_profile.length} aa, {src_txt})")
        else:
            st.write("Sequence Profile:", f"Unavailable ({seq_unavailable_reason})")
        if alignment is not None:
            st.write(
                "Pairwise Alignment:",
                f"Included (`{alignment.seq1_name}` vs `{alignment.seq2_name}`, {alignment.mode}, {alignment.matrix_name})",
            )
        else:
            st.write("Pairwise Alignment:", "Unavailable (not run for current structure session)")
    with c2:
        st.markdown("##### Structural analytics included")
        st.write("Ramachandran Plot:", f"Included ({rama.favored_pct:.1f}% favored)" if rama else "Unavailable")
        st.write("Secondary Structure:", f"Included ({ss.helix_pct:.1f}% H / {ss.sheet_pct:.1f}% E)" if ss else "Unavailable")
        st.write(
            "Contact Map:",
            f"Included (Chain {cmap.chain_id}, {cmap.atom_mode} ≤ {cmap.threshold_angstrom:.1f} Å, {cmap.total_contacts} contacts)"
            if cmap
            else "Unavailable",
        )

    st.divider()

    report_title = st.text_input(
        "Custom PDF report title (optional)",
        value="",
        placeholder="BioFold Workbench — Structural & Sequence Analysis Report",
        key="report_custom_title",
    )

    # Invalidate previously generated PDF bytes if ANY structure, chain, sequence, alignment, or title changes
    fp_parts = [
        current_struct_fp,
        str(report_chain or ""),
        seq_ctx.source,
        val.cleaned_sequence if val.is_valid else "",
        seq_ctx.label,
        f"{alignment.seq1_name}:{alignment.seq2_name}:{alignment.score:.2f}:{alignment.alignment_length}" if alignment else "NONE",
        f"{cmap.chain_id}:{cmap.atom_mode}:{cmap.threshold_angstrom}:{cmap.min_seq_sep}" if cmap else "NONE",
        report_title.strip(),
    ]
    current_fp = hashlib.sha256("|".join(fp_parts).encode("utf-8")).hexdigest()
    if st.session_state.get("generated_pdf_fp") != current_fp:
        st.session_state.pop("generated_pdf_bytes", None)

    if st.button("Generate PDF scientific report", type="primary", key="btn_gen_pdf"):
        with st.spinner("Rendering figures and compiling scientific PDF..."):
            pdf_bytes = generate_pdf_report(
                parsed=parsed,
                seq_profile=seq_profile,
                hydropathy=hydro,
                alignment=alignment,
                rama=rama,
                ss=ss,
                contact_map=cmap,
                report_title=report_title.strip() or None,
                sequence_name=seq_ctx.label,
                sequence_source=seq_ctx.source,
                sequence_mismatch_warning=seq_ctx.mismatch_warning,
                sequence_unavailable_reason=seq_unavailable_reason,
                report_chain=report_chain,
            )
            st.session_state["generated_pdf_bytes"] = pdf_bytes
            st.session_state["generated_pdf_fp"] = current_fp

    pdf_bytes = st.session_state.get("generated_pdf_bytes")
    d1, d2, d3, d4 = st.columns(4)
    if pdf_bytes:
        d1.download_button(
            "Download PDF report",
            data=pdf_bytes,
            file_name="biofold_scientific_report.pdf",
            mime="application/pdf",
            width="stretch",
        )
    else:
        d1.caption("Click **Generate PDF scientific report** to build the downloadable PDF.")

    txt_report = generate_text_report(
        parsed,
        seq_profile,
        alignment,
        rama,
        ss,
        cmap,
        sequence_name=seq_ctx.label,
        sequence_source=seq_ctx.source,
        sequence_mismatch_warning=seq_ctx.mismatch_warning,
        sequence_unavailable_reason=seq_unavailable_reason,
        report_chain=report_chain,
    )
    d2.download_button(
        "Download summary (.txt)",
        data=txt_report,
        file_name="biofold_analysis_report.txt",
        mime="text/plain",
        width="stretch",
    )

    if seq_profile is not None:
        wrapped = "\n".join(
            seq_profile.sequence[i : i + 70] for i in range(0, len(seq_profile.sequence), 70)
        )
        fasta_out = f">{seq_ctx.label} length={seq_profile.length}\n{wrapped}\n"
        d3.download_button(
            "Download cleaned FASTA",
            data=fasta_out,
            file_name="cleaned_sequence.fasta",
            mime="text/plain",
            width="stretch",
        )

    if parsed is not None:
        d4.download_button(
            f"Download {parsed.file_format.upper()}",
            data=parsed.raw_text,
            file_name=parsed.filename,
            mime="chemical/x-pdb" if parsed.file_format == "pdb" else "chemical/x-cif",
            width="stretch",
        )

    if alignment is not None:
        st.download_button(
            "Download alignment (.txt)",
            data=format_alignment_blocks(alignment, width=60),
            file_name="pairwise_alignment.txt",
            mime="text/plain",
        )
