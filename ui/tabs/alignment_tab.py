"""Pairwise Sequence Alignment tab with explicit session provenance."""

from __future__ import annotations

import streamlit as st

from bio3d.alignment import (
    SUPPORTED_MATRICES,
    AlignmentResult,
    format_alignment_blocks,
    format_alignment_html,
    parse_sequence_identity,
    run_pairwise_alignment,
)
from bio3d.config import HBA_HUMAN_SEQ, HBB_HUMAN_SEQ
from bio3d.sequence import parse_fasta_records
from ui.state import ActiveSequenceContext

DEFAULT_HBA_FASTA = f">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\n{HBA_HUMAN_SEQ}"
DEFAULT_HBB_FASTA = f">sp|P68871|HBB_HUMAN Hemoglobin subunit beta\n{HBB_HUMAN_SEQ}"


def _extract_first_record(raw_text: str, fallback_name: str) -> tuple[str, str, str | None]:
    records = parse_fasta_records(raw_text)
    if not records:
        return "", fallback_name, None
    name, acc = parse_sequence_identity(records[0].header, fallback_name=fallback_name)
    return records[0].sequence, name, acc


def render_alignment_tab(
    seq_ctx: ActiveSequenceContext,
    current_struct_fp: str,
) -> AlignmentResult | None:
    st.markdown("Compare two amino-acid sequences using global (Needleman–Wunsch) or local (Smith–Waterman) alignment.")

    st.session_state.setdefault("aln_seq1", DEFAULT_HBA_FASTA)
    st.session_state.setdefault("aln_seq2", DEFAULT_HBB_FASTA)
    st.session_state.setdefault("aln_mode", "global")
    st.session_state.setdefault("aln_matrix", "BLOSUM62")
    st.session_state.setdefault("aln_open_gap", -10.0)
    st.session_state.setdefault("aln_extend_gap", -0.5)

    active_sequence = seq_ctx.sequence
    preset_col1, preset_col2, _ = st.columns([1.2, 1.6, 1.8])
    if preset_col1.button("Load HBA vs HBB demo", width="stretch", key="btn_aln_demo"):
        st.session_state["aln_seq1"] = DEFAULT_HBA_FASTA
        st.session_state["aln_seq2"] = DEFAULT_HBB_FASTA
        st.rerun()
    if active_sequence and preset_col2.button(
        f"Use active sequence ({seq_ctx.label}) as Seq 1",
        width="stretch",
        key="btn_aln_use_active",
    ):
        header = seq_ctx.label.strip() or "Active_sequence"
        st.session_state["aln_seq1"] = f">{header}\n{active_sequence}"
        st.rerun()

    col_s1, col_s2 = st.columns(2, gap="medium")
    with col_s1:
        raw_s1 = st.text_area("Sequence 1 (FASTA or raw)", key="aln_seq1", height=130)
    with col_s2:
        raw_s2 = st.text_area("Sequence 2 (FASTA or raw)", key="aln_seq2", height=130)

    p1, p2, p3, p4 = st.columns([1.1, 1.1, 1.0, 1.0])
    mode_choice = p1.segmented_control("Alignment mode", ["global", "local"], key="aln_mode") or "global"
    matrix_name = p2.selectbox("Substitution matrix", list(SUPPORTED_MATRICES), key="aln_matrix")
    open_gap = p3.number_input("Gap open penalty", min_value=-30.0, max_value=0.0, step=0.5, key="aln_open_gap")
    extend_gap = p4.number_input("Gap extend penalty", min_value=-10.0, max_value=0.0, step=0.5, key="aln_extend_gap")

    seq1, name1, acc1 = _extract_first_record(raw_s1, "Sequence 1")
    seq2, name2, acc2 = _extract_first_record(raw_s2, "Sequence 2")

    act1, act2, _ = st.columns([1.3, 1.1, 2.2])
    run_clicked = act1.button(
        "Run Pairwise Alignment",
        type="primary",
        width="stretch",
        key="btn_run_alignment",
        disabled=not (seq1 and seq2),
    )
    clear_clicked = act2.button(
        "Clear Alignment",
        width="stretch",
        key="btn_clear_alignment",
    )

    if clear_clicked:
        for k in (
            "last_alignment",
            "alignment_active_in_session",
            "alignment_structure_fp",
            "alignment_inputs_snapshot",
        ):
            st.session_state.pop(k, None)
        st.rerun()

    if run_clicked:
        if not seq1 or not seq2:
            st.info("Provide both Sequence 1 and Sequence 2 to run pairwise alignment.")
            return None
        try:
            with st.spinner("Running pairwise sequence alignment..."):
                result = run_pairwise_alignment(
                    seq1=seq1,
                    seq2=seq2,
                    mode=mode_choice,  # type: ignore[arg-type]
                    matrix_name=matrix_name,
                    open_gap_score=float(open_gap),
                    extend_gap_score=float(extend_gap),
                    seq1_name=name1,
                    seq1_accession=acc1,
                    seq2_name=name2,
                    seq2_accession=acc2,
                )
            st.session_state["last_alignment"] = result
            st.session_state["alignment_active_in_session"] = True
            st.session_state["alignment_structure_fp"] = current_struct_fp
            st.session_state["alignment_inputs_snapshot"] = {
                "seq1": seq1,
                "seq2": seq2,
                "seq1_name": name1,
                "seq1_accession": acc1,
                "seq2_name": name2,
                "seq2_accession": acc2,
                "mode": mode_choice,
                "matrix_name": matrix_name,
                "open_gap": float(open_gap),
                "extend_gap": float(extend_gap),
            }
        except ValueError as exc:
            st.warning(str(exc))
            return None

    is_current_active = (
        st.session_state.get("alignment_active_in_session", False)
        and st.session_state.get("alignment_structure_fp") == current_struct_fp
    )
    result = st.session_state.get("last_alignment") if is_current_active else None

    if result is None:
        st.info(
            "Click **Run Pairwise Alignment** to execute an alignment for the current workspace session "
            "(pairwise results never leak automatically into PDF reports across structure switches)."
        )
        return None

    k1, k2, k3, k4, k5 = st.columns(5)
    k1.metric("Alignment score", f"{result.score:.1f}")
    k2.metric("Identity", f"{result.identity_pct:.1f}%", delta=f"{result.identical_count}/{result.alignment_length} cols")
    k3.metric("Similarity", f"{result.similarity_pct:.1f}%", delta=f"{result.similar_count}/{result.alignment_length} cols")
    k4.metric("Gaps", f"{result.gap_pct:.1f}%", delta=f"{result.gap_count} gap cols")
    k5.metric("Ungapped identity", f"{result.identity_ungapped_pct:.1f}%")

    st.info(result.interpretation)
    acc1_txt = f" [{result.seq1_accession}]" if result.seq1_accession else ""
    acc2_txt = f" [{result.seq2_accession}]" if result.seq2_accession else ""
    st.caption(
        f"**Seq 1:** {result.seq1_name}{acc1_txt} ({result.seq1_length} aa, span {result.seq1_span[0]}–{result.seq1_span[1]}) · "
        f"**Seq 2:** {result.seq2_name}{acc2_txt} ({result.seq2_length} aa, span {result.seq2_span[0]}–{result.seq2_span[1]}) · "
        f"Engine: {result.engine} · Symbols: `|` identical, `:` conservative (score > 0), `.` mismatch."
    )

    st.html(format_alignment_html(result, width=64))

    st.download_button(
        "Download alignment (.txt)",
        data=format_alignment_blocks(result, width=60),
        file_name=f"pairwise_alignment_{result.mode}_{result.matrix_name}.txt",
        mime="text/plain",
    )
    return result
