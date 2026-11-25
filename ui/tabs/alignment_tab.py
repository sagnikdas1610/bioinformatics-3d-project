"""Pairwise Sequence Alignment tab."""

from __future__ import annotations

import streamlit as st

from bio3d.alignment import (
    SUPPORTED_MATRICES,
    AlignmentResult,
    format_alignment_blocks,
    format_alignment_html,
    run_pairwise_alignment,
)
from bio3d.config import HBA_HUMAN_SEQ, HBB_HUMAN_SEQ
from bio3d.sequence import parse_fasta_records


def _extract_first_seq(raw_text: str) -> str:
    records = parse_fasta_records(raw_text)
    return records[0].sequence if records else ""


def render_alignment_tab(active_sequence: str) -> AlignmentResult | None:
    st.markdown("Compare two amino-acid sequences using global (Needleman–Wunsch) or local (Smith–Waterman) alignment.")

    if "aln_seq1" not in st.session_state:
        st.session_state["aln_seq1"] = HBA_HUMAN_SEQ
    if "aln_seq2" not in st.session_state:
        st.session_state["aln_seq2"] = HBB_HUMAN_SEQ

    preset_col1, preset_col2, _ = st.columns([1.2, 1.4, 2.0])
    if preset_col1.button("Load HBA vs HBB demo", width="stretch"):
        st.session_state["aln_seq1"] = HBA_HUMAN_SEQ
        st.session_state["aln_seq2"] = HBB_HUMAN_SEQ
        st.rerun()
    if active_sequence and preset_col2.button("Use sidebar sequence as Seq 1", width="stretch"):
        st.session_state["aln_seq1"] = active_sequence
        st.rerun()

    col_s1, col_s2 = st.columns(2, gap="medium")
    with col_s1:
        raw_s1 = st.text_area("Sequence 1 (FASTA or raw)", key="aln_seq1", height=130)
    with col_s2:
        raw_s2 = st.text_area("Sequence 2 (FASTA or raw)", key="aln_seq2", height=130)

    p1, p2, p3, p4 = st.columns([1.1, 1.1, 1.0, 1.0])
    mode_choice = p1.segmented_control("Alignment mode", ["global", "local"], default="global") or "global"
    matrix_name = p2.selectbox("Substitution matrix", list(SUPPORTED_MATRICES), index=0)
    open_gap = p3.number_input("Gap open penalty", min_value=-30.0, max_value=0.0, value=-10.0, step=0.5)
    extend_gap = p4.number_input("Gap extend penalty", min_value=-10.0, max_value=0.0, value=-0.5, step=0.5)

    seq1 = _extract_first_seq(raw_s1)
    seq2 = _extract_first_seq(raw_s2)

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
            )
    except ValueError as exc:
        st.warning(str(exc))
        return None

    st.session_state["last_alignment"] = result

    k1, k2, k3, k4, k5 = st.columns(5)
    k1.metric("Alignment score", f"{result.score:.1f}")
    k2.metric("Identity", f"{result.identity_pct:.1f}%", delta=f"{result.identical_count}/{result.alignment_length} cols")
    k3.metric("Similarity", f"{result.similarity_pct:.1f}%", delta=f"{result.similar_count}/{result.alignment_length} cols")
    k4.metric("Gaps", f"{result.gap_pct:.1f}%", delta=f"{result.gap_count} gap cols")
    k5.metric("Ungapped identity", f"{result.identity_ungapped_pct:.1f}%")

    st.info(result.interpretation)
    st.caption(
        f"Aligned spans: Seq 1 [{result.seq1_span[0]}–{result.seq1_span[1]} / {result.seq1_length} aa] · "
        f"Seq 2 [{result.seq2_span[0]}–{result.seq2_span[1]} / {result.seq2_length} aa] · "
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
