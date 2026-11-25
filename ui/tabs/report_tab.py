"""Scientific PDF & Data Export Center tab."""

from __future__ import annotations

import streamlit as st

from bio3d.alignment import AlignmentResult, format_alignment_blocks
from bio3d.reporting import generate_pdf_report, generate_text_report
from bio3d.sequence import SequenceProfile, kyte_doolittle_profile, validate_sequence
from bio3d.structure import ParsedStructure


def render_report_tab(
    parsed: ParsedStructure | None,
    sequence: str,
    sequence_label: str,
    cached_analyze_fn,
    cached_rama_fn,
    cached_ss_fn,
    cached_cmap_fn,
) -> None:
    st.markdown("Export publication-style PDF reports and reproducible analysis artifacts.")

    val = validate_sequence(sequence)
    seq_profile: SequenceProfile | None = cached_analyze_fn(val.cleaned_sequence) if val.is_valid else None
    hydro = kyte_doolittle_profile(val.cleaned_sequence, window=9) if val.is_valid else None
    alignment: AlignmentResult | None = st.session_state.get("last_alignment")

    rama = None
    ss = None
    cmap = None
    if parsed is not None:
        primary_chain = parsed.chains[0]
        rama = cached_rama_fn(parsed.raw_text, parsed.filename, "All")
        ss = cached_ss_fn(parsed.raw_text, parsed.filename, "All")
        cmap = cached_cmap_fn(parsed.raw_text, parsed.filename, primary_chain, "CA", 8.0, 6)

    # Checklist showing exactly what is available vs unavailable
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("##### Report contents status")
        st.write("Structure & Quality:", f"Included (`{parsed.filename}`)" if parsed else "Unavailable")
        st.write("Sequence Profile:", f"Included (`{sequence_label}`, {seq_profile.length} aa)" if seq_profile else "Unavailable")
        st.write("Pairwise Alignment:", f"Included ({alignment.mode}, {alignment.matrix_name})" if alignment else "Not run in session (marked Unavailable)")
    with c2:
        st.markdown("##### Structural analytics included")
        st.write("Ramachandran Plot:", f"Included ({rama.favored_pct:.1f}% favored)" if rama else "Unavailable")
        st.write("Secondary Structure:", f"Included ({ss.helix_pct:.1f}% H / {ss.sheet_pct:.1f}% E)" if ss else "Unavailable")
        st.write("Contact Map:", f"Included (Chain {cmap.chain_id}, {cmap.total_contacts} contacts)" if cmap else "Unavailable")

    st.divider()

    report_title = st.text_input(
        "Custom PDF report title (optional)",
        value="",
        placeholder="BioFold Workbench — Structural & Sequence Analysis Report",
    )

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
            )
            st.session_state["generated_pdf_bytes"] = pdf_bytes

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

    txt_report = generate_text_report(parsed, seq_profile, alignment, rama, ss, cmap)
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
        fasta_out = f">{sequence_label} length={seq_profile.length}\n{wrapped}\n"
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
