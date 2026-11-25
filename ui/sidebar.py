"""Sidebar workspace input controls for structure and sequence loading."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from bio3d.config import (
    BUNDLED_STRUCTURES,
    MAX_FASTA_BYTES,
    MAX_STRUCTURE_BYTES,
    UBIQUITIN_SEQUENCE,
)
from bio3d.prediction import ExternalServiceError, fetch_rcsb_structure
from bio3d.sequence import parse_fasta_records
from bio3d.structure import ParsedStructure


def render_sidebar(active_parsed: ParsedStructure | None) -> tuple[str, str]:
    """Render sidebar inputs and return (active_sequence_text, active_sequence_label)."""
    with st.sidebar:
        st.subheader("Workspace Inputs")

        # ── 1. Structure source ────────────────────────────────────────────────
        st.markdown("**1. 3D Structure Source**")
        source_mode = st.segmented_control(
            "Structure source",
            options=["Samples", "Upload", "RCSB ID"],
            default="Samples",
            label_visibility="collapsed",
        )
        if source_mode is None:
            source_mode = "Samples"

        if source_mode == "Samples":
            sample_keys = list(BUNDLED_STRUCTURES.keys())
            selected_sample = st.selectbox("Bundled reference structure", sample_keys, index=0)
            info = BUNDLED_STRUCTURES[selected_sample]
            if st.button("Load sample structure", type="primary", width="stretch", key="btn_load_sample"):
                path = Path(info["path"])
                if path.exists():
                    st.session_state["structure_raw"] = path.read_text(encoding="utf-8", errors="replace")
                    st.session_state["structure_filename"] = path.name
                    st.session_state["structure_origin_hint"] = "experimental"
                    st.session_state["structure_error"] = None
                    st.rerun()
                else:
                    st.session_state["structure_error"] = f"Bundled sample not found: {path.name}"

        elif source_mode == "Upload":
            uploaded = st.file_uploader(
                "Upload PDB or mmCIF",
                type=["pdb", "ent", "cif", "mmcif"],
                accept_multiple_files=False,
            )
            if uploaded is not None:
                sig = f"{uploaded.name}:{uploaded.size}"
                if st.session_state.get("last_uploaded_sig") != sig:
                    if uploaded.size > MAX_STRUCTURE_BYTES:
                        st.session_state["structure_error"] = (
                            "Uploaded structure exceeds the 20 MiB safety limit."
                        )
                    else:
                        text = uploaded.getvalue().decode("utf-8", errors="replace")
                        st.session_state["structure_raw"] = text
                        st.session_state["structure_filename"] = uploaded.name
                        st.session_state["structure_origin_hint"] = None
                        st.session_state["structure_error"] = None
                        st.session_state["last_uploaded_sig"] = sig
                        st.rerun()

        else:  # RCSB ID
            c1, c2 = st.columns([2, 1], vertical_alignment="bottom")
            pdb_id_in = c1.text_input("RCSB PDB ID", value="1UBQ", max_chars=4)
            rcsb_fmt = c2.selectbox("Format", ["pdb", "cif"], index=0)
            if st.button("Fetch from RCSB", type="primary", width="stretch", key="btn_fetch_rcsb"):
                with st.spinner(f"Fetching {pdb_id_in.upper()} from RCSB PDB..."):
                    try:
                        text, fname = fetch_rcsb_structure(pdb_id_in, fmt=rcsb_fmt)  # type: ignore[arg-type]
                        st.session_state["structure_raw"] = text
                        st.session_state["structure_filename"] = fname
                        st.session_state["structure_origin_hint"] = "experimental"
                        st.session_state["structure_error"] = None
                        st.rerun()
                    except (ValueError, ExternalServiceError) as exc:
                        st.session_state["structure_error"] = str(exc)

        if st.session_state.get("structure_error"):
            st.error(st.session_state["structure_error"])

        st.divider()

        # ── 2. Protein sequence source ─────────────────────────────────────────
        st.markdown("**2. Primary Sequence Input**")
        if active_parsed and active_parsed.chains:
            sync_chain = st.selectbox(
                "Extract sequence from loaded chain",
                options=active_parsed.chains,
                index=0,
                key="sidebar_sync_chain",
            )
            if st.button("Sync sequence from structure", width="stretch", key="btn_sync_seq"):
                chain_seq = active_parsed.chain_sequences.get(sync_chain, "")
                if chain_seq:
                    header = f">{active_parsed.pdb_id or active_parsed.filename}_chain_{sync_chain}"
                    st.session_state["sequence_text_area"] = f"{header}\n{chain_seq}"
                    st.rerun()

        fasta_file = st.file_uploader(
            "Upload FASTA (.fasta / .fa / .txt)",
            type=["fasta", "fa", "faa", "txt"],
            accept_multiple_files=False,
        )
        if fasta_file is not None:
            fsig = f"{fasta_file.name}:{fasta_file.size}"
            if st.session_state.get("last_fasta_sig") != fsig:
                if fasta_file.size > MAX_FASTA_BYTES:
                    st.error("FASTA file exceeds the 2 MiB limit.")
                else:
                    st.session_state["sequence_text_area"] = fasta_file.getvalue().decode(
                        "utf-8", errors="replace"
                    )
                    st.session_state["last_fasta_sig"] = fsig
                    st.rerun()

        if "sequence_text_area" not in st.session_state:
            st.session_state["sequence_text_area"] = f">1UBQ_Ubiquitin\n{UBIQUITIN_SEQUENCE}"

        seq_raw = st.text_area(
            "FASTA or raw amino-acid sequence",
            key="sequence_text_area",
            height=155,
            help="Supports raw 1-letter codes or single/multi-FASTA records.",
        )

        records = parse_fasta_records(seq_raw)
        selected_seq = ""
        selected_label = "Sequence"
        if len(records) > 1:
            options = [f"{idx + 1}. {r.header} ({len(r.sequence)} aa)" for idx, r in enumerate(records)]
            picked = st.selectbox("Multi-FASTA record selector", options, index=0)
            rec_idx = options.index(picked)
            selected_seq = records[rec_idx].sequence
            selected_label = records[rec_idx].header
            st.caption(f"Parsed {len(records)} FASTA entries separately (not concatenated).")
        elif len(records) == 1:
            selected_seq = records[0].sequence
            selected_label = records[0].header

        st.divider()
        st.caption(
            "External databases: [RCSB PDB](https://www.rcsb.org/) · "
            "[UniProt](https://www.uniprot.org/) · "
            "[AlphaFold DB](https://alphafold.ebi.ac.uk/) · "
            "[NCBI BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi)"
        )

    return selected_seq, selected_label
