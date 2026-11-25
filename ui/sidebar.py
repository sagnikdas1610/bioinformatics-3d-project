"""Sidebar workspace input controls for structure and structure-derived / custom sequence loading."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from bio3d.config import (
    BUNDLED_STRUCTURES,
    MAX_FASTA_BYTES,
    MAX_STRUCTURE_BYTES,
)
from bio3d.prediction import ExternalServiceError, fetch_rcsb_structure
from bio3d.sequence import parse_fasta_records, validate_sequence
from bio3d.structure import ParsedStructure
from ui.state import ActiveSequenceContext, set_active_structure


def _build_mismatch_warning(
    seq_label: str,
    cleaned_seq: str,
    active_parsed: ParsedStructure | None,
) -> tuple[bool | None, list[str], str | None]:
    if active_parsed is None or not cleaned_seq:
        return None, [], None
    matched = [ch for ch, ch_seq in active_parsed.chain_sequences.items() if ch_seq == cleaned_seq]
    if matched:
        return True, matched, None
    chain_desc = ", ".join(f"Chain {ch} ({len(s)} aa)" for ch, s in active_parsed.chain_sequences.items())
    struct_name = active_parsed.pdb_id or active_parsed.filename
    warning = (
        f"Sequence mismatch: custom sequence `{seq_label}` ({len(cleaned_seq)} aa) does not match "
        f"loaded structure `{struct_name}` ({chain_desc}). Analytics below reflect the custom FASTA sequence, "
        f"whereas 3D coordinates remain `{struct_name}`."
    )
    return False, [], warning


def render_sidebar(active_parsed: ParsedStructure | None) -> ActiveSequenceContext:
    """Render sidebar inputs and return explicit ActiveSequenceContext."""
    with st.sidebar:
        st.subheader("Workspace Inputs")

        # ── 1. Structure source ────────────────────────────────────────────────
        st.markdown("**1. 3D Structure Source**")
        st.session_state.setdefault("sidebar_structure_source", "Samples")
        source_mode = st.segmented_control(
            "Structure source",
            options=["Samples", "Upload", "RCSB ID"],
            label_visibility="collapsed",
            key="sidebar_structure_source",
        )
        if source_mode is None:
            source_mode = "Samples"

        if source_mode == "Samples":
            sample_keys = list(BUNDLED_STRUCTURES.keys())
            st.session_state.setdefault("sidebar_sample_select", sample_keys[0])
            if st.session_state["sidebar_sample_select"] not in sample_keys:
                st.session_state["sidebar_sample_select"] = sample_keys[0]
            selected_sample = st.selectbox(
                "Bundled reference structure",
                sample_keys,
                key="sidebar_sample_select",
            )
            info = BUNDLED_STRUCTURES[selected_sample]
            if st.button("Load sample structure", type="primary", width="stretch", key="btn_load_sample"):
                path = Path(info["path"])
                if path.exists():
                    set_active_structure(
                        path.read_text(encoding="utf-8", errors="replace"),
                        path.name,
                        "experimental",
                        reset_sequence_mode=True,
                    )
                    st.rerun()
                else:
                    st.session_state["structure_error"] = f"Bundled sample not found: {path.name}"

        elif source_mode == "Upload":
            uploaded = st.file_uploader(
                "Upload PDB or mmCIF",
                type=["pdb", "ent", "cif", "mmcif"],
                accept_multiple_files=False,
                key="sidebar_structure_upload",
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
                        set_active_structure(text, uploaded.name, None, reset_sequence_mode=True)
                        st.session_state["last_uploaded_sig"] = sig
                        st.rerun()

        else:  # RCSB ID
            c1, c2 = st.columns([2, 1], vertical_alignment="bottom")
            pdb_id_in = c1.text_input("RCSB PDB ID", value="1UBQ", max_chars=4, key="sidebar_rcsb_id")
            rcsb_fmt = c2.selectbox("Format", ["pdb", "cif"], index=0, key="sidebar_rcsb_fmt")
            if st.button("Fetch from RCSB", type="primary", width="stretch", key="btn_fetch_rcsb"):
                with st.spinner(f"Fetching {pdb_id_in.upper()} from RCSB PDB..."):
                    try:
                        text, fname = fetch_rcsb_structure(pdb_id_in, fmt=rcsb_fmt)  # type: ignore[arg-type]
                        set_active_structure(text, fname, "experimental", reset_sequence_mode=True)
                        st.rerun()
                    except (ValueError, ExternalServiceError) as exc:
                        st.session_state["structure_error"] = str(exc)

        if st.session_state.get("structure_error"):
            st.error(st.session_state["structure_error"])

        st.divider()

        # ── 2. Primary sequence source (structure-derived vs independent FASTA) ─
        st.markdown("**2. Primary Sequence Source**")
        st.session_state.setdefault("sequence_source_mode", "From Structure")
        seq_mode = st.segmented_control(
            "Sequence source mode",
            options=["From Structure", "Custom FASTA"],
            label_visibility="collapsed",
            key="sequence_source_mode",
        )
        if seq_mode is None:
            seq_mode = "From Structure"

        if seq_mode == "From Structure":
            if active_parsed and active_parsed.chains:
                if st.session_state.get("sidebar_structure_chain") not in active_parsed.chains:
                    st.session_state["sidebar_structure_chain"] = active_parsed.chains[0]
                selected_chain = st.selectbox(
                    "Active chain for sequence & report",
                    options=active_parsed.chains,
                    key="sidebar_structure_chain",
                )
                chain_seq = active_parsed.chain_sequences.get(selected_chain, "")
                struct_id = active_parsed.pdb_id or active_parsed.filename
                seq_label = f"{struct_id}_chain_{selected_chain}"
                extraction_note = active_parsed.chain_extraction_notes.get(selected_chain)

                st.caption(
                    f"Derived live from **{struct_id} · Chain {selected_chain}** ({len(chain_seq)} residues)."
                )
                if chain_seq:
                    st.code(chain_seq[:60] + ("..." if len(chain_seq) > 60 else ""), language=None)
                if active_parsed.unmapped_polymer_residues.get(selected_chain):
                    st.warning(extraction_note or "Unmapped residues detected in polymer chain.")
                elif active_parsed.modified_residues_mapped.get(selected_chain):
                    st.info(extraction_note)

                ctx = ActiveSequenceContext(
                    sequence=chain_seq,
                    label=seq_label,
                    source="structure",
                    chain_id=selected_chain,
                    matches_structure=True,
                    matched_chains=[selected_chain],
                    mismatch_warning=None,
                    extraction_note=extraction_note,
                )
            else:
                st.info("Load a 3D structure above to extract its chain sequence automatically.")
                ctx = ActiveSequenceContext(
                    sequence="",
                    label="Unavailable",
                    source="unavailable",
                    chain_id=None,
                    matches_structure=None,
                    matched_chains=[],
                    mismatch_warning=None,
                    extraction_note="No 3D structure is currently loaded.",
                )
        else:  # Custom FASTA
            st.session_state.setdefault("custom_fasta_text", "")
            if active_parsed and active_parsed.chains:
                default_ch = (
                    st.session_state.get("sidebar_structure_chain")
                    if st.session_state.get("sidebar_structure_chain") in active_parsed.chains
                    else active_parsed.chains[0]
                )
                if st.button("Copy active structure chain", width="stretch", key="btn_copy_struct_to_custom"):
                    cseq = active_parsed.chain_sequences.get(default_ch, "")
                    if cseq:
                        st.session_state["custom_fasta_text"] = (
                            f">{active_parsed.pdb_id or active_parsed.filename}_chain_{default_ch}\n{cseq}"
                        )
                        st.rerun()

            fasta_file = st.file_uploader(
                "Upload FASTA (.fasta / .fa / .txt)",
                type=["fasta", "fa", "faa", "txt"],
                accept_multiple_files=False,
                key="sidebar_fasta_upload",
            )
            if fasta_file is not None:
                fsig = f"{fasta_file.name}:{fasta_file.size}"
                if st.session_state.get("last_fasta_sig") != fsig:
                    if fasta_file.size > MAX_FASTA_BYTES:
                        st.error("FASTA file exceeds the 2 MiB limit.")
                    else:
                        st.session_state["custom_fasta_text"] = fasta_file.getvalue().decode(
                            "utf-8", errors="replace"
                        )
                        st.session_state["last_fasta_sig"] = fsig
                        st.rerun()

            seq_raw = st.text_area(
                "Independent FASTA or raw amino-acid sequence",
                key="custom_fasta_text",
                height=145,
                help="Analyzed independently of the 3D structure unless it exactly matches a chain.",
            )

            records = parse_fasta_records(seq_raw)
            selected_seq = ""
            selected_label = "Custom_sequence"
            if len(records) > 1:
                options = [f"{idx + 1}. {r.header} ({len(r.sequence)} aa)" for idx, r in enumerate(records)]
                picked = st.selectbox("Multi-FASTA record selector", options, index=0, key="sidebar_multifasta_select")
                rec_idx = options.index(picked)
                selected_seq = records[rec_idx].sequence
                selected_label = records[rec_idx].header
                st.caption(f"Parsed {len(records)} FASTA entries separately (not concatenated).")
            elif len(records) == 1:
                selected_seq = records[0].sequence
                selected_label = records[0].header

            val = validate_sequence(selected_seq)
            check_seq = val.cleaned_sequence if val.is_valid else selected_seq.strip()
            matches, matched_chains, mismatch_warning = _build_mismatch_warning(
                selected_label, check_seq, active_parsed
            )
            if mismatch_warning:
                st.warning(mismatch_warning)

            ctx = ActiveSequenceContext(
                sequence=selected_seq,
                label=selected_label,
                source="custom_fasta" if selected_seq.strip() else "unavailable",
                chain_id=matched_chains[0] if len(matched_chains) == 1 else None,
                matches_structure=matches,
                matched_chains=matched_chains,
                mismatch_warning=mismatch_warning,
                extraction_note=None if selected_seq.strip() else "No custom FASTA sequence was entered.",
            )

        st.divider()
        st.caption(
            "External databases: [RCSB PDB](https://www.rcsb.org/) · "
            "[UniProt](https://www.uniprot.org/) · "
            "[AlphaFold DB](https://alphafold.ebi.ac.uk/) · "
            "[NCBI BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi)"
        )

    return ctx
