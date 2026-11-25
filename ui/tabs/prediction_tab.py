"""AI-Assisted Protein Structure Prediction & AlphaFold DB Lookup tab."""

from __future__ import annotations

import streamlit as st

from bio3d.charts import plot_plddt_profile
from bio3d.config import MAX_ESMFOLD_RESIDUES, MIN_ESMFOLD_RESIDUES
from bio3d.prediction import (
    ExternalServiceError,
    fetch_alphafold_structure,
    predict_with_esmfold,
)
from bio3d.sequence import validate_sequence
from bio3d.structure import parse_structure


def _load_into_workspace(pdb_text: str, filename: str, origin_hint: str) -> None:
    st.session_state["structure_raw"] = pdb_text
    st.session_state["structure_filename"] = filename
    st.session_state["structure_origin_hint"] = origin_hint
    st.session_state["structure_error"] = None
    st.session_state["main_nav_section"] = "3D Explorer"


def render_prediction_tab(active_sequence: str) -> None:
    st.warning(
        "Computational predictions (AlphaFold DB and ESMFold) are structural hypotheses rather than experimentally determined coordinates. "
        "Always inspect residue-level pLDDT confidence scores and verify stereochemistry before drawing biological conclusions.",
        icon=":material/science:",
    )

    col_af, col_esm = st.columns(2, gap="large")

    # ── 1. AlphaFold Protein Structure Database ───────────────────────────────
    with col_af:
        with st.container(border=True):
            st.markdown("#### AlphaFold DB Lookup (UniProt Accession)")
            st.caption(
                "Retrieves precomputed AlphaFold structures directly from EMBL-EBI AlphaFold DB without server-side GPU compute."
            )
            uniprot_in = st.text_input(
                "UniProt accession",
                value="P0CG48",
                help="Examples: P0CG48 (Polyubiquitin), P69905 (Hemoglobin subunit alpha), P04637 (TP53).",
            )
            if st.button("Fetch AlphaFold model", type="primary", width="stretch", key="btn_af_fetch"):
                with st.spinner(f"Querying AlphaFold DB for {uniprot_in.upper()}..."):
                    try:
                        af_res = fetch_alphafold_structure(uniprot_in)
                        af_parsed = parse_structure(
                            af_res["pdb_text"],
                            filename=af_res["filename"],
                            origin_hint="predicted_alphafold",
                        )
                        st.session_state["af_result"] = af_res
                        st.session_state["af_parsed"] = af_parsed
                    except (ValueError, ExternalServiceError) as exc:
                        st.error(str(exc))

            af_res = st.session_state.get("af_result")
            af_parsed = st.session_state.get("af_parsed")
            if af_res and af_parsed:
                st.success(
                    f"Loaded **{af_res['entry_id']}** ({af_res['gene']} · {af_res['organism']})"
                )
                if af_parsed.plddt_stats:
                    ps = af_parsed.plddt_stats
                    m1, m2 = st.columns(2)
                    m1.metric("Mean pLDDT", f"{ps['mean']:.1f}")
                    m2.metric("Confident (≥ 70)", f"{ps['pct_very_high'] + ps['pct_confident']:.1f}%")
                    st.plotly_chart(plot_plddt_profile(ps), width="stretch")

                b1, b2 = st.columns(2)
                b1.button(
                    "Open in 3D Explorer",
                    width="stretch",
                    key="btn_af_open",
                    on_click=_load_into_workspace,
                    args=(af_res["pdb_text"], af_res["filename"], "predicted_alphafold"),
                )
                b2.download_button(
                    "Download PDB",
                    data=af_res["pdb_text"],
                    file_name=af_res["filename"],
                    mime="chemical/x-pdb",
                    width="stretch",
                )

    # ── 2. Live ESMFold API Folding ───────────────────────────────────────────
    with col_esm:
        with st.container(border=True):
            st.markdown("#### Live Single-Sequence Prediction (ESMFold API)")
            st.caption(
                f"Submits a single chain ({MIN_ESMFOLD_RESIDUES}–{MAX_ESMFOLD_RESIDUES} standard amino acids) to the public ESM Atlas API. "
                "Typical runtime is 10–45 seconds; public endpoints may rate-limit."
            )
            if "esm_seq_input" not in st.session_state:
                st.session_state["esm_seq_input"] = active_sequence or ""

            if active_sequence and st.button("Copy active sidebar sequence", width="stretch", key="btn_copy_esm"):
                st.session_state["esm_seq_input"] = active_sequence
                st.rerun()

            esm_seq_raw = st.text_area("Input sequence (10–400 standard aa)", key="esm_seq_input", height=115)
            val = validate_sequence(esm_seq_raw)
            can_run = val.is_valid and MIN_ESMFOLD_RESIDUES <= val.length <= MAX_ESMFOLD_RESIDUES
            if esm_seq_raw.strip() and not val.is_valid:
                st.warning(val.message)
            elif val.is_valid and not can_run:
                st.warning(
                    f"Length is {val.length} aa. Live ESMFold API supports {MIN_ESMFOLD_RESIDUES}–{MAX_ESMFOLD_RESIDUES} residues; "
                    "use ColabFold below for longer sequences."
                )

            if st.button("Run ESMFold prediction", type="primary", width="stretch", disabled=not can_run, key="btn_esm_run"):
                with st.spinner("Folding sequence via ESMFold Atlas API (10–45s)..."):
                    try:
                        esm_res = predict_with_esmfold(val.cleaned_sequence)
                        esm_parsed = parse_structure(
                            str(esm_res["pdb_text"]),
                            filename=str(esm_res["filename"]),
                            origin_hint="predicted_esmfold",
                        )
                        st.session_state["esm_result"] = esm_res
                        st.session_state["esm_parsed"] = esm_parsed
                    except (ValueError, ExternalServiceError) as exc:
                        st.error(str(exc))

            esm_res = st.session_state.get("esm_result")
            esm_parsed = st.session_state.get("esm_parsed")
            if esm_res and esm_parsed:
                st.success(f"Prediction complete ({esm_res['length']} residues via {esm_res['service']}).")
                if esm_parsed.plddt_stats:
                    ps = esm_parsed.plddt_stats
                    m1, m2 = st.columns(2)
                    m1.metric("Mean pLDDT", f"{ps['mean']:.1f}")
                    m2.metric("Confident (≥ 70)", f"{ps['pct_very_high'] + ps['pct_confident']:.1f}%")
                    st.plotly_chart(plot_plddt_profile(ps), width="stretch")

                b1, b2 = st.columns(2)
                b1.button(
                    "Open in 3D Explorer",
                    width="stretch",
                    key="btn_esm_open",
                    on_click=_load_into_workspace,
                    args=(str(esm_res["pdb_text"]), str(esm_res["filename"]), "predicted_esmfold"),
                )
                b2.download_button(
                    "Download PDB",
                    data=str(esm_res["pdb_text"]),
                    file_name=str(esm_res["filename"]),
                    mime="chemical/x-pdb",
                    width="stretch",
                )

    st.divider()
    with st.expander("pLDDT confidence guide & ColabFold workflow for >400 aa / multi-chain complexes", expanded=False):
        st.markdown(
            "- **pLDDT ≥ 90 (Dark Blue):** Very high backbone and side-chain confidence.\n"
            "- **70 ≤ pLDDT < 90 (Light Blue):** Confident backbone fold.\n"
            "- **50 ≤ pLDDT < 70 (Yellow):** Low confidence; interpret local geometry cautiously.\n"
            "- **pLDDT < 50 (Orange):** Very low confidence; often corresponds to intrinsically disordered regions.\n\n"
            "For proteins longer than 400 residues or multi-chain complexes: "
            "[ColabFold Notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) · "
            "[Repository Notebook (`AI_Structure_Generator.ipynb`)](https://colab.research.google.com/github/sagnikdas1610/bioinformatics-3d-project/blob/main/AI_Structure_Generator.ipynb)"
        )
