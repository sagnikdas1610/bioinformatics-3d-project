"""Primary Sequence Analytics tab."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from bio3d.charts import plot_aa_composition, plot_charge_curve, plot_hydropathy_profile
from bio3d.sequence import (
    SequenceProfile,
    charge_vs_ph_curve,
    composition_rows,
    kyte_doolittle_profile,
    validate_sequence,
)
from bio3d.structure import ParsedStructure


def render_sequence_tab(
    sequence: str,
    sequence_label: str,
    parsed: ParsedStructure | None,
    cached_analyze_fn,
) -> SequenceProfile | None:
    validation = validate_sequence(sequence)
    if not validation.is_valid:
        st.warning(validation.message)
        return None

    if validation.stripped_stop_codon:
        st.info(validation.message)

    # Warn cleanly if active sequence does not match any chain in the currently loaded 3D structure
    if parsed is not None:
        matched_chains = [
            ch for ch, ch_seq in parsed.chain_sequences.items() if ch_seq == validation.cleaned_sequence
        ]
        if matched_chains:
            st.caption(
                f"Sequence `{sequence_label}` matches loaded structure chain(s): **{', '.join(matched_chains)}**."
            )
        else:
            st.info(
                f"Active sequence (`{sequence_label}`, {validation.length} aa) differs from the chains in "
                f"`{parsed.filename}`. Use **Sync sequence from structure** in the sidebar to analyze a PDB chain directly."
            )

    profile: SequenceProfile = cached_analyze_fn(validation.cleaned_sequence)

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Sequence length", f"{profile.length:,} aa")
    c2.metric("Molecular weight", f"{profile.molecular_weight_da:,.1f} Da", delta=f"{profile.molecular_weight_kda:.2f} kDa")
    c3.metric("Theoretical pI", f"{profile.theoretical_pi:.2f}")
    c4.metric("Net charge (pH 7.0)", f"{profile.net_charge_ph7:+.2f}")

    c5, c6, c7, c8 = st.columns(4)
    c5.metric("GRAVY hydropathy", f"{profile.gravy:.3f}")
    c6.metric("Instability index", f"{profile.instability_index:.2f}", delta=profile.instability_class)
    c7.metric("Aromaticity", f"{profile.aromaticity:.3f}")
    c8.metric("ε280 (red / ox)", f"{profile.extinction_reduced:,} / {profile.extinction_oxidized:,}")

    st.plotly_chart(plot_aa_composition(profile), width="stretch")

    col_left, col_right = st.columns(2, gap="medium")
    with col_left:
        max_win = min(21, max(3, profile.length))
        default_win = min(9, max_win)
        window = st.slider("Kyte-Doolittle sliding window", min_value=3, max_value=max_win, value=default_win, step=2)
        pos, scores = kyte_doolittle_profile(profile.sequence, window=window)
        st.plotly_chart(plot_hydropathy_profile(pos, scores, window), width="stretch")
    with col_right:
        phs, charges = charge_vs_ph_curve(profile.sequence)
        st.plotly_chart(plot_charge_curve(phs, charges, profile.theoretical_pi), width="stretch")

    with st.expander("Amino-acid composition table & biological interpretation", expanded=False):
        df = pd.DataFrame(composition_rows(profile))
        st.dataframe(df, width="stretch", hide_index=True)
        st.markdown(
            "- **Molecular weight & ε280:** Calculated from average isotopic residue masses assuming unmodified N/C termini; molar extinction coefficients at 280 nm reflect reduced cysteines vs. cystine disulfide bridges.\n"
            "- **GRAVY (Grand Average of Hydropathicity):** Positive values indicate net hydrophobic character (common in transmembrane segments or buried cores); negative values indicate net hydrophilic surface tendency.\n"
            "- **Instability index:** Guruprasad et al. (1990) dipeptide heuristic (< 40 is classified as stable in vitro); it is a statistical sequence heuristic and does not guarantee cellular half-life or folding stability."
        )

    return profile
