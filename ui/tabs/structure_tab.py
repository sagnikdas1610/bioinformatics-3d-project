"""Advanced Structural Analysis tab: Quality, Ramachandran, Secondary Structure, and Contact Map."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from bio3d.charts import (
    plot_bfactor_profile,
    plot_contact_map,
    plot_plddt_profile,
    plot_ramachandran,
    plot_secondary_structure,
)
from bio3d.config import ORIGIN_LABELS
from bio3d.structural_analysis import RAMA_CLASSES
from bio3d.structure import ParsedStructure


def render_structure_tab(
    parsed: ParsedStructure | None,
    cached_rama_fn,
    cached_ss_fn,
    cached_cmap_fn,
) -> None:
    if parsed is None:
        st.warning("Load a structure in the sidebar to run structural analytics.")
        return

    top_c1, top_c2 = st.columns([1.3, 2.7], vertical_alignment="bottom")
    chain_options = ["All"] + parsed.chains
    selected_chain = top_c1.selectbox("Chain scope (Ramachandran & Secondary Structure)", chain_options, index=0)
    sub_view = top_c2.segmented_control(
        "Analysis module",
        options=["Structure Quality", "Ramachandran Plot", "Secondary Structure", "Contact Map"],
        default="Structure Quality",
    ) or "Structure Quality"

    st.divider()

    if sub_view == "Structure Quality":
        q1, q2, q3, q4 = st.columns(4)
        q1.metric("Polymer chains", f"{len(parsed.chains)}")
        q2.metric("Polymer residues", f"{parsed.polymer_residue_count:,}")
        q3.metric("Polymer atoms (Model 1)", f"{parsed.polymer_atom_count:,}")
        q4.metric(
            "Resolution",
            f"{parsed.resolution_angstrom:.2f} Å" if parsed.resolution_angstrom is not None else "Unavailable",
        )

        st.markdown(f"**Origin:** {ORIGIN_LABELS.get(parsed.structure_origin, parsed.structure_origin)}")
        st.markdown(f"**Experimental method:** {parsed.experimental_method or 'Not reported in file header'}")
        if parsed.title:
            st.markdown(f"**Title:** {parsed.title}")
        st.markdown(
            f"**Non-polymer records:** {parsed.hetero_ligand_count} hetero ligand residues "
            f"({', '.join(parsed.hetero_ligand_names) if parsed.hetero_ligand_names else 'none'}) · "
            f"{parsed.water_count} water molecules · {parsed.model_count} model(s) in file."
        )

        st.markdown("##### Missing residues & chain continuity")
        st.info(parsed.missing_residue_note)
        if parsed.missing_residues:
            st.dataframe(pd.DataFrame(parsed.missing_residues), width="stretch", hide_index=True)

        # Strictly separate predicted pLDDT vs. experimental B-factor
        if parsed.is_predicted and parsed.plddt_stats:
            ps = parsed.plddt_stats
            st.markdown("##### AI prediction confidence (pLDDT)")
            p1, p2, p3, p4 = st.columns(4)
            p1.metric("Mean pLDDT", f"{ps['mean']:.1f}")
            p2.metric("Very high (≥ 90)", f"{ps['pct_very_high']:.1f}%")
            p3.metric("Confident (70–90)", f"{ps['pct_confident']:.1f}%")
            p4.metric("Low (< 70)", f"{ps['pct_low'] + ps['pct_very_low']:.1f}%")
            st.caption(f"Normalized to 0–100 pLDDT scale from raw B-factor column ({ps['raw_scale']}).")
            st.plotly_chart(plot_plddt_profile(ps), width="stretch")
        elif parsed.bfactor_stats:
            bs = parsed.bfactor_stats
            st.markdown("##### Crystallographic / experimental B-factor distribution (Å²)")
            b1, b2, b3, b4 = st.columns(4)
            b1.metric("Mean B-factor", f"{bs['mean']:.2f} Å²")
            b2.metric("Median B-factor", f"{bs['median']:.2f} Å²")
            b3.metric("Min", f"{bs['min']:.2f} Å²")
            b4.metric("Max", f"{bs['max']:.2f} Å²")
            st.caption(
                "Experimental atomic displacement parameters (B-factors in Å²) reflect local thermal motion or static disorder; "
                "they are distinct from AlphaFold/ESMFold pLDDT confidence scores."
            )
            st.plotly_chart(plot_bfactor_profile(bs), width="stretch")

    elif sub_view == "Ramachandran Plot":
        with st.spinner("Computing backbone φ/ψ dihedral angles..."):
            rama = cached_rama_fn(parsed.raw_text, parsed.filename, selected_chain)

        r1, r2, r3, r4 = st.columns(4)
        r1.metric("Evaluated residues", f"{rama.total_evaluated}", delta=f"{rama.excluded_count} termini/breaks excluded")
        r2.metric("Favored", f"{rama.favored_pct:.1f}%", delta=f"{rama.favored_count} residues")
        r3.metric("Allowed", f"{rama.allowed_pct:.1f}%", delta=f"{rama.allowed_count} residues")
        r4.metric("Outliers", f"{rama.outlier_pct:.1f}%", delta=f"{rama.outlier_count} residues")

        st.info(rama.interpretation)
        st.caption(f"Reference dataset: **{rama.reference_source}** (Favored ≥ 0.02, Allowed ≥ 0.0005).")
        bg_class = st.selectbox(
            "Reference contour background (MolProbity Top8000)",
            list(RAMA_CLASSES),
            index=0,
        )
        st.plotly_chart(plot_ramachandran(rama, bg_class=bg_class), width="stretch")

        if rama.outliers:
            st.markdown("##### Residues in outlier regions")
            out_df = pd.DataFrame(
                [
                    {
                        "Chain": o.chain,
                        "Residue": f"{o.resname} {o.resseq}{o.icode}",
                        "Rama Class": o.rama_class,
                        "φ (deg)": round(o.phi, 2),
                        "ψ (deg)": round(o.psi, 2),
                    }
                    for o in rama.outliers
                ]
            )
            st.dataframe(out_df, width="stretch", hide_index=True)
        else:
            st.success("Zero Ramachandran outliers detected under Lovell et al. (2003) 4-class criteria.")

        with st.expander("Biological relevance of Ramachandran analysis", expanded=False):
            st.markdown(
                "Steric clashes between backbone carbonyl/amide atoms and side-chain Cβ atoms restrict protein φ/ψ torsion angles to specific favored basins "
                "(right-handed α-helix, β-sheet/PPII, and left-handed α_L). Glycine lacks a side chain and samples symmetric conformations, while Proline's pyrrolidine ring "
                "locks φ near −65° and restricts the preceding (Pre-Proline) residue. High favored percentages (> 95%) and minimal outliers (< 1–2%) indicate clean stereochemical geometry."
            )

    elif sub_view == "Secondary Structure":
        with st.spinner("Computing Kabsch-Sander backbone hydrogen-bond secondary structure..."):
            ss = cached_ss_fn(parsed.raw_text, parsed.filename, selected_chain)

        s1, s2, s3, s4 = st.columns(4)
        s1.metric("Residues analyzed", f"{ss.total_residues}")
        s2.metric("Helix (H / G / I)", f"{ss.helix_pct:.1f}%", delta=f"{ss.helix_count} residues")
        s3.metric("Sheet (E / B)", f"{ss.sheet_pct:.1f}%", delta=f"{ss.sheet_count} residues")
        s4.metric("Coil & Turn (T / C)", f"{ss.coil_pct:.1f}%", delta=f"{ss.coil_count} residues")

        st.caption(f"Method: **{ss.method}**")
        if ss.header_available and ss.header_helix_pct is not None:
            st.info(
                f"File header `HELIX`/`SHEET` records comparison — "
                f"Header Helix: **{ss.header_helix_pct:.1f}%** · "
                f"Header Sheet: **{ss.header_sheet_pct:.1f}%** · "
                f"Header Coil/Turn: **{ss.header_coil_pct:.1f}%**."
            )

        if ss.residues:
            st.plotly_chart(plot_secondary_structure(ss), width="stretch")
            with st.expander("6-state breakdown & per-residue assignments", expanded=False):
                st.write("6-state residue counts:", ss.state6_counts)
                ss_df = pd.DataFrame(
                    [
                        {
                            "Chain": r.chain,
                            "ResSeq": r.resseq,
                            "Residue": r.resname,
                            "DSSP Code": r.dssp_code,
                            "3-State": r.state3,
                            "Header Annotation": r.header_state3 or "—",
                        }
                        for r in ss.residues
                    ]
                )
                st.dataframe(ss_df, width="stretch", hide_index=True)

    else:  # Contact Map
        c_col1, c_col2, c_col3, c_col4 = st.columns(4)
        cm_chain = c_col1.selectbox("Chain", parsed.chains, index=0, key="cm_chain_select")
        atom_mode = c_col2.segmented_control("Atom anchor", ["CA", "CB"], default="CA") or "CA"
        cutoff = c_col3.slider("Contact threshold (Å)", min_value=4.0, max_value=15.0, value=8.0, step=0.5)
        min_sep = c_col4.slider("Min seq separation |i-j|", min_value=1, max_value=24, value=6, step=1)

        with st.spinner("Calculating residue-residue distance matrix..."):
            cmap = cached_cmap_fn(parsed.raw_text, parsed.filename, cm_chain, atom_mode, float(cutoff), int(min_sep))

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Residues mapped", f"{cmap.residue_count}")
        m2.metric("Total contacts", f"{cmap.total_contacts:,}")
        m3.metric("Long-range (|i-j| ≥ 12)", f"{cmap.long_range_contacts:,}")
        m4.metric("Contact density", f"{cmap.contact_density_pct:.2f}%")

        if cmap.atom_mode == "CB":
            st.caption(
                f"Distance threshold: **≤ {cmap.threshold_angstrom:.1f} Å** (`CB` mode; "
                f"Glycine lacks Cβ, so **Cα fallback** was applied cleanly to {cmap.glycine_fallback_count} Glycine residue(s)"
                + (f" and {cmap.missing_cb_fallback_count} residue(s) missing CB" if cmap.missing_cb_fallback_count else "")
                + ")."
            )
        else:
            st.caption(f"Distance threshold: **≤ {cmap.threshold_angstrom:.1f} Å** between backbone **Cα** atoms (`|i - j| ≥ {cmap.min_seq_sep}`).")

        if cmap.truncated_from:
            st.warning(f"Chain {cmap.chain_id} has {cmap.truncated_from} residues; capped at the first {cmap.residue_count} residues for memory safety.")

        st.plotly_chart(plot_contact_map(cmap), width="stretch")

        if cmap.top_contacts:
            with st.expander("Explore residue-residue contact pairs (sorted by sequence separation)", expanded=False):
                st.dataframe(pd.DataFrame(cmap.top_contacts), width="stretch", hide_index=True)
