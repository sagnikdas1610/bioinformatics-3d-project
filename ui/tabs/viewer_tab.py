"""Interactive 3D Molecular Explorer tab."""

from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components

from bio3d.config import ORIGIN_LABELS
from bio3d.structure import ParsedStructure
from bio3d.viewer import ColorMode, RenderStyle, ViewerConfig, build_viewer_html


def _parse_residue_selection(spec: str) -> list[int]:
    """Parse comma-separated numbers and ranges like '25, 40-44, 70' into residue IDs."""
    ids: list[int] = []
    for token in (spec or "").replace(";", ",").split(","):
        tok = token.strip()
        if not tok:
            continue
        if "-" in tok and not tok.startswith("-"):
            parts = tok.split("-", 1)
            if parts[0].strip().isdigit() and parts[1].strip().isdigit():
                a, b = int(parts[0].strip()), int(parts[1].strip())
                lo, hi = min(a, b), max(a, b)
                if hi - lo <= 150:
                    ids.extend(range(lo, hi + 1))
        elif tok.lstrip("-").isdigit():
            ids.append(int(tok))
    return sorted(set(ids))[:150]


def render_viewer_tab(parsed: ParsedStructure | None) -> None:
    if parsed is None:
        st.warning("Load a sample structure, upload a PDB/mmCIF file, or fetch an ID from RCSB PDB in the sidebar.")
        return

    viewer_col, ctrl_col = st.columns([2.35, 1.0], gap="large", vertical_alignment="top")

    with ctrl_col:
        st.markdown("#### Viewer Controls")
        render_style: RenderStyle = st.selectbox(
            "Representation",
            ["Cartoon", "Stick", "Ball & Stick", "Sphere"],
            index=0,
        )  # type: ignore[assignment]

        # Strictly gate color modes by structure origin so B-factors are never mislabeled as pLDDT
        if parsed.is_predicted:
            color_options: list[ColorMode] = [
                "Predicted pLDDT",
                "Secondary Structure",
                "By Chain",
                "N->C Spectrum",
            ]
        else:
            color_options = [
                "Secondary Structure",
                "By Chain",
                "N->C Spectrum",
                "Experimental B-factor",
            ]
        color_mode: ColorMode = st.selectbox("Color scheme", color_options, index=0)  # type: ignore[assignment]

        chain_choices = ["All"] + parsed.chains
        if st.session_state.get("viewer_focus_chain") not in chain_choices:
            st.session_state["viewer_focus_chain"] = "All"
        selected_chain = st.selectbox("Focus chain", chain_choices, key="viewer_focus_chain")

        c_surf, c_lig = st.columns(2)
        show_surface = c_surf.toggle("Surface", value=False)
        show_ligands = c_lig.toggle("Ligands", value=True)
        surface_opacity = 0.65
        if show_surface:
            if parsed.polymer_atom_count > 20000:
                st.caption("Surface disabled automatically for structures >20,000 atoms to preserve browser responsiveness.")
            else:
                surface_opacity = st.slider("Surface opacity", 0.15, 0.95, 0.65, 0.05)

        st.markdown("#### Residue Highlighting")
        st.session_state.setdefault("viewer_hi_text", "")
        hi_text = st.text_input(
            "Residue numbers (e.g. 25, 44-48)",
            placeholder="e.g. 16, 44-48",
            key="viewer_hi_text",
        )
        hi_ids = _parse_residue_selection(hi_text)
        c_lbl, c_zm = st.columns(2)
        show_labels = c_lbl.checkbox("Show labels", value=True, disabled=not hi_ids)
        zoom_sel = c_zm.checkbox("Zoom selection", value=False, disabled=not hi_ids)

        st.divider()
        st.markdown("#### Structure Metadata")
        m1, m2 = st.columns(2)
        m1.metric("Polymer residues", f"{parsed.polymer_residue_count:,}")
        m2.metric("Chains", f"{len(parsed.chains)}")
        m3, m4 = st.columns(2)
        m3.metric("Atoms (Model 1)", f"{parsed.polymer_atom_count:,}")
        m4.metric("Waters", f"{parsed.water_count:,}")

        st.caption(f"**Origin:** {ORIGIN_LABELS.get(parsed.structure_origin, parsed.structure_origin)}")
        if parsed.experimental_method:
            st.caption(f"**Method:** {parsed.experimental_method}")
        if parsed.resolution_angstrom is not None:
            st.caption(f"**Resolution:** {parsed.resolution_angstrom:.2f} Å")
        if parsed.hetero_ligand_names:
            st.caption(f"**Ligands:** {', '.join(parsed.hetero_ligand_names[:10])}")
        for w in parsed.warnings:
            st.info(w)

        st.download_button(
            f"Download {parsed.file_format.upper()} file",
            data=parsed.raw_text,
            file_name=parsed.filename,
            mime="chemical/x-pdb" if parsed.file_format == "pdb" else "chemical/x-cif",
            width="stretch",
        )

    b_min = float(parsed.bfactor_stats["p5"]) if parsed.bfactor_stats else 5.0
    b_max = float(parsed.bfactor_stats["p95"]) if parsed.bfactor_stats else 80.0
    plddt_0_1 = bool(parsed.plddt_stats and str(parsed.plddt_stats["raw_scale"]).startswith("0-1"))

    config = ViewerConfig(
        render_style=render_style,
        color_mode=color_mode,
        file_format=parsed.file_format,
        source_label=parsed.pdb_id or parsed.filename,
        is_predicted=parsed.is_predicted,
        plddt_scale_0_1=plddt_0_1,
        bfactor_min=b_min,
        bfactor_max=b_max,
        selected_chain=selected_chain,
        highlight_residues=hi_ids,
        highlight_chain=selected_chain if selected_chain != "All" else None,
        show_surface=show_surface,
        surface_opacity=surface_opacity,
        show_ligands=show_ligands,
        show_labels=show_labels,
        zoom_to_selection=zoom_sel,
        atom_count=parsed.polymer_atom_count,
    )

    with viewer_col:
        html = build_viewer_html(parsed.raw_text, config)
        if hasattr(st, "iframe"):
            st.iframe(html, height=640)
        else:
            components.html(html, height=640, scrolling=False)
