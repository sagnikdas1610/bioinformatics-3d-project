"""BioFold Workbench — Bioinformatics 3D Protein Structure & Analytics Platform.

Preserved root entry point for Streamlit Community Cloud deployment.
"""

from __future__ import annotations

from html import escape
from pathlib import Path

import streamlit as st

from bio3d.config import APP_NAME, APP_VERSION, BUNDLED_STRUCTURES, ORIGIN_LABELS
from bio3d.sequence import SequenceProfile, analyze_sequence
from bio3d.structural_analysis import (
    ContactMapResult,
    RamachandranResult,
    SecondaryStructureResult,
    compute_contact_map,
    compute_ramachandran,
    compute_secondary_structure,
)
from bio3d.structure import ParsedStructure, parse_structure
from ui.sidebar import render_sidebar
from ui.styles import inject_scientific_styles
from ui.tabs.alignment_tab import render_alignment_tab
from ui.tabs.prediction_tab import render_prediction_tab
from ui.tabs.report_tab import render_report_tab
from ui.tabs.sequence_tab import render_sequence_tab
from ui.tabs.structure_tab import render_structure_tab
from ui.tabs.viewer_tab import render_viewer_tab

st.set_page_config(
    page_title=f"{APP_NAME} | Protein Structure & Analytics",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_data(show_spinner=False, max_entries=16)
def cached_parse_structure(raw_text: str, filename: str, origin_hint: str | None) -> ParsedStructure:
    return parse_structure(raw_text, filename=filename, origin_hint=origin_hint)


@st.cache_data(show_spinner=False, max_entries=32)
def cached_analyze_sequence(seq: str) -> SequenceProfile:
    return analyze_sequence(seq)


@st.cache_data(show_spinner=False, max_entries=16)
def cached_ramachandran(raw_text: str, filename: str, chain_id: str) -> RamachandranResult:
    parsed = cached_parse_structure(raw_text, filename, None)
    return compute_ramachandran(parsed, chain_id=None if chain_id == "All" else chain_id)


@st.cache_data(show_spinner=False, max_entries=16)
def cached_secondary_structure(raw_text: str, filename: str, chain_id: str) -> SecondaryStructureResult:
    parsed = cached_parse_structure(raw_text, filename, None)
    return compute_secondary_structure(parsed, chain_id=None if chain_id == "All" else chain_id)


@st.cache_data(show_spinner=False, max_entries=24)
def cached_contact_map(
    raw_text: str,
    filename: str,
    chain_id: str,
    atom_mode: str,
    threshold: float,
    min_seq_sep: int,
) -> ContactMapResult:
    parsed = cached_parse_structure(raw_text, filename, None)
    return compute_contact_map(
        parsed,
        chain_id=chain_id,
        atom_mode=atom_mode,  # type: ignore[arg-type]
        threshold_angstrom=threshold,
        min_seq_sep=min_seq_sep,
    )


def _init_default_structure() -> None:
    if "structure_raw" in st.session_state:
        return
    default_entry = next(iter(BUNDLED_STRUCTURES.values()))
    default_path = Path(default_entry["path"])
    if default_path.exists():
        st.session_state["structure_raw"] = default_path.read_text(encoding="utf-8", errors="replace")
        st.session_state["structure_filename"] = default_path.name
        st.session_state["structure_origin_hint"] = "experimental"
        st.session_state["structure_error"] = None
    else:
        st.session_state["structure_raw"] = ""
        st.session_state["structure_filename"] = ""
        st.session_state["structure_origin_hint"] = None
        st.session_state["structure_error"] = f"Missing bundled sample structure: {default_path.name}"


inject_scientific_styles()
_init_default_structure()

parsed_structure: ParsedStructure | None = None
raw_struct = st.session_state.get("structure_raw", "")
fname_struct = st.session_state.get("structure_filename", "structure.pdb")
origin_hint = st.session_state.get("structure_origin_hint")

if raw_struct:
    try:
        parsed_structure = cached_parse_structure(raw_struct, fname_struct, origin_hint)
    except Exception as exc:
        st.session_state["structure_error"] = f"Structure parse failed: {exc}"
        parsed_structure = None

active_sequence, active_seq_label = render_sidebar(parsed_structure)

# ── Compact scientific header ──────────────────────────────────────────────────
if parsed_structure is not None:
    badge_cls = "pred" if parsed_structure.is_predicted else "exp"
    origin_text = escape(ORIGIN_LABELS.get(parsed_structure.structure_origin, parsed_structure.structure_origin))
    struct_badge = escape(
        f"{parsed_structure.pdb_id or parsed_structure.filename} · "
        f"{parsed_structure.polymer_residue_count} aa · "
        f"{len(parsed_structure.chains)} chain(s)"
    )
    badges_html = (
        f'<span class="bf-badge {badge_cls}">{origin_text}</span>'
        f'<span class="bf-badge">{struct_badge}</span>'
    )
else:
    badges_html = '<span class="bf-badge warn">No structure loaded</span>'

st.markdown(
    f"""
    <div class="bf-header">
      <div>
        <div class="bf-title">{APP_NAME} <span style="font-size:0.75rem;color:#64748b;font-weight:500;">v{APP_VERSION}</span></div>
        <div class="bf-subtitle">Interactive 3D structural biology, sequence analytics, pairwise alignment, and stereochemical diagnostics.</div>
      </div>
      <div class="bf-badge-row">{badges_html}</div>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Primary workspace navigation (renders only the active section) ─────────────
NAV_SECTIONS = [
    "3D Explorer",
    "Sequence Analytics",
    "Pairwise Alignment",
    "Structural Analysis",
    "AI Structure Prediction",
    "Reports & Exports",
]

active_section = st.segmented_control(
    "Workspace section",
    options=NAV_SECTIONS,
    default="3D Explorer",
    label_visibility="collapsed",
    key="main_nav_section",
)
if active_section is None:
    active_section = "3D Explorer"

if active_section == "3D Explorer":
    render_viewer_tab(parsed_structure)
elif active_section == "Sequence Analytics":
    render_sequence_tab(active_sequence, active_seq_label, parsed_structure, cached_analyze_sequence)
elif active_section == "Pairwise Alignment":
    render_alignment_tab(active_sequence)
elif active_section == "Structural Analysis":
    render_structure_tab(parsed_structure, cached_ramachandran, cached_secondary_structure, cached_contact_map)
elif active_section == "AI Structure Prediction":
    render_prediction_tab(active_sequence)
else:
    render_report_tab(
        parsed_structure,
        active_sequence,
        active_seq_label,
        cached_analyze_sequence,
        cached_ramachandran,
        cached_secondary_structure,
        cached_contact_map,
    )
