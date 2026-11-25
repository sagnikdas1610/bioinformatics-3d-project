"""Centralized session-state invalidation, content fingerprinting, and active sequence context."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Literal

import streamlit as st

_DEPENDENT_STATE_KEYS: tuple[str, ...] = (
    "generated_pdf_bytes",
    "generated_pdf_fp",
    "generated_txt_report",
    "active_cm_params",
    "cm_chain_select",
    "struct_chain_scope",
    "viewer_focus_chain",
    "viewer_hi_text",
    "sidebar_structure_chain",
    "last_alignment",
    "alignment_active_in_session",
    "alignment_structure_fp",
    "alignment_inputs_snapshot",
)


@dataclass(frozen=True)
class ActiveSequenceContext:
    """Explicit provenance context for the active sequence used in analytics and reporting."""

    sequence: str
    label: str
    source: Literal["structure", "custom_fasta", "unavailable"]
    chain_id: str | None
    matches_structure: bool | None
    matched_chains: list[str]
    mismatch_warning: str | None
    extraction_note: str | None


def compute_structure_fingerprint(raw_text: str, filename: str, origin_hint: str | None) -> str:
    """Compute a deterministic SHA-256 content fingerprint for the active 3D structure."""
    payload = f"{filename}|{origin_hint or ''}|{len(raw_text)}|{raw_text}"
    return hashlib.sha256(payload.encode("utf-8", errors="replace")).hexdigest()


def compute_sequence_fingerprint(seq: str) -> str:
    """Compute a deterministic SHA-256 content fingerprint for an amino-acid sequence."""
    return hashlib.sha256(seq.encode("utf-8", errors="replace")).hexdigest()


def invalidate_structure_dependents() -> None:
    """Clear all structure-dependent derived artifacts, chain selections, and alignment session state."""
    for key in _DEPENDENT_STATE_KEYS:
        st.session_state.pop(key, None)


def set_active_structure(
    raw_text: str,
    filename: str,
    origin_hint: str | None,
    *,
    nav_target: str | None = None,
    reset_sequence_mode: bool = True,
) -> bool:
    """Set active 3D structure and purge stale derived analytics if structure content changed."""
    new_fp = compute_structure_fingerprint(raw_text, filename, origin_hint)
    prev_fp = st.session_state.get("structure_fingerprint")
    changed = new_fp != prev_fp

    st.session_state["structure_raw"] = raw_text
    st.session_state["structure_filename"] = filename
    st.session_state["structure_origin_hint"] = origin_hint
    st.session_state["structure_fingerprint"] = new_fp
    st.session_state["structure_error"] = None

    if changed:
        invalidate_structure_dependents()
        if reset_sequence_mode:
            st.session_state["sequence_source_mode"] = "From Structure"

    if nav_target is not None:
        st.session_state["main_nav_section"] = nav_target

    return changed


def sync_structure_fingerprint_guard() -> str:
    """Detect out-of-band structure state edits (e.g., tests or callbacks) and purge stale dependents."""
    raw_text = st.session_state.get("structure_raw", "")
    filename = st.session_state.get("structure_filename", "")
    origin_hint = st.session_state.get("structure_origin_hint")
    live_fp = compute_structure_fingerprint(raw_text, filename, origin_hint)
    prev_fp = st.session_state.get("structure_fingerprint")

    if live_fp != prev_fp:
        invalidate_structure_dependents()
        st.session_state["structure_fingerprint"] = live_fp
        st.session_state["sequence_source_mode"] = "From Structure"

    return live_fp
