"""Minimal scientific UI polish complementing .streamlit/config.toml."""

from __future__ import annotations

import streamlit as st


def inject_scientific_styles() -> None:
    """Apply subtle, restraint-first CSS for badges and scientific metric spacing."""
    st.markdown(
        """
        <style>
          .bf-header {
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            gap: 0.75rem;
            padding: 0.75rem 1.1rem;
            margin-bottom: 0.75rem;
            border: 1px solid #1e293b;
            border-radius: 12px;
            background: #111827;
          }
          .bf-title {
            font-size: 1.35rem;
            font-weight: 700;
            color: #f8fafc;
            letter-spacing: -0.02em;
            margin: 0;
          }
          .bf-subtitle {
            font-size: 0.84rem;
            color: #94a3b8;
            margin: 0.1rem 0 0;
          }
          .bf-badge-row {
            display: flex;
            align-items: center;
            gap: 0.5rem;
            flex-wrap: wrap;
          }
          .bf-badge {
            font-size: 0.74rem;
            font-weight: 600;
            padding: 0.25rem 0.65rem;
            border-radius: 999px;
            border: 1px solid #334155;
            background: #1e293b;
            color: #cbd5e1;
          }
          .bf-badge.exp {
            border-color: rgba(56, 189, 248, 0.4);
            background: rgba(56, 189, 248, 0.12);
            color: #7dd3fc;
          }
          .bf-badge.pred {
            border-color: rgba(167, 139, 250, 0.45);
            background: rgba(167, 139, 250, 0.14);
            color: #ddd6fe;
          }
          .bf-badge.warn {
            border-color: rgba(251, 191, 36, 0.4);
            background: rgba(251, 191, 36, 0.12);
            color: #fde68a;
          }
          footer { visibility: hidden; }
        </style>
        """,
        unsafe_allow_html=True,
    )
