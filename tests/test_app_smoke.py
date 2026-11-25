"""Headless Streamlit AppTest smoke tests across all workspace views."""

from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest

APP_PATH = Path(__file__).resolve().parent.parent / "app.py"


def test_app_boots_and_navigates_all_sections():
    at = AppTest.from_file(str(APP_PATH), default_timeout=20)
    at.run()
    assert not at.exception

    for section in (
        "Sequence Analytics",
        "Pairwise Alignment",
        "Structural Analysis",
        "AI Structure Prediction",
        "Reports & Exports",
        "3D Explorer",
    ):
        nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
        nav.set_value(section).run()
        assert not at.exception, f"Exception on section {section}: {at.exception}"

    # Verify all Structural Analysis sub-views render cleanly
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Structural Analysis").run()
    for sub in ("Structure Quality", "Ramachandran Plot", "Secondary Structure", "Contact Map"):
        sub_ctrl = next(w for w in at.segmented_control if w.key != "main_nav_section" and w.label == "Analysis module")
        sub_ctrl.set_value(sub).run()
        assert not at.exception, f"Exception in Structural Analysis/{sub}: {at.exception}"

    # Verify PDF report button generates a valid PDF in session state
    nav = next(w for w in at.segmented_control if w.key == "main_nav_section")
    nav.set_value("Reports & Exports").run()
    gen_btn = next(b for b in at.button if b.key == "btn_gen_pdf")
    gen_btn.click().run()
    assert not at.exception
    pdf_bytes = at.session_state["generated_pdf_bytes"]
    assert isinstance(pdf_bytes, bytes)
    assert pdf_bytes.startswith(b"%PDF-")
