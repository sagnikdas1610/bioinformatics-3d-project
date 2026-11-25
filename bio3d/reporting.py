"""Scientific PDF and plain-text report generation using ReportLab and Matplotlib (Agg)."""

from __future__ import annotations

from datetime import datetime, timezone
from io import BytesIO
import os
from pathlib import Path
import tempfile
from xml.sax.saxutils import escape as xml_escape

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "biofold_mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import LETTER  # noqa: E402
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet  # noqa: E402
from reportlab.lib.units import inch  # noqa: E402
from reportlab.pdfbase import pdfmetrics  # noqa: E402
from reportlab.pdfbase.pdfmetrics import registerFontFamily  # noqa: E402
from reportlab.pdfbase.ttfonts import TTFont  # noqa: E402
from reportlab.platypus import (  # noqa: E402
    HRFlowable,
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from bio3d.alignment import AlignmentResult
from bio3d.config import APP_NAME, APP_VERSION, ORIGIN_LABELS, STANDARD_AA
from bio3d.sequence import SequenceProfile
from bio3d.structural_analysis import (
    ContactMapResult,
    RamachandranResult,
    SecondaryStructureResult,
    rama_background_grid,
)
from bio3d.structure import ParsedStructure

USABLE_WIDTH = 7.20 * inch  # LETTER (8.5 in) - 2 * 0.65 in margin


def _register_unicode_fonts() -> tuple[str, str]:
    """Register Unicode-capable DejaVu Sans fonts bundled with Matplotlib."""
    font_dir = Path(matplotlib.__file__).resolve().parent / "mpl-data" / "fonts" / "ttf"
    reg_ttf = font_dir / "DejaVuSans.ttf"
    bold_ttf = font_dir / "DejaVuSans-Bold.ttf"
    if reg_ttf.exists() and bold_ttf.exists():
        if "BioSans" not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont("BioSans", str(reg_ttf)))
        if "BioSans-Bold" not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont("BioSans-Bold", str(bold_ttf)))
        registerFontFamily("BioSans", normal="BioSans", bold="BioSans-Bold")
        return "BioSans", "BioSans-Bold"
    return "Helvetica", "Helvetica-Bold"


FONT_REG, FONT_BOLD = _register_unicode_fonts()


def _build_pdf_styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "h1": ParagraphStyle(
            "BFH1",
            parent=base["Heading1"],
            fontName=FONT_BOLD,
            fontSize=15.5,
            leading=19,
            textColor=colors.HexColor("#0f172a"),
            spaceAfter=4,
        ),
        "h2": ParagraphStyle(
            "BFH2",
            parent=base["Heading2"],
            fontName=FONT_BOLD,
            fontSize=11,
            leading=14.5,
            textColor=colors.HexColor("#0369a1"),
            spaceBefore=10,
            spaceAfter=5,
            keepWithNext=True,
        ),
        "body": ParagraphStyle(
            "BFBody",
            parent=base["Normal"],
            fontName=FONT_REG,
            fontSize=8.8,
            leading=12.5,
            textColor=colors.HexColor("#1e293b"),
        ),
        "warn": ParagraphStyle(
            "BFWarn",
            parent=base["Normal"],
            fontName=FONT_BOLD,
            fontSize=8.5,
            leading=12.0,
            textColor=colors.HexColor("#b45309"),
        ),
        "muted": ParagraphStyle(
            "BFMuted",
            parent=base["Normal"],
            fontName=FONT_REG,
            fontSize=8.0,
            leading=11.5,
            textColor=colors.HexColor("#475569"),
        ),
        "th": ParagraphStyle(
            "BFTableHeader",
            parent=base["Normal"],
            fontName=FONT_BOLD,
            fontSize=8.2,
            leading=10.8,
            textColor=colors.white,
        ),
        "td_label": ParagraphStyle(
            "BFTableLabel",
            parent=base["Normal"],
            fontName=FONT_BOLD,
            fontSize=8.1,
            leading=11.0,
            textColor=colors.HexColor("#0f172a"),
            wordWrap="LTR",
        ),
        "td_val": ParagraphStyle(
            "BFTableValue",
            parent=base["Normal"],
            fontName=FONT_REG,
            fontSize=8.1,
            leading=11.0,
            textColor=colors.HexColor("#1e293b"),
            wordWrap="LTR",
        ),
        "td_prov": ParagraphStyle(
            "BFTableProv",
            parent=base["Normal"],
            fontName=FONT_BOLD,
            fontSize=7.6,
            leading=10.5,
            textColor=colors.HexColor("#334155"),
            wordWrap="LTR",
        ),
    }


def _safe_xml(text: object) -> str:
    return xml_escape(str(text))


def _styled_table(
    rows: list[list[str]],
    col_widths: list[float],
    styles: dict[str, ParagraphStyle],
    markup_cols: set[int] | None = None,
) -> Table:
    """Create a ReportLab table where every cell is a wrapped Paragraph to prevent overlap."""
    markup_cols = markup_cols or set()
    wrapped_rows: list[list[Paragraph]] = []
    for r_idx, row in enumerate(rows):
        wrapped_row: list[Paragraph] = []
        for c_idx, cell in enumerate(row):
            if r_idx == 0:
                p_style = styles["th"]
                content = _safe_xml(cell)
            else:
                if c_idx == 0:
                    p_style = styles["td_label"]
                elif c_idx == len(row) - 1 and len(row) == 3:
                    p_style = styles["td_prov"]
                else:
                    p_style = styles["td_val"]
                content = str(cell) if c_idx in markup_cols else _safe_xml(cell)
            wrapped_row.append(Paragraph(content, p_style))
        wrapped_rows.append(wrapped_row)

    table = Table(wrapped_rows, colWidths=col_widths, repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0f172a")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8fafc")]),
                ("GRID", (0, 0), (-1, -1), 0.45, colors.HexColor("#cbd5e1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    return table


def _fig_to_png_bytes(fig: plt.Figure, dpi: int = 160) -> BytesIO:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return buf


def _make_backbone_3d_figure(parsed: ParsedStructure) -> BytesIO | None:
    """Render a deterministic 3D C-alpha backbone trace figure for the PDF report."""
    has_ca = any("CA" in r.atoms for chain_res in parsed.polymer_residues.values() for r in chain_res)
    if not has_ca:
        return None

    fig = plt.figure(figsize=(5.4, 3.3))
    ax = fig.add_subplot(111, projection="3d")
    palette = ["#0284c7", "#7c3aed", "#059669", "#d97706", "#dc2626", "#0d9488"]

    plotted_chains = 0
    for idx, ch in enumerate(parsed.chains):
        pts = [r.atoms["CA"].coord for r in parsed.polymer_residues[ch] if "CA" in r.atoms]
        if len(pts) < 2:
            continue
        arr = np.vstack(pts)
        ax.plot(arr[:, 0], arr[:, 1], arr[:, 2], linewidth=1.8, color=palette[idx % len(palette)], label=f"Chain {ch}")
        ax.scatter(arr[0, 0], arr[0, 1], arr[0, 2], color="#16a34a", s=24)
        plotted_chains += 1

    ax.set_title(f"Cα Backbone Trace — {parsed.pdb_id or parsed.filename}", fontsize=10, fontweight="bold")
    ax.set_xlabel("X (Å)", fontsize=8)
    ax.set_ylabel("Y (Å)", fontsize=8)
    ax.set_zlabel("Z (Å)", fontsize=8)
    ax.tick_params(labelsize=7)
    if 0 < plotted_chains <= 6:
        ax.legend(fontsize=7, loc="upper right")
    return _fig_to_png_bytes(fig)


def _make_sequence_figure(
    profile: SequenceProfile,
    hydropathy: tuple[list[int], list[float]] | None,
) -> BytesIO:
    """Render amino-acid composition and Kyte-Doolittle hydropathy figures."""
    ncols = 2 if hydropathy and len(hydropathy[0]) > 1 else 1
    fig, axes = plt.subplots(1, ncols, figsize=(7.0 if ncols == 2 else 5.2, 2.55))
    ax1 = axes[0] if ncols == 2 else axes

    aas = list(STANDARD_AA)
    counts = [profile.aa_counts[a] for a in aas]
    ax1.bar(aas, counts, color="#0284c7", width=0.72)
    ax1.set_title("Amino-Acid Composition", fontsize=9.5, fontweight="bold")
    ax1.set_ylabel("Count", fontsize=8)
    ax1.tick_params(axis="x", labelsize=7)
    ax1.tick_params(axis="y", labelsize=7)
    ax1.grid(axis="y", linestyle=":", alpha=0.45)

    if ncols == 2 and hydropathy is not None:
        pos, scores = hydropathy
        ax2 = axes[1]
        ax2.axhline(0.0, color="#64748b", linestyle=":", linewidth=0.9)
        ax2.plot(pos, scores, color="#0d9488", linewidth=1.4)
        ax2.fill_between(pos, scores, 0.0, alpha=0.15, color="#0d9488")
        ax2.set_title("Kyte-Doolittle Hydropathy", fontsize=9.5, fontweight="bold")
        ax2.set_xlabel("Residue Position", fontsize=8)
        ax2.set_ylabel("Score", fontsize=8)
        ax2.tick_params(labelsize=7)
        ax2.grid(True, linestyle=":", alpha=0.35)

    fig.tight_layout()
    return _fig_to_png_bytes(fig)


def _make_ramachandran_figure(rama: RamachandranResult) -> BytesIO:
    """Render Matplotlib Ramachandran plot with MolProbity Top8000 contours."""
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    phis_g, psis_g, z_g = rama_background_grid("General")
    ax.contourf(phis_g, psis_g, z_g, levels=[-0.5, 0.5, 1.5, 2.5], colors=["#f8fafc", "#dbeafe", "#93c5fd"], alpha=0.8)
    ax.axhline(0, color="#94a3b8", linewidth=0.6)
    ax.axvline(0, color="#94a3b8", linewidth=0.6)

    for reg, col, marker, sz in (("Favored", "#0369a1", "o", 18), ("Allowed", "#d97706", "D", 24), ("Outlier", "#dc2626", "x", 36)):
        pts = [r for r in rama.residues if r.region == reg]
        if pts:
            ax.scatter([p.phi for p in pts], [p.psi for p in pts], c=col, marker=marker, s=sz, label=f"{reg} ({len(pts)})", zorder=4)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_xlabel("φ (deg)", fontsize=8.5)
    ax.set_ylabel("ψ (deg)", fontsize=8.5)
    ax.set_title("Ramachandran Plot (MolProbity Top8000 General)", fontsize=9.2, fontweight="bold")
    ax.tick_params(labelsize=7.5)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    return _fig_to_png_bytes(fig)


def _make_ss_and_contact_figure(
    ss: SecondaryStructureResult | None,
    cmap: ContactMapResult | None,
) -> BytesIO | None:
    """Render secondary structure distribution and/or contact map for the PDF."""
    if ss is None and cmap is None:
        return None
    ncols = 2 if (ss is not None and cmap is not None) else 1
    fig, axes = plt.subplots(1, ncols, figsize=(7.0 if ncols == 2 else 4.8, 2.9))

    idx = 0
    if ss is not None:
        ax_ss = axes[idx] if ncols == 2 else axes
        idx += 1
        cats = ["Helix\n(H/G/I)", "Sheet\n(E/B)", "Coil/Turn\n(T/C)"]
        vals = [ss.helix_pct, ss.sheet_pct, ss.coil_pct]
        bars = ax_ss.bar(cats, vals, color=["#db2777", "#d97706", "#0284c7"], width=0.55)
        for bar, v in zip(bars, vals):
            ax_ss.text(bar.get_x() + bar.get_width() / 2, v + 1.5, f"{v:.1f}%", ha="center", fontsize=8)
        ax_ss.set_ylim(0, max(vals + [20.0]) * 1.18)
        ax_ss.set_ylabel("Percentage (%)", fontsize=8)
        ax_ss.set_title("Secondary Structure Distribution", fontsize=9.5, fontweight="bold")
        ax_ss.tick_params(labelsize=8)

    if cmap is not None:
        ax_cm = axes[idx] if ncols == 2 else axes
        masked = np.where(cmap.contact_matrix, cmap.distance_matrix, np.nan)
        im = ax_cm.imshow(masked, cmap="viridis_r", vmin=2.5, vmax=cmap.threshold_angstrom, origin="lower")
        fig.colorbar(im, ax=ax_cm, fraction=0.046, pad=0.04, label="Dist (Å)")
        ax_cm.set_title(f"Contact Map ({cmap.atom_mode}, ≤{cmap.threshold_angstrom:.1f} Å)", fontsize=9.5, fontweight="bold")
        ax_cm.set_xlabel("Residue Index", fontsize=8)
        ax_cm.set_ylabel("Residue Index", fontsize=8)
        ax_cm.tick_params(labelsize=7)

    fig.tight_layout()
    return _fig_to_png_bytes(fig)


def _draw_page_footer(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont(FONT_REG, 7.5)
    canvas.setFillColor(colors.HexColor("#64748b"))
    canvas.drawString(0.65 * inch, 0.40 * inch, f"{APP_NAME} v{APP_VERSION} — Computational Structural & Sequence Report")
    canvas.drawRightString(LETTER[0] - 0.65 * inch, 0.40 * inch, f"Page {doc.page}")
    canvas.restoreState()


def _format_seq_source_label(sequence_source: str | None, report_chain: str | None) -> str:
    if sequence_source == "structure":
        ch = f" (Chain {report_chain})" if report_chain else ""
        return f"Derived from active 3D structure coordinates{ch}"
    if sequence_source == "custom_fasta":
        return "Independent uploaded / manual FASTA input"
    if sequence_source == "unavailable":
        return "Unavailable"
    return "Active session sequence"


def generate_pdf_report(
    parsed: ParsedStructure | None,
    seq_profile: SequenceProfile | None,
    hydropathy: tuple[list[int], list[float]] | None = None,
    alignment: AlignmentResult | None = None,
    rama: RamachandranResult | None = None,
    ss: SecondaryStructureResult | None = None,
    contact_map: ContactMapResult | None = None,
    report_title: str | None = None,
    sequence_name: str | None = None,
    sequence_source: str | None = None,
    sequence_mismatch_warning: str | None = None,
    sequence_unavailable_reason: str | None = None,
    report_chain: str | None = None,
) -> bytes:
    """Build a publication-style scientific PDF report with wrapped tables and Unicode typography."""
    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=LETTER,
        leftMargin=0.65 * inch,
        rightMargin=0.65 * inch,
        topMargin=0.60 * inch,
        bottomMargin=0.65 * inch,
        title=f"{APP_NAME} Scientific Report",
    )

    styles = _build_pdf_styles()
    col3_widths = [1.85 * inch, 3.95 * inch, 1.40 * inch]  # sum = 7.20 in

    story: list[object] = []
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    story.append(Paragraph(_safe_xml(report_title or f"{APP_NAME} — Structural & Sequence Analysis Report"), styles["h1"]))
    story.append(
        Paragraph(
            f"Generated: {timestamp} | Platform v{APP_VERSION} | Provenance column distinguishes Experimental (Measured), Calculated, AI-Predicted, Reference, and Unavailable data.",
            styles["muted"],
        )
    )
    story.append(Spacer(1, 5))
    story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#cbd5e1")))
    story.append(Spacer(1, 8))

    # 1. Structure Identification & Quality
    sec1_head = Paragraph("1. Structure Identification & Quality Summary", styles["h2"])
    if parsed is not None:
        origin_label = ORIGIN_LABELS.get(parsed.structure_origin, parsed.structure_origin)
        origin_prov = "AI-Predicted" if parsed.is_predicted else ("Experimental (Measured)" if parsed.experimental_method else "File Record")
        chain_breakdown = ", ".join(f"{ch}: {len(res)} aa" for ch, res in parsed.polymer_residues.items())
        struct_rows = [
            ["Field", "Value", "Provenance"],
            ["Source File", f"{_safe_xml(parsed.filename)} ({parsed.file_format.upper()})", "File Record"],
            ["Origin Classification", _safe_xml(origin_label), origin_prov],
            [
                "PDB / Entry ID",
                _safe_xml(parsed.pdb_id) if parsed.pdb_id else "Unavailable in file header",
                "Experimental (Measured)" if (parsed.pdb_id and not parsed.is_predicted) else ("AI-Predicted" if parsed.pdb_id else "Unavailable"),
            ],
            [
                "Experimental Method",
                _safe_xml(parsed.experimental_method) if parsed.experimental_method else ("Computational model (not experimental)" if parsed.is_predicted else "Unavailable"),
                "Experimental (Measured)" if parsed.experimental_method else "Unavailable",
            ],
            [
                "Resolution",
                f"{parsed.resolution_angstrom:.2f} Å" if parsed.resolution_angstrom is not None else "Unavailable (not reported / non-diffraction)",
                "Experimental (Measured)" if parsed.resolution_angstrom is not None else "Unavailable",
            ],
            [
                "Polymer Chains / Residues",
                f"{len(parsed.chains)} ({_safe_xml(chain_breakdown)}) | Total: {parsed.polymer_residue_count} aa",
                "Calculated (Structure)",
            ],
            ["Polymer Atoms (Model 1)", f"{parsed.polymer_atom_count:,} protein atoms (models in file: {parsed.model_count})", "Calculated (Structure)"],
            ["Water Molecules (HOH)", f"{parsed.water_count:,} waters", "Calculated (Structure)"],
            [
                "Hetero Ligands (Non-Water)",
                f"{parsed.hetero_ligand_count:,} ({_safe_xml(', '.join(parsed.hetero_ligand_names) or 'none')})",
                "Calculated (Structure)",
            ],
            ["Missing Residues Check", _safe_xml(parsed.missing_residue_note), "File / Calculated"],
        ]
        if parsed.modified_residues_mapped or parsed.unmapped_polymer_residues:
            notes = "; ".join(parsed.chain_extraction_notes.values())
            struct_rows.append(["Polymer Residue Mapping", _safe_xml(notes), "Calculated (Structure)"])

        if parsed.is_predicted and parsed.plddt_stats:
            ps = parsed.plddt_stats
            struct_rows.append([
                "pLDDT Confidence (0–100)",
                f"Mean {ps['mean']:.1f} | Median {ps['median']:.1f} | ≥90: {ps['pct_very_high']:.1f}% | &lt;50: {ps['pct_very_low']:.1f}% (raw scale: {_safe_xml(ps['raw_scale'])})",
                "AI-Predicted",
            ])
        elif parsed.bfactor_stats:
            bs = parsed.bfactor_stats
            b_summary = (
                f"Cα: mean {bs['mean']:.2f}, median {bs['median']:.2f}, range {bs['min']:.2f}–{bs['max']:.2f} Å<sup>2</sup>"
            )
            if "all_atom_mean" in bs:
                b_summary += (
                    f" | All polymer atoms ({bs['all_atom_count']:,}): mean {bs['all_atom_mean']:.2f}, "
                    f"range {bs['all_atom_min']:.2f}–{bs['all_atom_max']:.2f} Å<sup>2</sup> "
                    f"(Occupancy cols 55–60: {bs['occupancy_min']:.2f}–{bs['occupancy_max']:.2f})"
                )
            struct_rows.append([
                "Experimental B-factor",
                b_summary,
                "Experimental (Measured)",
            ])
        sec1_table = _styled_table(struct_rows, col3_widths, styles, markup_cols={1})
        story.append(KeepTogether([sec1_head, sec1_table]))
        backbone_png = _make_backbone_3d_figure(parsed)
        if backbone_png:
            story.append(Spacer(1, 8))
            story.append(Image(backbone_png, width=4.5 * inch, height=2.75 * inch))
    else:
        story.append(KeepTogether([sec1_head, Paragraph("[Unavailable] No 3D structure was loaded during this session.", styles["body"])]))

    # 2. Primary Sequence Statistics
    sec2_head = Paragraph("2. Primary Sequence & Molecular Properties", styles["h2"])
    src_desc = _format_seq_source_label(sequence_source, report_chain)
    if seq_profile is not None:
        match_status = "No structure loaded"
        if parsed is not None:
            matched = [ch for ch, s in parsed.chain_sequences.items() if s == seq_profile.sequence]
            match_status = f"Exact match to Chain {', '.join(matched)}" if matched else f"Independent sequence (differs from {parsed.filename} chains)"

        seq_rows = [
            ["Descriptor", "Value", "Provenance"],
            ["Sequence Source", _safe_xml(src_desc), "Session Provenance"],
            ["Sequence Label / Match", f"{_safe_xml(sequence_name or 'Active sequence')} — {_safe_xml(match_status)}", "File / Calculated"],
            ["Length", f"{seq_profile.length} amino acids", "Calculated (Sequence)"],
            ["Molecular Weight (avg, red. Cys)", f"{seq_profile.molecular_weight_da:,.2f} Da ({seq_profile.molecular_weight_kda:.2f} kDa)", "Calculated (Sequence)"],
            ["Theoretical Isoelectric Point (pI)", f"{seq_profile.theoretical_pi:.2f}", "Calculated (Sequence)"],
            ["GRAVY Hydropathy (Kyte-Doolittle)", f"{seq_profile.gravy:.3f}", "Calculated (Sequence)"],
            ["Instability Index (Guruprasad)", f"{seq_profile.instability_index:.2f} ({_safe_xml(seq_profile.instability_class)})", "Calculated (Sequence)"],
            ["Aromaticity (Phe+Trp+Tyr)", f"{seq_profile.aromaticity:.3f}", "Calculated (Sequence)"],
            ["Molar Extinction (280 nm, red / ox)", f"{seq_profile.extinction_reduced:,} / {seq_profile.extinction_oxidized:,} M<sup>-1</sup> cm<sup>-1</sup>", "Calculated (Sequence)"],
            ["Estimated Net Charge at pH 7.0", f"{seq_profile.net_charge_ph7:+.2f}", "Calculated (Sequence)"],
        ]
        if sequence_mismatch_warning:
            seq_rows.insert(3, ["Sequence–Structure Warning", f"<b>{_safe_xml(sequence_mismatch_warning)}</b>", "Integrity Check"])

        sec2_table = _styled_table(seq_rows, col3_widths, styles, markup_cols={1})
        story.append(KeepTogether([sec2_head, sec2_table]))
        story.append(Spacer(1, 8))
        seq_fig = _make_sequence_figure(seq_profile, hydropathy)
        story.append(Image(seq_fig, width=6.6 * inch, height=2.40 * inch))
    else:
        reason = sequence_unavailable_reason or "No valid standard amino-acid sequence was analyzed."
        story.append(KeepTogether([sec2_head, Paragraph(f"[UNAVAILABLE] {_safe_xml(reason)}", styles["body"])]))

    # 3. Pairwise Alignment
    sec3_head = Paragraph("3. Pairwise Sequence Alignment", styles["h2"])
    if alignment is not None:
        s1_acc = alignment.seq1_accession or "None (unannotated input)"
        s2_acc = alignment.seq2_accession or "None (unannotated input)"
        aln_rows = [
            ["Metric", "Value", "Provenance"],
            ["Sequence A", f"{alignment.seq1_name} | Accession: {s1_acc} | Length: {alignment.seq1_length} aa (aligned span {alignment.seq1_span[0]}–{alignment.seq1_span[1]})", "Input Metadata"],
            ["Sequence B", f"{alignment.seq2_name} | Accession: {s2_acc} | Length: {alignment.seq2_length} aa (aligned span {alignment.seq2_span[0]}–{alignment.seq2_span[1]})", "Input Metadata"],
            ["Alignment Mode", f"{alignment.mode.capitalize()} ({'Needleman–Wunsch' if alignment.mode == 'global' else 'Smith–Waterman'})", "Calculated (Alignment)"],
            ["Scoring Matrix & Gap Penalties", f"{alignment.matrix_name} | Gap open = {alignment.open_gap_score} | Gap extend = {alignment.extend_gap_score}", "Calculated (Alignment)"],
            ["Alignment Score", f"{alignment.score:.1f} (Engine: {alignment.engine})", "Calculated (Alignment)"],
            ["Aligned Length", f"{alignment.alignment_length} columns", "Calculated (Alignment)"],
            ["Sequence Identity (aligned cols)", f"{alignment.identity_pct:.2f}% ({alignment.identical_count}/{alignment.alignment_length}; ungapped: {alignment.identity_ungapped_pct:.2f}%)", "Calculated (Alignment)"],
            ["Sequence Similarity (score > 0)", f"{alignment.similarity_pct:.2f}% ({alignment.similar_count}/{alignment.alignment_length})", "Calculated (Alignment)"],
            ["Gap Count & Percentage", f"{alignment.gap_count} gap columns ({alignment.gap_pct:.2f}%)", "Calculated (Alignment)"],
        ]
        sec3_table = _styled_table(aln_rows, col3_widths, styles)
        story.append(KeepTogether([sec3_head, sec3_table, Spacer(1, 4), Paragraph(_safe_xml(alignment.interpretation), styles["muted"])]))
    else:
        story.append(
            KeepTogether(
                [
                    sec3_head,
                    Paragraph(
                        "[UNAVAILABLE] No pairwise sequence alignment was performed for the current structure/report session.",
                        styles["body"],
                    ),
                ]
            )
        )

    # 4. Advanced Structural Analysis (start cleanly on new page when full structural data is present)
    if rama is not None or ss is not None or contact_map is not None:
        story.append(PageBreak())

    chain_scope_suffix = f" (Chain {report_chain})" if report_chain and report_chain != "All" else ""
    sec4_head = Paragraph(f"4. Structural Stereochemistry, Secondary Structure & Contacts{_safe_xml(chain_scope_suffix)}", styles["h2"])
    if rama is not None:
        rama_rows = [
            ["Ramachandran Metric", "Value", "Provenance"],
            ["Reference Distribution", rama.reference_source, "Reference (MolProbity)"],
            ["Chain Scope & Evaluated Residues", f"{rama.total_evaluated} internal residues evaluated (excluded termini/breaks: {rama.excluded_count})", "Calculated (Structure)"],
            ["Favored / Allowed / Outliers", f"{rama.favored_pct:.1f}% ({rama.favored_count}) / {rama.allowed_pct:.1f}% ({rama.allowed_count}) / {rama.outlier_pct:.1f}% ({rama.outlier_count})", "Calculated (Structure)"],
            ["Stereochemical Assessment", rama.interpretation, "Calculated (Structure)"],
        ]
        rama_table = _styled_table(rama_rows, col3_widths, styles)
        story.append(KeepTogether([sec4_head, rama_table, Spacer(1, 6), Image(_make_ramachandran_figure(rama), width=3.95 * inch, height=3.40 * inch)]))
        if rama.outliers:
            out_rows = [["Chain", "Residue", "MolProbity Class", "φ (deg)", "ψ (deg)"]]
            for o in rama.outliers[:15]:
                out_rows.append([o.chain, f"{o.resname} {o.resseq}", o.rama_class, f"{o.phi:.1f}", f"{o.psi:.1f}"])
            out_table = _styled_table(
                out_rows,
                [0.75 * inch, 1.55 * inch, 1.80 * inch, 1.55 * inch, 1.55 * inch],
                styles,
            )
            story.append(Spacer(1, 6))
            story.append(KeepTogether([Paragraph("Ramachandran Outliers", styles["td_label"]), Spacer(1, 3), out_table]))
    else:
        story.append(KeepTogether([sec4_head, Paragraph("[UNAVAILABLE] Ramachandran analysis was not computed.", styles["body"])]))

    if ss is not None:
        story.append(Spacer(1, 8))
        ss_rows = [
            ["Secondary Structure Metric", "Value", "Provenance"],
            ["Method", ss.method, "Calculated (Structure)"],
            ["Residues Analyzed", f"{ss.total_residues} aa", "Calculated (Structure)"],
            ["Helix (H/G/I) / Sheet (E/B) / Coil & Turn", f"{ss.helix_pct:.1f}% ({ss.helix_count}) / {ss.sheet_pct:.1f}% ({ss.sheet_count}) / {ss.coil_pct:.1f}% ({ss.coil_count})", "Calculated (Structure)"],
        ]
        if ss.header_available and ss.header_helix_pct is not None:
            ss_rows.append([
                "File Header HELIX/SHEET Annotation",
                f"Helix {ss.header_helix_pct:.1f}% / Sheet {ss.header_sheet_pct:.1f}% / Coil {ss.header_coil_pct:.1f}%",
                "Experimental (Measured)",
            ])
        story.append(_styled_table(ss_rows, col3_widths, styles))

    if contact_map is not None:
        story.append(Spacer(1, 8))
        gly_note = (
            f" (Glycine Cα fallback applied to {contact_map.glycine_fallback_count} Gly residues)"
            if contact_map.atom_mode == "CB"
            else ""
        )
        cm_rows = [
            ["Contact Map Metric", "Value", "Provenance"],
            ["Chain & Mode", f"Chain {contact_map.chain_id} ({contact_map.residue_count} residues) | {contact_map.atom_mode} atoms{gly_note}", "Calculated (Structure)"],
            ["Distance Cutoff & Seq Separation", f"<= {contact_map.threshold_angstrom:.1f} Å | |i - j| >= {contact_map.min_seq_sep} (polymer index)", "Calculated (Structure)"],
            ["Total Contacts / Long-Range (|i-j| >= 12)", f"{contact_map.total_contacts} / {contact_map.long_range_contacts} (density {contact_map.contact_density_pct:.2f}%)", "Calculated (Structure)"],
        ]
        story.append(_styled_table(cm_rows, col3_widths, styles))

    sscm_img = _make_ss_and_contact_figure(ss, contact_map)
    if sscm_img is not None:
        story.append(Spacer(1, 8))
        w = 6.6 * inch if (ss is not None and contact_map is not None) else 4.5 * inch
        story.append(Image(sscm_img, width=w, height=2.75 * inch))

    # 5. Methodology & Limitations
    sec5_head = Paragraph("5. Methodology, Software & Limitations", styles["h2"])
    sec5_body = Paragraph(
        "• <b>Sequence Descriptors:</b> Biopython ProtParam (average isotopic mass, reduced cysteines, Guruprasad et al. 1990 instability index, Kyte-Doolittle 1982 hydropathy).<br/>"
        "• <b>Pairwise Alignment:</b> Biopython PairwiseAligner / Gotoh affine-gap dynamic programming with NCBI BLOSUM/PAM matrices.<br/>"
        "• <b>Structure Parsing &amp; Analytics:</b> Model-1 primary-conformer parsing (PDB &amp; mmCIF), Richardson Lab MolProbity Top8000 6-class Ramachandran distributions (CC BY 4.0; Williams et al. 2018; Lovell et al. 2003), Kabsch-Sander (1983) electrostatic H-bond secondary structure, and Euclidean Cα/Cβ contact maps (Glycine Cα fallback in Cβ mode).<br/>"
        "• <b>Limitations:</b> Predicted structures (AlphaFold DB / ESMFold) are computational hypotheses rather than experimental coordinates. Sequence-based heuristics (e.g., instability index) do not guarantee in vitro stability. Full crystallographic refinement validation requires external tools such as MolProbity or RCSB validation reports.",
        styles["muted"],
    )
    story.append(KeepTogether([sec5_head, sec5_body]))

    doc.build(story, onFirstPage=_draw_page_footer, onLaterPages=_draw_page_footer)
    return buffer.getvalue()


def generate_text_report(
    parsed: ParsedStructure | None,
    seq_profile: SequenceProfile | None,
    alignment: AlignmentResult | None = None,
    rama: RamachandranResult | None = None,
    ss: SecondaryStructureResult | None = None,
    contact_map: ContactMapResult | None = None,
    sequence_name: str | None = None,
    sequence_source: str | None = None,
    sequence_mismatch_warning: str | None = None,
    sequence_unavailable_reason: str | None = None,
    report_chain: str | None = None,
) -> str:
    """Build a structured plain-text summary report preserving legacy .txt export."""
    lines = [
        f"{APP_NAME} v{APP_VERSION} — Scientific Analysis Summary",
        f"Generated (UTC): {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "PROVENANCE NOTICE: Computational descriptors are marked [CALCULATED], experimental file records [MEASURED], and AI predictions [AI-PREDICTED].",
        "",
    ]
    if parsed is not None:
        lines += [
            "STRUCTURE SUMMARY",
            f"File: {parsed.filename} ({parsed.file_format.upper()})",
            f"Origin: {ORIGIN_LABELS.get(parsed.structure_origin, parsed.structure_origin)}",
            f"PDB ID: {parsed.pdb_id or 'Unavailable'}",
            f"Method: {parsed.experimental_method or 'Unavailable'}",
            f"Resolution: {f'{parsed.resolution_angstrom:.2f} A' if parsed.resolution_angstrom is not None else 'Unavailable'}",
            f"Chains: {', '.join(parsed.chains)} ({parsed.polymer_residue_count} polymer residues, {parsed.polymer_atom_count} protein atoms in Model 1)",
            f"Waters: {parsed.water_count} | Hetero ligands: {parsed.hetero_ligand_count} ({', '.join(parsed.hetero_ligand_names) or 'none'})",
            f"Missing residues: {parsed.missing_residue_note}",
        ]
        if parsed.is_predicted and parsed.plddt_stats:
            lines.append(
                f"[AI-PREDICTED] Mean pLDDT: {parsed.plddt_stats['mean']:.2f} "
                f"(scale: {parsed.plddt_stats['raw_scale']})"
            )
        elif parsed.bfactor_stats:
            bs = parsed.bfactor_stats
            line = f"[MEASURED] Mean CA B-factor: {bs['mean']:.2f} A^2"
            if "all_atom_mean" in bs:
                line += (
                    f" | All-atom ({bs['all_atom_count']} atoms): mean {bs['all_atom_mean']:.2f} A^2 "
                    f"(range {bs['all_atom_min']:.2f}-{bs['all_atom_max']:.2f})"
                    f" | Occupancy: {bs['occupancy_min']:.2f}-{bs['occupancy_max']:.2f}"
                )
            lines.append(line)
        lines.append("")

    if seq_profile is not None:
        lines += [
            "SEQUENCE METRICS [CALCULATED]",
            f"Sequence Source: {_format_seq_source_label(sequence_source, report_chain)}",
            f"Sequence Label: {sequence_name or 'Active sequence'}",
        ]
        if sequence_mismatch_warning:
            lines.append(f"WARNING: {sequence_mismatch_warning}")
        lines += [
            f"Length: {seq_profile.length} aa",
            f"Molecular Weight: {seq_profile.molecular_weight_da:.2f} Da",
            f"Theoretical pI: {seq_profile.theoretical_pi:.2f}",
            f"Aromaticity: {seq_profile.aromaticity:.3f}",
            f"Instability Index: {seq_profile.instability_index:.2f} ({seq_profile.instability_class})",
            f"GRAVY: {seq_profile.gravy:.3f}",
            f"Extinction (280 nm, red/ox): {seq_profile.extinction_reduced} / {seq_profile.extinction_oxidized} M^-1 cm^-1",
            f"Sequence: {seq_profile.sequence}",
            "",
        ]
    else:
        reason = sequence_unavailable_reason or "No valid standard amino-acid sequence was analyzed."
        lines += [
            "SEQUENCE METRICS [UNAVAILABLE]",
            f"Reason: {reason}",
            "",
        ]

    if alignment is not None:
        lines += [
            "PAIRWISE ALIGNMENT [CALCULATED]",
            f"Seq A: {alignment.seq1_name} (Accession: {alignment.seq1_accession or 'None'}, {alignment.seq1_length} aa)",
            f"Seq B: {alignment.seq2_name} (Accession: {alignment.seq2_accession or 'None'}, {alignment.seq2_length} aa)",
            f"Mode: {alignment.mode} | Matrix: {alignment.matrix_name} | Open={alignment.open_gap_score} | Extend={alignment.extend_gap_score}",
            f"Score: {alignment.score:.1f} | Aligned Length: {alignment.alignment_length} cols | Gap Count: {alignment.gap_count}",
            f"Identity: {alignment.identity_pct:.2f}% | Similarity: {alignment.similarity_pct:.2f}% | Gaps: {alignment.gap_pct:.2f}%",
            "",
        ]
    else:
        lines += [
            "PAIRWISE ALIGNMENT [UNAVAILABLE]",
            "No pairwise sequence alignment was performed for the current structure/report session.",
            "",
        ]

    if rama is not None:
        lines += [
            "RAMACHANDRAN ANALYSIS [CALCULATED]",
            f"Evaluated: {rama.total_evaluated} | Favored: {rama.favored_pct:.1f}% | Allowed: {rama.allowed_pct:.1f}% | Outliers: {rama.outlier_pct:.1f}%",
            "",
        ]

    if ss is not None:
        lines += [
            "SECONDARY STRUCTURE [CALCULATED]",
            f"Method: {ss.method}",
            f"Helix: {ss.helix_pct:.1f}% | Sheet: {ss.sheet_pct:.1f}% | Coil/Turn: {ss.coil_pct:.1f}%",
            "",
        ]

    if contact_map is not None:
        lines += [
            "CONTACT MAP [CALCULATED]",
            f"Chain {contact_map.chain_id} ({contact_map.atom_mode}, <= {contact_map.threshold_angstrom:.1f} A, |i-j| >= {contact_map.min_seq_sep}): "
            f"{contact_map.total_contacts} contacts ({contact_map.long_range_contacts} long-range)",
            "",
        ]

    return "\n".join(lines)
