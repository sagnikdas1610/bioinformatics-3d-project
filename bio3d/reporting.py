"""Scientific PDF and plain-text report generation using ReportLab and Matplotlib (Agg)."""

from __future__ import annotations

from datetime import datetime, timezone
from io import BytesIO
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "biofold_mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import LETTER  # noqa: E402
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet  # noqa: E402
from reportlab.lib.units import inch  # noqa: E402
from reportlab.platypus import (  # noqa: E402
    HRFlowable,
    Image,
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

    fig = plt.figure(figsize=(5.4, 3.6))
    ax = fig.add_subplot(111, projection="3d")
    palette = ["#0284c7", "#7c3aed", "#059669", "#d97706", "#dc2626", "#0d9488"]

    for idx, ch in enumerate(parsed.chains):
        pts = [r.atoms["CA"].coord for r in parsed.polymer_residues[ch] if "CA" in r.atoms]
        if len(pts) < 2:
            continue
        arr = np.vstack(pts)
        ax.plot(arr[:, 0], arr[:, 1], arr[:, 2], linewidth=1.8, color=palette[idx % len(palette)], label=f"Chain {ch}")
        ax.scatter(arr[0, 0], arr[0, 1], arr[0, 2], color="#16a34a", s=24)

    ax.set_title(f"Cα Backbone Trace — {parsed.pdb_id or parsed.filename}", fontsize=10, fontweight="bold")
    ax.set_xlabel("X (Å)", fontsize=8)
    ax.set_ylabel("Y (Å)", fontsize=8)
    ax.set_zlabel("Z (Å)", fontsize=8)
    ax.tick_params(labelsize=7)
    if len(parsed.chains) <= 6:
        ax.legend(fontsize=7, loc="upper right")
    return _fig_to_png_bytes(fig)


def _make_sequence_figure(
    profile: SequenceProfile,
    hydropathy: tuple[list[int], list[float]] | None,
) -> BytesIO:
    """Render amino-acid composition and Kyte-Doolittle hydropathy figures."""
    ncols = 2 if hydropathy and len(hydropathy[0]) > 1 else 1
    fig, axes = plt.subplots(1, ncols, figsize=(7.0 if ncols == 2 else 5.2, 2.7))
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
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
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
    fig, axes = plt.subplots(1, ncols, figsize=(7.0 if ncols == 2 else 4.8, 3.0))

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


def _styled_table(rows: list[list[str]], col_widths: list[float]) -> Table:
    table = Table(rows, colWidths=col_widths)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0f172a")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8fafc")]),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cbd5e1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def generate_pdf_report(
    parsed: ParsedStructure | None,
    seq_profile: SequenceProfile | None,
    hydropathy: tuple[list[int], list[float]] | None = None,
    alignment: AlignmentResult | None = None,
    rama: RamachandranResult | None = None,
    ss: SecondaryStructureResult | None = None,
    contact_map: ContactMapResult | None = None,
    report_title: str | None = None,
) -> bytes:
    """Build a publication-style scientific PDF report distinguishing measured/calculated/predicted data."""
    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=LETTER,
        leftMargin=0.65 * inch,
        rightMargin=0.65 * inch,
        topMargin=0.65 * inch,
        bottomMargin=0.65 * inch,
        title=f"{APP_NAME} Scientific Report",
    )

    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("BFH1", parent=styles["Heading1"], fontSize=16, textColor=colors.HexColor("#0f172a"), spaceAfter=4)
    h2 = ParagraphStyle("BFH2", parent=styles["Heading2"], fontSize=11.5, textColor=colors.HexColor("#0369a1"), spaceBefore=10, spaceAfter=4)
    body = ParagraphStyle("BFBody", parent=styles["Normal"], fontSize=9, leading=12.5, textColor=colors.HexColor("#1e293b"))
    muted = ParagraphStyle("BFMuted", parent=styles["Normal"], fontSize=8, leading=11, textColor=colors.HexColor("#475569"))

    story: list[object] = []
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    story.append(Paragraph(report_title or f"{APP_NAME} — Structural & Sequence Analysis Report", h1))
    story.append(Paragraph(f"Generated: {timestamp} | Platform v{APP_VERSION} | Provenance tags distinguish [MEASURED], [CALCULATED], [AI-PREDICTED], and [UNAVAILABLE] fields.", muted))
    story.append(Spacer(1, 6))
    story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#cbd5e1")))
    story.append(Spacer(1, 8))

    # 1. Structure Identification & Quality
    story.append(Paragraph("1. Structure Identification & Quality Summary", h2))
    if parsed is not None:
        origin_label = ORIGIN_LABELS.get(parsed.structure_origin, parsed.structure_origin)
        prov_tag = "[AI-PREDICTED]" if parsed.is_predicted else "[MEASURED / FILE RECORD]"
        struct_rows = [
            ["Field", "Value", "Provenance"],
            ["Source File", parsed.filename, "[FILE]"],
            ["Origin Classification", origin_label, prov_tag],
            ["PDB / Entry ID", parsed.pdb_id or "Unavailable in file header", "[MEASURED]" if parsed.pdb_id else "[UNAVAILABLE]"],
            ["Experimental Method", parsed.experimental_method or ("Computational prediction" if parsed.is_predicted else "Unavailable"), "[MEASURED]" if parsed.experimental_method else "[UNAVAILABLE]"],
            ["Resolution", f"{parsed.resolution_angstrom:.2f} Å" if parsed.resolution_angstrom is not None else "Unavailable (not reported / non-diffraction)", "[MEASURED]" if parsed.resolution_angstrom is not None else "[UNAVAILABLE]"],
            ["Polymer Chains / Residues", f"{len(parsed.chains)} ({', '.join(parsed.chains)}) / {parsed.polymer_residue_count} aa", "[CALCULATED]"],
            ["Polymer Atoms (Model 1)", f"{parsed.polymer_atom_count:,} (models in file: {parsed.model_count})", "[CALCULATED]"],
            ["Hetero Ligands / Waters", f"{parsed.hetero_ligand_count} ({', '.join(parsed.hetero_ligand_names) or 'none'}) / {parsed.water_count} waters", "[CALCULATED]"],
            ["Missing Residues Check", parsed.missing_residue_note, "[FILE / CALCULATED]"],
        ]
        if parsed.is_predicted and parsed.plddt_stats:
            ps = parsed.plddt_stats
            struct_rows.append([
                "pLDDT Confidence (0–100)",
                f"Mean {ps['mean']:.1f} | Median {ps['median']:.1f} | ≥90: {ps['pct_very_high']:.1f}% | <50: {ps['pct_very_low']:.1f}% (raw scale: {ps['raw_scale']})",
                "[AI-PREDICTED]",
            ])
        elif parsed.bfactor_stats:
            bs = parsed.bfactor_stats
            struct_rows.append([
                "Experimental B-factor (Å²)",
                f"Mean {bs['mean']:.2f} | Median {bs['median']:.2f} | Range {bs['min']:.2f}–{bs['max']:.2f} Å²",
                "[MEASURED]",
            ])
        story.append(_styled_table(struct_rows, [1.65 * inch, 4.15 * inch, 1.35 * inch]))
        backbone_png = _make_backbone_3d_figure(parsed)
        if backbone_png:
            story.append(Spacer(1, 6))
            story.append(Image(backbone_png, width=4.6 * inch, height=3.05 * inch))
    else:
        story.append(Paragraph("[UNAVAILABLE] No 3D structure was loaded during this session.", body))

    # 2. Primary Sequence Statistics
    story.append(Paragraph("2. Primary Sequence & Molecular Properties", h2))
    if seq_profile is not None:
        seq_rows = [
            ["Descriptor", "Value", "Provenance"],
            ["Length", f"{seq_profile.length} amino acids", "[CALCULATED]"],
            ["Molecular Weight (avg mass, reduced Cys)", f"{seq_profile.molecular_weight_da:,.2f} Da ({seq_profile.molecular_weight_kda:.2f} kDa)", "[CALCULATED]"],
            ["Theoretical Isoelectric Point (pI)", f"{seq_profile.theoretical_pi:.2f}", "[CALCULATED]"],
            ["GRAVY Hydropathy (Kyte-Doolittle)", f"{seq_profile.gravy:.3f}", "[CALCULATED]"],
            ["Instability Index (Guruprasad heuristic)", f"{seq_profile.instability_index:.2f} ({seq_profile.instability_class})", "[CALCULATED]"],
            ["Aromaticity (Phe+Trp+Tyr fraction)", f"{seq_profile.aromaticity:.3f}", "[CALCULATED]"],
            ["Molar Extinction (280 nm, red / ox)", f"{seq_profile.extinction_reduced:,} / {seq_profile.extinction_oxidized:,} M⁻¹cm⁻¹", "[CALCULATED]"],
            ["Estimated Net Charge at pH 7.0", f"{seq_profile.net_charge_ph7:+.2f}", "[CALCULATED]"],
        ]
        story.append(_styled_table(seq_rows, [2.3 * inch, 3.5 * inch, 1.35 * inch]))
        story.append(Spacer(1, 6))
        seq_fig = _make_sequence_figure(seq_profile, hydropathy)
        story.append(Image(seq_fig, width=6.7 * inch, height=2.55 * inch))
    else:
        story.append(Paragraph("[UNAVAILABLE] No valid standard amino-acid sequence was analyzed.", body))

    # 3. Pairwise Alignment
    story.append(Paragraph("3. Pairwise Sequence Alignment", h2))
    if alignment is not None:
        aln_rows = [
            ["Metric", "Value", "Provenance"],
            ["Mode & Scoring Matrix", f"{alignment.mode.capitalize()} | {alignment.matrix_name} (gap open={alignment.open_gap_score}, ext={alignment.extend_gap_score})", "[CALCULATED]"],
            ["Alignment Score", f"{alignment.score:.1f} (Engine: {alignment.engine})", "[CALCULATED]"],
            ["Sequence Identity (aligned columns)", f"{alignment.identity_pct:.2f}% ({alignment.identical_count}/{alignment.alignment_length})", "[CALCULATED]"],
            ["Sequence Similarity (score > 0)", f"{alignment.similarity_pct:.2f}% ({alignment.similar_count}/{alignment.alignment_length})", "[CALCULATED]"],
            ["Gaps", f"{alignment.gap_pct:.2f}% ({alignment.gap_count} positions)", "[CALCULATED]"],
        ]
        story.append(_styled_table(aln_rows, [2.2 * inch, 3.6 * inch, 1.35 * inch]))
        story.append(Spacer(1, 4))
        story.append(Paragraph(alignment.interpretation, muted))
    else:
        story.append(Paragraph("[UNAVAILABLE] Pairwise sequence alignment was not executed in this session.", body))

    # 4. Advanced Structural Analysis
    story.append(Paragraph("4. Structural Stereochemistry, Secondary Structure & Contacts", h2))
    if rama is not None:
        rama_rows = [
            ["Ramachandran Metric (MolProbity Top8000)", "Value", "Provenance"],
            ["Reference Distribution", rama.reference_source, "[REFERENCE]"],
            ["Evaluated Internal Residues", f"{rama.total_evaluated} (excluded termini/breaks: {rama.excluded_count})", "[CALCULATED]"],
            ["Favored / Allowed / Outliers", f"{rama.favored_pct:.1f}% ({rama.favored_count}) / {rama.allowed_pct:.1f}% ({rama.allowed_count}) / {rama.outlier_pct:.1f}% ({rama.outlier_count})", "[CALCULATED]"],
            ["Stereochemical Assessment", rama.interpretation, "[CALCULATED]"],
        ]
        story.append(_styled_table(rama_rows, [2.2 * inch, 3.6 * inch, 1.35 * inch]))
        story.append(Spacer(1, 6))
        story.append(Image(_make_ramachandran_figure(rama), width=4.0 * inch, height=3.6 * inch))
        if rama.outliers:
            out_rows = [["Chain", "Residue", "Class", "φ (deg)", "ψ (deg)"]]
            for o in rama.outliers[:12]:
                out_rows.append([o.chain, f"{o.resname} {o.resseq}", o.rama_class, f"{o.phi:.1f}", f"{o.psi:.1f}"])
            story.append(Spacer(1, 4))
            story.append(_styled_table(out_rows, [0.9 * inch, 1.6 * inch, 1.6 * inch, 1.5 * inch, 1.55 * inch]))
    else:
        story.append(Paragraph("[UNAVAILABLE] Ramachandran analysis was not computed.", body))

    if ss is not None:
        story.append(Spacer(1, 6))
        ss_rows = [
            ["Secondary Structure Metric", "Value", "Provenance"],
            ["Method", ss.method, "[CALCULATED]"],
            ["Helix (H/G/I) / Sheet (E/B) / Coil & Turn", f"{ss.helix_pct:.1f}% ({ss.helix_count}) / {ss.sheet_pct:.1f}% ({ss.sheet_count}) / {ss.coil_pct:.1f}% ({ss.coil_count})", "[CALCULATED]"],
        ]
        if ss.header_available and ss.header_helix_pct is not None:
            ss_rows.append([
                "File Header HELIX/SHEET Annotation",
                f"Helix {ss.header_helix_pct:.1f}% / Sheet {ss.header_sheet_pct:.1f}% / Coil {ss.header_coil_pct:.1f}%",
                "[MEASURED / FILE RECORD]",
            ])
        story.append(_styled_table(ss_rows, [2.2 * inch, 3.6 * inch, 1.35 * inch]))

    if contact_map is not None:
        story.append(Spacer(1, 6))
        gly_note = (
            f" (Glycine Cα fallback applied to {contact_map.glycine_fallback_count} Gly residues)"
            if contact_map.atom_mode == "CB"
            else ""
        )
        cm_rows = [
            ["Contact Map Metric", "Value", "Provenance"],
            ["Chain & Mode", f"Chain {contact_map.chain_id} | {contact_map.atom_mode} atoms{gly_note}", "[CALCULATED]"],
            ["Distance Cutoff & Seq Separation", f"≤ {contact_map.threshold_angstrom:.1f} Å | |i - j| ≥ {contact_map.min_seq_sep}", "[CALCULATED]"],
            ["Total Contacts / Long-Range (|i-j| ≥ 12)", f"{contact_map.total_contacts} / {contact_map.long_range_contacts} (density {contact_map.contact_density_pct:.2f}%)", "[CALCULATED]"],
        ]
        story.append(_styled_table(cm_rows, [2.2 * inch, 3.6 * inch, 1.35 * inch]))

    sscm_img = _make_ss_and_contact_figure(ss, contact_map)
    if sscm_img is not None:
        story.append(Spacer(1, 6))
        w = 6.7 * inch if (ss is not None and contact_map is not None) else 4.6 * inch
        story.append(Image(sscm_img, width=w, height=2.85 * inch))

    # 5. Methodology & Limitations
    story.append(Paragraph("5. Methodology, Software & Limitations", h2))
    story.append(
        Paragraph(
            "• <b>Sequence Descriptors:</b> Biopython ProtParam (average isotopic mass, reduced cysteines, Guruprasad et al. 1990 instability index, Kyte-Doolittle 1982 hydropathy).<br/>"
            "• <b>Pairwise Alignment:</b> Biopython PairwiseAligner / Gotoh affine-gap dynamic programming with NCBI BLOSUM/PAM matrices.<br/>"
            "• <b>Structure Parsing & Analytics:</b> Model-1 primary-conformer parsing (PDB & mmCIF), Richardson Lab MolProbity Top8000 6-class Ramachandran distributions (CC BY 4.0; Williams et al. 2018; Lovell et al. 2003), Kabsch-Sander (1983) electrostatic H-bond secondary structure, and Euclidean Cα/Cβ contact maps (Glycine Cα fallback in Cβ mode).<br/>"
            "• <b>Limitations:</b> Predicted structures (AlphaFold DB / ESMFold) are computational hypotheses rather than experimental coordinates. Sequence-based heuristics (e.g., instability index) do not guarantee in vitro stability. Full crystallographic refinement validation requires external tools such as MolProbity or RCSB validation reports.",
            muted,
        )
    )

    doc.build(story)
    return buffer.getvalue()


def generate_text_report(
    parsed: ParsedStructure | None,
    seq_profile: SequenceProfile | None,
    alignment: AlignmentResult | None = None,
    rama: RamachandranResult | None = None,
    ss: SecondaryStructureResult | None = None,
    contact_map: ContactMapResult | None = None,
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
            f"Chains: {', '.join(parsed.chains)} ({parsed.polymer_residue_count} polymer residues, {parsed.polymer_atom_count} atoms in Model 1)",
            f"Hetero ligands: {parsed.hetero_ligand_count} | Waters: {parsed.water_count}",
            f"Missing residues: {parsed.missing_residue_note}",
        ]
        if parsed.is_predicted and parsed.plddt_stats:
            lines.append(
                f"[AI-PREDICTED] Mean pLDDT: {parsed.plddt_stats['mean']:.2f} "
                f"(scale: {parsed.plddt_stats['raw_scale']})"
            )
        elif parsed.bfactor_stats:
            lines.append(f"[MEASURED] Mean B-factor: {parsed.bfactor_stats['mean']:.2f} A^2")
        lines.append("")

    if seq_profile is not None:
        lines += [
            "SEQUENCE METRICS [CALCULATED]",
            f"Length: {seq_profile.length} aa",
            f"Molecular Weight: {seq_profile.molecular_weight_da:.2f} Da",
            f"Theoretical pI: {seq_profile.theoretical_pi:.2f}",
            f"Aromaticity: {seq_profile.aromaticity:.3f}",
            f"Instability Index: {seq_profile.instability_index:.2f} ({seq_profile.instability_class})",
            f"GRAVY: {seq_profile.gravy:.3f}",
            f"Sequence: {seq_profile.sequence}",
            "",
        ]

    if alignment is not None:
        lines += [
            "PAIRWISE ALIGNMENT [CALCULATED]",
            f"Mode: {alignment.mode} | Matrix: {alignment.matrix_name} | Score: {alignment.score:.1f}",
            f"Identity: {alignment.identity_pct:.2f}% | Similarity: {alignment.similarity_pct:.2f}% | Gaps: {alignment.gap_pct:.2f}%",
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
