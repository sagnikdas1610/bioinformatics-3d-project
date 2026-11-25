"""Interactive Plotly charts styled cleanly for the BioFold Workbench dark theme."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from bio3d.config import AA_GROUPS, GROUP_COLORS, STANDARD_AA
from bio3d.sequence import SequenceProfile
from bio3d.structural_analysis import (
    ContactMapResult,
    RamachandranResult,
    SecondaryStructureResult,
    rama_background_grid,
)

_PLOT_LAYOUT = dict(
    template="plotly_dark",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="#0f172a",
    font=dict(color="#e2e8f0", size=12),
    margin=dict(l=52, r=24, t=42, b=46),
)


def plot_aa_composition(profile: SequenceProfile) -> go.Figure:
    aas = list(STANDARD_AA)
    counts = [profile.aa_counts[aa] for aa in aas]
    pcts = [profile.aa_percents[aa] for aa in aas]
    groups = [AA_GROUPS[aa] for aa in aas]
    colors = [GROUP_COLORS[g] for g in groups]

    fig = go.Figure(
        data=go.Bar(
            x=aas,
            y=counts,
            marker_color=colors,
            customdata=list(zip(pcts, groups)),
            hovertemplate="<b>%{x}</b> (%{customdata[1]})<br>Count: %{y}<br>Share: %{customdata[0]:.2f}%<extra></extra>",
        )
    )
    fig.update_layout(
        **_PLOT_LAYOUT,
        title="Amino-acid composition by physicochemical class",
        xaxis_title="Amino acid (1-letter code)",
        yaxis_title="Residue count",
        height=340,
    )
    return fig


def plot_hydropathy_profile(positions: list[int], scores: list[float], window: int) -> go.Figure:
    fig = go.Figure()
    fig.add_hline(y=0.0, line_dash="dot", line_color="#64748b")
    fig.add_hline(
        y=1.6,
        line_dash="dash",
        line_color="#f59e0b",
        annotation_text="Transmembrane heuristic (+1.6)",
        annotation_position="top right",
    )
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=scores,
            mode="lines",
            line=dict(color="#38bdf8", width=2.2),
            fill="tozeroy",
            fillcolor="rgba(56, 189, 248, 0.14)",
            hovertemplate="Residue %{x}<br>Mean hydropathy: %{y:.2f}<extra></extra>",
            name="Kyte-Doolittle",
        )
    )
    fig.update_layout(
        **_PLOT_LAYOUT,
        title=f"Kyte-Doolittle hydropathy profile (window = {window})",
        xaxis_title="Sequence position (window center)",
        yaxis_title="Hydropathy score",
        height=320,
        showlegend=False,
    )
    return fig


def plot_charge_curve(phs: list[float], charges: list[float], pi: float) -> go.Figure:
    fig = go.Figure()
    fig.add_hline(y=0.0, line_dash="dot", line_color="#64748b")
    fig.add_vline(
        x=pi,
        line_dash="dash",
        line_color="#2dd4bf",
        annotation_text=f"pI = {pi:.2f}",
        annotation_position="top right",
    )
    fig.add_trace(
        go.Scatter(
            x=phs,
            y=charges,
            mode="lines+markers",
            line=dict(color="#a78bfa", width=2.3),
            marker=dict(size=5),
            hovertemplate="pH %{x:.1f}<br>Net charge: %{y:+.2f}<extra></extra>",
        )
    )
    fig.update_layout(
        **_PLOT_LAYOUT,
        title="Theoretical net charge vs. pH",
        xaxis_title="pH",
        yaxis_title="Estimated net charge",
        height=320,
        showlegend=False,
    )
    return fig


def plot_ramachandran(
    rama: RamachandranResult,
    bg_class: str = "General",
) -> go.Figure:
    phis_grid, psis_grid, z_grid = rama_background_grid(bg_class)  # type: ignore[arg-type]
    fig = go.Figure()
    fig.add_trace(
        go.Contour(
            x=phis_grid,
            y=psis_grid,
            z=z_grid,
            contours=dict(start=0.5, end=2.5, size=1.0, coloring="fill"),
            colorscale=[
                [0.0, "#0f172a"],
                [0.33, "#0f172a"],
                [0.5, "#1e3a5f"],
                [1.0, "#1d4ed8"],
            ],
            opacity=0.48,
            showscale=False,
            hoverinfo="skip",
            line=dict(width=0.6, color="rgba(148, 163, 184, 0.35)"),
        )
    )

    region_styles = [
        ("Favored", "#38bdf8", "circle", 7),
        ("Allowed", "#fbbf24", "diamond", 8),
        ("Outlier", "#f87171", "x", 10),
    ]
    for region_name, color, symbol, size in region_styles:
        pts = [r for r in rama.residues if r.region == region_name]
        if not pts:
            continue
        fig.add_trace(
            go.Scatter(
                x=[p.phi for p in pts],
                y=[p.psi for p in pts],
                mode="markers",
                name=f"{region_name} ({len(pts)})",
                marker=dict(color=color, symbol=symbol, size=size, line=dict(width=0.8, color="#0b1120")),
                customdata=[[p.chain, p.resname, p.resseq, p.rama_class, p.prob] for p in pts],
                hovertemplate=(
                    "<b>%{customdata[0]}:%{customdata[1]}%{customdata[2]}</b> (%{customdata[3]})<br>"
                    "φ = %{x:.1f}°, ψ = %{y:.1f}° · Top8000 p = %{customdata[4]:.4f}<extra>" + region_name + "</extra>"
                ),
            )
        )

    fig.add_hline(y=0, line_width=0.6, line_color="rgba(148,163,184,0.3)")
    fig.add_vline(x=0, line_width=0.6, line_color="rgba(148,163,184,0.3)")
    fig.update_layout(
        **_PLOT_LAYOUT,
        title=f"Ramachandran plot (MolProbity Top8000 contour: {bg_class})",
        xaxis=dict(title="φ (deg)", range=[-180, 180], dtick=60),
        yaxis=dict(title="ψ (deg)", range=[-180, 180], dtick=60),
        height=500,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1.0),
    )
    return fig


def plot_secondary_structure(ss: SecondaryStructureResult) -> go.Figure:
    code_num = {"H": 2, "G": 2, "I": 2, "E": 1, "B": 1, "T": 0, "C": 0}
    x_vals = [f"{r.chain}:{r.resseq}" for r in ss.residues]
    y_vals = [code_num[r.dssp_code] for r in ss.residues]
    colors = [
        "#ec4899" if r.state3 == "Helix" else ("#f59e0b" if r.state3 == "Sheet" else "#38bdf8")
        for r in ss.residues
    ]
    custom = [[r.resname, r.dssp_code, r.state3] for r in ss.residues]

    fig = go.Figure(
        data=go.Bar(
            x=x_vals,
            y=[1] * len(x_vals),
            marker_color=colors,
            customdata=custom,
            hovertemplate="<b>%{x}</b> (%{customdata[0]})<br>DSSP: %{customdata[1]} (%{customdata[2]})<extra></extra>",
        )
    )
    fig.update_layout(
        **_PLOT_LAYOUT,
        title="Per-residue secondary structure track (pink = Helix, amber = Sheet, cyan = Coil/Turn)",
        xaxis_title="Residue",
        yaxis=dict(visible=False),
        bargap=0.05,
        height=200,
    )
    return fig


def plot_contact_map(cmap: ContactMapResult) -> go.Figure:
    z = np.where(cmap.contact_matrix, cmap.distance_matrix, np.nan)
    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=cmap.residue_labels,
            y=cmap.residue_labels,
            colorscale="Teal_r",
            zmin=2.5,
            zmax=cmap.threshold_angstrom,
            colorbar=dict(title="Dist (Å)"),
            hovertemplate="Res 1: %{y}<br>Res 2: %{x}<br>Distance: %{z:.2f} Å<extra></extra>",
        )
    )
    fig.update_layout(
        **_PLOT_LAYOUT,
        title=(
            f"Chain {cmap.chain_id} contact map ({cmap.atom_mode} mode, "
            f"cutoff ≤ {cmap.threshold_angstrom:.1f} Å, |i-j| ≥ {cmap.min_seq_sep})"
        ),
        xaxis_title="Residue",
        yaxis=dict(title="Residue", scaleanchor="x", scaleratio=1),
        height=520,
    )
    return fig


def plot_plddt_profile(plddt_stats: dict[str, object]) -> go.Figure:
    rows = plddt_stats["per_residue"]  # type: ignore[index]
    x = [f"{r['chain']}:{r['resseq']}" for r in rows]
    y = [float(r["value"]) for r in rows]
    colors = [
        "#0053d6" if v >= 90 else ("#65cbf3" if v >= 70 else ("#ffdb13" if v >= 50 else "#ff7d45"))
        for v in y
    ]
    fig = go.Figure(
        data=go.Bar(
            x=x,
            y=y,
            marker_color=colors,
            hovertemplate="<b>%{x}</b><br>pLDDT: %{y:.1f}<extra></extra>",
        )
    )
    fig.add_hline(y=90, line_dash="dot", line_color="#60a5fa")
    fig.add_hline(y=70, line_dash="dot", line_color="#fbbf24")
    fig.update_layout(
        **_PLOT_LAYOUT,
        title="Per-residue pLDDT confidence profile (0–100 scale)",
        xaxis_title="Residue",
        yaxis=dict(title="pLDDT confidence", range=[0, 100]),
        height=290,
    )
    return fig


def plot_bfactor_profile(bfactor_stats: dict[str, object]) -> go.Figure:
    rows = bfactor_stats["per_residue"]  # type: ignore[index]
    x = [f"{r['chain']}:{r['resseq']}" for r in rows]
    y = [float(r["value"]) for r in rows]
    fig = go.Figure(
        data=go.Scatter(
            x=x,
            y=y,
            mode="lines",
            line=dict(color="#f59e0b", width=2.0),
            fill="tozeroy",
            fillcolor="rgba(245, 158, 11, 0.14)",
            hovertemplate="<b>%{x}</b><br>B-factor: %{y:.2f} Å²<extra></extra>",
        )
    )
    fig.update_layout(
        **_PLOT_LAYOUT,
        title="Experimental Cα crystallographic B-factor profile (Å²)",
        xaxis_title="Residue",
        yaxis_title="B-factor (Å²)",
        height=290,
    )
    return fig
