from __future__ import annotations

import base64
from collections import Counter
from datetime import datetime, timezone
from html import escape
import json
from pathlib import Path
import re

import pandas as pd
import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import streamlit.components.v1 as components

APP_NAME = "BioFold Workbench"
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
SAMPLE_SEQUENCE = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)
MAX_PDB_BYTES = 20 * 1024 * 1024
DEMO_PDB_PATH = Path(__file__).resolve().parent / "data" / "example_ubiquitin_fragment.pdb"

st.set_page_config(
    page_title=f"{APP_NAME} | 3D Protein Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Glassmorphism + gradient theme. All core analysis functionality remains in Python/Streamlit.
st.markdown(
    r"""
    <style>
      :root {
        --bg-0: #070b18;
        --bg-1: #10172b;
        --text-main: #edf4ff;
        --text-muted: #9eacc8;
        --line: rgba(174, 195, 255, .17);
        --glass: rgba(20, 29, 54, .68);
        --aqua: #67e8f9;
        --violet: #a78bfa;
      }
      .stApp {
        color: var(--text-main);
        background:
          radial-gradient(ellipse at 12% 7%, rgba(54, 95, 205, .27), transparent 35%),
          radial-gradient(ellipse at 88% 16%, rgba(139, 92, 246, .22), transparent 33%),
          radial-gradient(ellipse at 62% 75%, rgba(14, 165, 233, .10), transparent 38%),
          linear-gradient(145deg, var(--bg-0) 0%, var(--bg-1) 55%, #0a1122 100%);
        background-attachment: fixed;
      }
      [data-testid="stHeader"] { background: rgba(7, 11, 24, .55); }
      .main .block-container { max-width: 1680px; padding-top: 1.25rem; padding-bottom: 3rem; }
      h1, h2, h3, h4, p, label, [data-testid="stMarkdownContainer"] { color: var(--text-main); }
      h1, h2, h3 { letter-spacing: -.025em; }
      .hero-glass {
        position: relative; overflow: hidden; padding: 1.7rem 1.9rem 1.55rem;
        border: 1px solid rgba(188, 202, 255, .24); border-radius: 26px;
        background: linear-gradient(120deg, rgba(34, 49, 91, .66), rgba(28, 26, 62, .44) 53%, rgba(13, 86, 112, .28));
        box-shadow: 0 18px 60px rgba(0, 0, 0, .22), inset 0 1px 0 rgba(255,255,255,.07);
        backdrop-filter: blur(24px); -webkit-backdrop-filter: blur(24px); margin-bottom: 1.15rem;
      }
      .hero-glass:after { content: ""; position:absolute; width:280px; height:280px; border-radius:50%;
        right:-85px; top:-155px; background:rgba(103,232,249,.13); filter:blur(8px); pointer-events:none; }
      .eyebrow { color: #9cecf7; font-size: .72rem; font-weight: 800; letter-spacing: .18em; text-transform: uppercase; }
      .hero-title { margin: .4rem 0 .5rem; font-weight: 800; font-size: clamp(2rem, 3vw, 3.1rem); line-height: 1.08;
        color: #f7f9ff; letter-spacing: -.045em; }
      .gradient-text { background: linear-gradient(90deg, #a5f3fc, #c4b5fd 55%, #f0abfc); -webkit-background-clip:text; background-clip:text; color:transparent; }
      .hero-subtitle { color: #c1cce5; font-size: 1rem; line-height: 1.65; max-width: 900px; margin:0; }
      .status-row { display:flex; align-items:center; gap:.55rem; flex-wrap:wrap; margin-top:1rem; }
      .status-pill { display:inline-flex; align-items:center; gap:.45rem; padding:.38rem .72rem; border-radius:999px;
        border:1px solid rgba(192,205,255,.19); background:rgba(8,14,32,.36); color:#d7e3ff; font-size:.72rem; font-weight:700; letter-spacing:.035em; }
      .status-dot { width:7px; height:7px; border-radius:50%; display:inline-block; background:#67e8f9; box-shadow:0 0 12px rgba(103,232,249,.8); }
      .section-kicker { color:#9de8f2; font-size:.72rem; letter-spacing:.15em; font-weight:800; text-transform:uppercase; margin:.3rem 0 .2rem; }
      .section-note { color:#9eacc8; font-size:.84rem; margin-bottom:.65rem; }
      .soft-rule { height:1px; width:100%; margin:1.05rem 0; background:linear-gradient(90deg, transparent, rgba(176,194,255,.28), transparent); }
      div[data-testid="stMetric"] { border:1px solid var(--line); background:linear-gradient(145deg, rgba(31,43,74,.69), rgba(15,22,43,.58));
        padding:.85rem 1rem; border-radius:18px; box-shadow:0 12px 35px rgba(0,0,0,.13); backdrop-filter:blur(15px); }
      div[data-testid="stMetricLabel"] { color:#aab8d5 !important; }
      div[data-testid="stMetricValue"] { color:#f2f7ff !important; }
      div[data-testid="stMetricDelta"] { color:#a5f3fc !important; }
      [data-testid="stSidebar"] { border-right:1px solid rgba(178,195,255,.14); background:rgba(7,12,27,.72); }
      [data-testid="stSidebar"] > div:first-child { background:linear-gradient(180deg, rgba(19,28,54,.88), rgba(8,13,29,.92)); }
      [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p { color:#bdc9e2; }
      [data-testid="stTextArea"] textarea, [data-testid="stTextInput"] input {
        background:rgba(7,13,29,.72) !important; color:#eff6ff !important; border:1px solid rgba(160,180,244,.25) !important; border-radius:14px !important;
      }
      div[data-baseweb="select"] > div { background:rgba(10,17,37,.76) !important; border:1px solid rgba(160,180,244,.25) !important; border-radius:13px !important; color:#eff6ff !important; }
      [data-testid="stFileUploader"] section { background:rgba(10,17,37,.56); border:1px dashed rgba(148,171,245,.4); border-radius:16px; }
      [data-testid="stFileUploader"] small, [data-testid="stFileUploader"] span { color:#b7c5e2 !important; }
      .stButton > button, .stDownloadButton > button { border-radius:12px; border:1px solid rgba(160,180,244,.28); color:#eff6ff;
        background:linear-gradient(115deg, rgba(59,130,246,.25), rgba(139,92,246,.22)); transition:all .18s ease; }
      .stButton > button:hover, .stDownloadButton > button:hover { border-color:rgba(103,232,249,.7); color:white; transform:translateY(-1px); }
      div[data-testid="stTabs"] button { color:#aebbd7; }
      div[data-testid="stTabs"] button[aria-selected="true"] { color:#a5f3fc; border-bottom-color:#67e8f9; }
      div[data-testid="stAlert"] { border-radius:15px; border:1px solid rgba(160,180,244,.17); }
      .tiny-muted { color:#9eacc8; font-size:.79rem; line-height:1.55; }
      .viewer-instructions { margin:.3rem .15rem .1rem; color:#9eacc8; font-size:.78rem; line-height:1.6; }
      footer { visibility:hidden; }
      @media (max-width: 900px) {
        .hero-glass { padding:1.2rem; border-radius:20px; }
        .main .block-container { padding-left:1rem; padding-right:1rem; }
      }
    </style>
    """,
    unsafe_allow_html=True,
)


def clean_sequence(raw_text: str) -> str:
    """Remove FASTA headers and whitespace, and normalize to uppercase."""
    pieces: list[str] = []
    for line in raw_text.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        pieces.append(re.sub(r"\s+", "", line))
    return "".join(pieces).upper()


def sequence_metrics(sequence: str) -> dict[str, float | int]:
    analysis = ProteinAnalysis(sequence)
    return {
        "Length (residues)": len(sequence),
        "Molecular weight (Da)": analysis.molecular_weight(),
        "Theoretical pI": analysis.isoelectric_point(),
        "Aromaticity": analysis.aromaticity(),
        "Instability index": analysis.instability_index(),
        "GRAVY hydropathy": analysis.gravy(),
    }


def pdb_summary(pdb_text: str) -> dict[str, object]:
    atom_lines = [line for line in pdb_text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    residues: set[tuple[str, str, str]] = set()
    chains: set[str] = set()
    b_factors: list[float] = []
    for line in atom_lines:
        if len(line) >= 27:
            chain = line[21].strip() or "(blank)"
            residues.add((chain, line[22:26].strip(), line[26].strip()))
            chains.add(chain)
        if len(line) >= 66:
            try:
                b_factors.append(float(line[60:66].strip()))
            except ValueError:
                pass
    return {
        "Atom records": len(atom_lines),
        "Observed residues": len(residues),
        "Chains": sorted(chains),
        "Mean B-factor field": sum(b_factors) / len(b_factors) if b_factors else None,
        "File size (KiB)": len(pdb_text.encode("utf-8")) / 1024,
    }


def make_viewer_html(pdb_text: str, render_style: str, color_mode: str, source_label: str) -> str:
    """Build a self-contained viewer iframe; the molecular library loads from its public CDN."""
    encoded_pdb = base64.b64encode(pdb_text.encode("utf-8")).decode("ascii")
    style_map = {"Cartoon": "cartoon", "Stick": "stick", "Sphere": "sphere"}
    color_map = {"Rainbow": "rainbow", "By chain": "chain", "B-factor field": "bfactor"}
    mode_json = json.dumps(style_map[render_style])
    color_json = json.dumps(color_map[color_mode])
    source_html = escape(source_label)
    pdb_json = json.dumps(encoded_pdb)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  *{{box-sizing:border-box}} html,body{{width:100%;height:100%;margin:0;padding:0;background:transparent;font-family:Inter,ui-sans-serif,system-ui,-apple-system,Segoe UI,sans-serif;color:#eef4ff;overflow:hidden}}
  .viewer-card{{height:100%;display:flex;flex-direction:column;overflow:hidden;border-radius:24px;border:1px solid rgba(187,202,255,.24);
    background:linear-gradient(145deg,rgba(22,32,60,.97),rgba(10,17,36,.97) 62%,rgba(15,33,57,.97));
    box-shadow:0 20px 55px rgba(0,0,0,.27),inset 0 1px 0 rgba(255,255,255,.065)}}
  .viewer-head{{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:17px 19px;border-bottom:1px solid rgba(172,190,245,.13);background:linear-gradient(100deg,rgba(64,89,157,.15),rgba(45,212,191,.055))}}
  .kicker{{font-size:10px;font-weight:800;letter-spacing:.17em;text-transform:uppercase;color:#8ee8f5;margin-bottom:4px}}
  .title{{font-size:17px;font-weight:750;letter-spacing:-.025em;color:#f3f6ff}}
  .source{{max-width:44%;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:10px;font-weight:750;letter-spacing:.07em;text-transform:uppercase;
    padding:8px 10px;border:1px solid rgba(153,176,255,.2);border-radius:999px;background:rgba(5,11,28,.45);color:#c8d9fb}}
  #viewer{{position:relative;flex:1;min-height:360px;width:100%;background:radial-gradient(ellipse at 50% 42%,rgba(57,91,156,.2),transparent 42%),linear-gradient(150deg,#0c142a,#080f20)}}
  .viewer-foot{{display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap;padding:11px 17px;border-top:1px solid rgba(172,190,245,.13);background:rgba(6,12,27,.47);color:#93a6cb;font-size:11px}}
  .foot-chip{{display:inline-flex;align-items:center;gap:6px}} .dot{{width:6px;height:6px;border-radius:50%;background:#5eead4;box-shadow:0 0 9px #5eead488;display:inline-block}}
  .loadmsg{{padding:24px;font-size:13px;color:#c7d4f1}}
</style></head><body>
<div class="viewer-card">
  <div class="viewer-head"><div><div class="kicker">Molecular workspace / 01</div><div class="title">Interactive 3D structure</div></div><div class="source" title="{source_html}">{source_html}</div></div>
  <div id="viewer"><div class="loadmsg">Preparing molecular scene…</div></div>
  <div class="viewer-foot"><span class="foot-chip"><i class="dot"></i> Interactive molecular canvas</span><span>Drag to rotate · Scroll to zoom · Right-drag to pan</span></div>
</div>
<script>
(function(){{
  const root=document.getElementById('viewer');
  if(!window.$3Dmol){{root.innerHTML='<div class="loadmsg">The 3D viewer library could not load. Check your internet connection and refresh the page.</div>';return;}}
  const pdbText=atob({pdb_json});
  const mode={mode_json}; const color={color_json};
  root.innerHTML='';
  const viewer=$3Dmol.createViewer(root,{{backgroundColor:'#0b1428',antialias:true}});
  viewer.addModel(pdbText,'pdb');
  let style={{}};
  if(mode==='cartoon') style.cartoon={{}};
  else if(mode==='stick') style.stick={{radius:0.17}};
  else style.sphere={{scale:0.28}};
  if(color==='rainbow') style[mode].color='spectrum';
  else if(color==='chain') style[mode].colorscheme='chain';
  else style[mode].colorscheme={{prop:'b',gradient:'roygb',min:0,max:100}};
  viewer.setStyle({{}},style);
  viewer.zoomTo(); viewer.render();
  window.addEventListener('resize',()=>viewer.resize());
}})();
</script></body></html>"""


def build_report(sequence: str, metrics: dict[str, float | int] | None,
                 pdb_info: dict[str, object] | None, filename: str | None,
                 source_label: str) -> str:
    lines = [
        f"{APP_NAME} — Analysis Summary",
        f"Generated (UTC): {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}",
        f"Displayed structure source: {source_label}",
        "",
        "IMPORTANT: These outputs are computational summaries, not experimental evidence or a full structural validation.",
        "",
    ]
    if sequence:
        lines += ["SEQUENCE", f"Length: {len(sequence)} amino-acid residues", "", "Sequence:", sequence, ""]
    if metrics:
        lines += ["SEQUENCE METRICS"]
        for key, value in metrics.items():
            lines.append(f"{key}: {value:.4f}" if isinstance(value, float) else f"{key}: {value}")
        lines += ["Interpretation note: The instability index is a sequence-derived heuristic; it does not establish real-world protein stability.", ""]
    if pdb_info:
        lines += ["PDB FILE SUMMARY", f"File: {filename or 'displayed structure'}"]
        for key, value in pdb_info.items():
            if key == "Mean B-factor field" and isinstance(value, float):
                lines.append(f"{key}: {value:.3f}")
            else:
                lines.append(f"{key}: {value}")
        lines += ["Interpretation note: PDB B-factor fields have structure-file-specific meanings; they are not universally equivalent to confidence scores.", ""]
    lines += [
        "Suggested next steps:",
        "1. Check sequence similarity with NCBI BLASTp.",
        "2. Compare available experimental structures in RCSB PDB.",
        "3. Run an appropriate structural validation workflow (e.g., Ramachandran/PROCHECK, Verify3D, ERRAT) before drawing conclusions.",
    ]
    return "\n".join(lines)


# ── Page header ────────────────────────────────────────────────────────────────
st.markdown(
    """
    <section class="hero-glass">
      <div class="eyebrow">Bioinformatics · Protein intelligence workspace</div>
      <div class="hero-title">See the structure.<br><span class="gradient-text">Understand the protein.</span></div>
      <p class="hero-subtitle">Explore molecular structures in an interactive 3D canvas, inspect protein sequence properties, and export clear analysis summaries — all in one research workspace.</p>
      <div class="status-row">
        <span class="status-pill"><i class="status-dot"></i> 3D structure explorer</span>
        <span class="status-pill">FASTA sequence analytics</span>
        <span class="status-pill">PDB diagnostics & exports</span>
      </div>
    </section>
    """,
    unsafe_allow_html=True,
)

# ── User inputs ────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown('<div class="eyebrow">WORKSPACE INPUTS</div>', unsafe_allow_html=True)
    st.markdown("### Load your protein")
    st.caption("Upload a PDB structure to replace the demo model. You can also paste a raw or FASTA-formatted amino-acid sequence for analysis.")
    uploaded_pdb = st.file_uploader("Upload structure (.pdb / .ent)", type=["pdb", "ent"], accept_multiple_files=False)
    st.markdown('<div class="soft-rule"></div>', unsafe_allow_html=True)
    sequence_input = st.text_area(
        "Protein sequence (FASTA or raw)", value=SAMPLE_SEQUENCE, height=175,
        help="FASTA headers beginning with > and whitespace are removed automatically.",
    )
    st.markdown('<div class="tiny-muted">Prediction is performed by a separate Colab or external service. This website focuses on visualization and analysis.</div>', unsafe_allow_html=True)
    st.markdown('<div class="soft-rule"></div>', unsafe_allow_html=True)
    st.markdown("**External research tools**")
    st.markdown("[RCSB PDB](https://www.rcsb.org/) · [UniProt](https://www.uniprot.org/) · [NCBI BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi)")

sequence = clean_sequence(sequence_input)
invalid = sorted(set(sequence) - set(STANDARD_AA))
seq_valid = bool(sequence) and not invalid
metrics: dict[str, float | int] | None = None
composition: Counter[str] = Counter(sequence) if seq_valid else Counter()
if seq_valid:
    try:
        metrics = sequence_metrics(sequence)
    except Exception as exc:  # shown as a user-facing error
        st.error(f"Could not calculate sequence properties: {exc}")

# Read a user upload when present; otherwise keep the homepage useful with a bundled demo fragment.
pdb_text: str | None = None
pdb_info: dict[str, object] | None = None
pdb_filename: str | None = None
pdb_error: str | None = None
source_label = "Demo · ubiquitin fragment (17–33)"
using_uploaded = False
if uploaded_pdb is not None:
    pdb_filename = uploaded_pdb.name
    source_label = uploaded_pdb.name
    if uploaded_pdb.size > MAX_PDB_BYTES:
        pdb_error = "The uploaded file is larger than the 20 MiB safety limit. The demo structure remains available."
    else:
        try:
            candidate = uploaded_pdb.getvalue().decode("utf-8", errors="replace")
            if not any(line.startswith(("ATOM  ", "HETATM")) for line in candidate.splitlines()):
                pdb_error = "This file does not contain PDB ATOM/HETATM records. The demo structure remains available."
            else:
                pdb_text = candidate
                pdb_info = pdb_summary(pdb_text)
                using_uploaded = True
        except Exception as exc:
            pdb_error = f"Could not read the uploaded file: {exc}. The demo structure remains available."

if not using_uploaded:
    source_label = "Demo · ubiquitin fragment (17–33)"
    try:
        pdb_text = DEMO_PDB_PATH.read_text(encoding="utf-8", errors="replace")
        pdb_info = pdb_summary(pdb_text)
        pdb_filename = "example_ubiquitin_fragment.pdb"
    except OSError:
        pdb_error = (pdb_error + " " if pdb_error else "") + "The bundled demo structure file is missing. Add data/example_ubiquitin_fragment.pdb to the repository."

if invalid:
    st.warning("Sequence contains unsupported symbols: " + ", ".join(invalid) + ". Use the 20 standard one-letter amino-acid codes (ACDEFGHIKLMNPQRSTVWY).")
elif not sequence:
    st.info("Enter an amino-acid sequence in the sidebar to calculate sequence properties.")
if pdb_error:
    st.warning(pdb_error)

source_is_demo = not using_uploaded

# ── Main molecular workspace ───────────────────────────────────────────────────
st.markdown('<div class="section-kicker">01 / Molecular workspace</div>', unsafe_allow_html=True)
st.markdown('<div class="section-note">Your structure is the centrepiece. Use the controls beside the canvas to explore different molecular representations.</div>', unsafe_allow_html=True)
viewer_col, control_col = st.columns([2.35, 0.95], gap="large", vertical_alignment="top")
with control_col:
    st.markdown('<div class="section-kicker">VIEWER CONTROLS</div>', unsafe_allow_html=True)
    render_style = st.selectbox("Representation", ["Cartoon", "Stick", "Sphere"], index=0)
    color_mode = st.selectbox("Color mapping", ["Rainbow", "By chain", "B-factor field"], index=0)
    st.markdown('<div class="soft-rule"></div>', unsafe_allow_html=True)
    st.markdown('<div class="section-kicker">STRUCTURE SNAPSHOT</div>', unsafe_allow_html=True)
    if pdb_info:
        st.metric("Atom records", f"{pdb_info['Atom records']:,}")
        st.metric("Observed residues", f"{pdb_info['Observed residues']:,}")
        st.metric("Protein chains", f"{len(pdb_info['Chains']):,}")
        st.caption("Counts are parsed from the uploaded or demo PDB file.")
    if source_is_demo:
        st.info("Demo mode: this bundled file is an experimental ubiquitin fragment (residues 17–33), not the full protein. Upload your own PDB to explore another structure.")
    else:
        st.success("Custom PDB loaded")

with viewer_col:
    if pdb_text:
        components.html(make_viewer_html(pdb_text, render_style, color_mode, source_label), height=690, scrolling=False)
    else:
        st.error("No valid PDB file is available. Add the demo PDB file to the repository or upload a valid structure.")

# Small quick-read metrics under the 3D viewer
st.markdown('<div class="soft-rule"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-kicker">02 / Quick readout</div>', unsafe_allow_html=True)
m1, m2, m3, m4 = st.columns(4)
m1.metric("Sequence length", f"{len(sequence):,} aa" if seq_valid else "—")
m2.metric("Molecular weight", f"{metrics['Molecular weight (Da)']:,.1f} Da" if metrics else "—")
m3.metric("Theoretical pI", f"{metrics['Theoretical pI']:.2f}" if metrics else "—")
m4.metric("Displayed model", "Demo fragment" if source_is_demo else "Uploaded PDB")

# ── Detailed tools remain available below the viewer ────────────────────────────
st.markdown('<div class="soft-rule"></div>', unsafe_allow_html=True)
tab_sequence, tab_quality, tab_workflow, tab_export = st.tabs(
    ["Sequence analytics", "Model diagnostics", "Workflow & resources", "Downloads & report"]
)

with tab_sequence:
    st.markdown('<div class="section-kicker">SEQUENCE LAB</div>', unsafe_allow_html=True)
    st.subheader("Primary sequence profile")
    if not sequence:
        st.info("Paste a sequence in the sidebar to begin.")
    elif invalid:
        st.warning("Fix the unsupported sequence symbols before running the calculations.")
    elif metrics:
        a, b, c = st.columns(3)
        a.metric("Length", f"{metrics['Length (residues)']:,} residues")
        b.metric("Molecular weight", f"{metrics['Molecular weight (Da)']:,.2f} Da")
        c.metric("Theoretical pI", f"{metrics['Theoretical pI']:.2f}")
        d, e, f = st.columns(3)
        d.metric("Aromaticity", f"{metrics['Aromaticity']:.3f}")
        e.metric("Instability index", f"{metrics['Instability index']:.2f}")
        f.metric("GRAVY hydropathy", f"{metrics['GRAVY hydropathy']:.3f}")
        st.caption("These are sequence-derived descriptors. The instability index is a heuristic; none of these values proves that a protein is stable or biologically functional.")
        st.markdown("#### Amino-acid composition")
        comp_df = pd.DataFrame({"Residues": [composition.get(aa, 0) for aa in STANDARD_AA]}, index=list(STANDARD_AA))
        st.bar_chart(comp_df, height=300)

with tab_quality:
    st.markdown('<div class="section-kicker">STRUCTURE QUALITY / BASIC CHECKS</div>', unsafe_allow_html=True)
    st.subheader("PDB file diagnostics")
    if pdb_text and pdb_info:
        q1, q2, q3 = st.columns(3)
        q1.metric("Atom records present", "Yes")
        q2.metric("Observed residues", f"{pdb_info['Observed residues']:,}")
        q3.metric("Chains found", f"{len(pdb_info['Chains']):,}")
        st.write("**Chain labels:**", ", ".join(pdb_info["Chains"]) if pdb_info["Chains"] else "No chain labels found")
        st.write("**PDB file size:**", f"{pdb_info['File size (KiB)']:.1f} KiB")
        if pdb_info["Mean B-factor field"] is not None:
            st.metric("Mean B-factor field", f"{pdb_info['Mean B-factor field']:.2f}")
            st.caption("B-factor values are reported as present in the PDB file. They are not universally equivalent to prediction-confidence scores.")
        st.warning("These are basic file diagnostics, not full structural validation. This app does not calculate Ramachandran favored-region percentages, ERRAT, Verify3D, clash scores, or experimental accuracy.")
    else:
        st.info("Upload a PDB file to run basic file diagnostics.")
    st.markdown("#### Suggested external validation workflow")
    st.markdown(
        "1. **SAVES / PROCHECK** — review stereochemistry and Ramachandran regions.\n"
        "2. **ERRAT and Verify3D** — examine non-bonded interactions and sequence–structure compatibility.\n"
        "3. **Model-specific confidence** — use the metric supplied by the prediction method; do not infer it from viewer color alone.\n"
        "4. **Biological context** — compare templates, sequence identity, missing residues, oligomeric state, and relevant literature."
    )
    st.markdown("[UCLA SAVES](https://saves.mbi.ucla.edu/) · [RCSB PDB](https://www.rcsb.org/) · [AlphaFold DB](https://alphafold.ebi.ac.uk/)")

with tab_workflow:
    st.markdown('<div class="section-kicker">RESEARCH PIPELINE</div>', unsafe_allow_html=True)
    st.subheader("From sequence to structural insight")
    step1, step2, step3 = st.columns(3)
    with step1:
        st.markdown("### 01 · Sequence")
        st.write("Retrieve a protein sequence from UniProt or NCBI Protein. Paste it into the sidebar to calculate sequence-derived properties.")
    with step2:
        st.markdown("### 02 · Structure")
        st.write("Obtain an experimental structure from RCSB PDB or run a compatible prediction workflow in Google Colab. Upload the resulting PDB file.")
    with step3:
        st.markdown("### 03 · Review")
        st.write("Inspect the 3D model, review basic PDB diagnostics, export results, and use dedicated validation tools before drawing conclusions.")
    st.info("The public website does not run a large folding model on its hosting server. Structure prediction depends on separate external services and may be rate-limited or unavailable.")
    st.markdown("#### Quick links")
    st.markdown(
        "- [RCSB Protein Data Bank](https://www.rcsb.org/) — browse experimental structures.\n"
        "- [UniProt](https://www.uniprot.org/) — obtain curated protein sequence information.\n"
        "- [NCBI BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi) — look for sequence similarities.\n"
        "- [ESM Atlas](https://esmatlas.com/resources?action=fold) — external sequence-to-structure prediction interface.\n"
        "- [ColabFold notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb) — alternative notebook workflow.\n"
        "- [UCLA SAVES](https://saves.mbi.ucla.edu/) — structural validation resources."
    )

with tab_export:
    st.markdown('<div class="section-kicker">EXPORT CENTRE</div>', unsafe_allow_html=True)
    st.subheader("Download your analysis outputs")
    if seq_valid:
        fasta_text = f">BioFold_Workbench_sequence length={len(sequence)}\n" + "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70)) + "\n"
        st.download_button("Download cleaned FASTA", data=fasta_text, file_name="cleaned_sequence.fasta", mime="text/plain")
    if pdb_text:
        st.download_button("Download displayed PDB", data=pdb_text, file_name=pdb_filename or "displayed_structure.pdb", mime="chemical/x-pdb")
    report = build_report(sequence if seq_valid else "", metrics, pdb_info, pdb_filename, source_label)
    st.download_button("Download analysis report (.txt)", data=report, file_name="biofold_analysis_report.txt", mime="text/plain")
    st.caption("Use external structures and prediction outputs according to their terms and attribution requirements. Predicted structures are hypotheses and should not be presented as experimentally verified models.")

st.markdown('<div class="soft-rule"></div>', unsafe_allow_html=True)
st.markdown('<div class="tiny-muted">BioFold Workbench · Educational bioinformatics project · Glassmorphism UI · Structure prediction and external validation services are not guaranteed to be available at all times.</div>', unsafe_allow_html=True)
