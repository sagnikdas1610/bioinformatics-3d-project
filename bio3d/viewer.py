"""Interactive 3Dmol.js molecular viewer HTML builder."""

from __future__ import annotations

import base64
from dataclasses import dataclass, field
from html import escape
import json
from typing import Literal

RenderStyle = Literal["Cartoon", "Stick", "Ball & Stick", "Sphere"]
ColorMode = Literal[
    "Secondary Structure",
    "By Chain",
    "N->C Spectrum",
    "Experimental B-factor",
    "Predicted pLDDT",
]


@dataclass(frozen=True)
class ViewerConfig:
    """Configuration options for the interactive 3Dmol.js viewport."""

    render_style: RenderStyle = "Cartoon"
    color_mode: ColorMode = "Secondary Structure"
    file_format: Literal["pdb", "cif"] = "pdb"
    source_label: str = "Structure"
    is_predicted: bool = False
    plddt_scale_0_1: bool = False
    bfactor_min: float = 5.0
    bfactor_max: float = 80.0
    selected_chain: str = "All"
    highlight_residues: list[int] = field(default_factory=list)
    highlight_chain: str | None = None
    show_surface: bool = False
    surface_opacity: float = 0.68
    show_ligands: bool = True
    show_labels: bool = True
    zoom_to_selection: bool = False
    atom_count: int = 0


def _legend_html(config: ViewerConfig) -> str:
    if config.color_mode == "Secondary Structure":
        return (
            '<div class="legend">'
            '<span class="leg"><i style="background:#ec4899"></i>α-Helix</span>'
            '<span class="leg"><i style="background:#f59e0b"></i>β-Sheet</span>'
            '<span class="leg"><i style="background:#38bdf8"></i>Coil / Turn</span>'
            "</div>"
        )
    if config.color_mode == "Predicted pLDDT":
        return (
            '<div class="legend">'
            '<span class="leg"><i style="background:#0053d6"></i>≥90 Very high</span>'
            '<span class="leg"><i style="background:#65cbf3"></i>70–90 Confident</span>'
            '<span class="leg"><i style="background:#ffdb13"></i>50–70 Low</span>'
            '<span class="leg"><i style="background:#ff7d45"></i>&lt;50 Very low</span>'
            "</div>"
        )
    if config.color_mode == "Experimental B-factor":
        return (
            '<div class="legend">'
            f'<span>B-factor ({config.bfactor_min:.0f}–{config.bfactor_max:.0f} Å²):</span>'
            '<span class="grad-b"></span>'
            '<span>Low (ordered) → High (flexible)</span>'
            "</div>"
        )
    if config.color_mode == "N->C Spectrum":
        return (
            '<div class="legend">'
            '<span>N-terminus</span><span class="grad-spec"></span><span>C-terminus</span>'
            "</div>"
        )
    return '<div class="legend"><span>Colored by polymer chain ID</span></div>'


def build_viewer_html(raw_text: str, config: ViewerConfig) -> str:
    """Build a self-contained UTF-8 safe 3Dmol.js HTML canvas with interactive controls."""
    encoded_b64 = base64.b64encode(raw_text.encode("utf-8")).decode("ascii")
    payload = {
        "b64": encoded_b64,
        "fmt": config.file_format,
        "style": config.render_style,
        "color": config.color_mode,
        "chain": config.selected_chain,
        "highlightResidues": [int(r) for r in config.highlight_residues],
        "highlightChain": config.highlight_chain,
        "showSurface": bool(config.show_surface and config.atom_count <= 20000),
        "surfaceOpacity": float(config.surface_opacity),
        "showLigands": bool(config.show_ligands),
        "showLabels": bool(config.show_labels),
        "zoomToSelection": bool(config.zoom_to_selection),
        "bMin": float(config.bfactor_min),
        "bMax": float(max(config.bfactor_max, config.bfactor_min + 1.0)),
        "plddtScale100": bool(config.plddt_scale_0_1),
        "filename": config.source_label,
    }
    cfg_json = json.dumps(payload)
    source_html = escape(config.source_label)
    legend_html = _legend_html(config)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"
        onerror="var s=document.createElement('script');s.src='https://3dmol.org/build/3Dmol-min.js';document.head.appendChild(s);"></script>
<style>
  * {{ box-sizing: border-box; }}
  html, body {{
    width: 100%; height: 100%; margin: 0; padding: 0;
    background: #0b1120; color: #e2e8f0;
    font-family: Inter, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif;
    overflow: hidden;
  }}
  .shell {{
    height: 100%; display: flex; flex-direction: column;
    border: 1px solid #1e293b; border-radius: 14px; overflow: hidden;
    background: #0b1120;
  }}
  .topbar {{
    display: flex; align-items: center; justify-content: space-between; gap: 10px;
    padding: 10px 14px; background: #111827; border-bottom: 1px solid #1e293b;
    flex-wrap: wrap;
  }}
  .meta-left {{ display: flex; align-items: center; gap: 10px; min-width: 0; }}
  .badge {{
    font-size: 11px; font-weight: 600; padding: 3px 9px; border-radius: 999px;
    background: rgba(56, 189, 248, 0.14); color: #7dd3fc; border: 1px solid rgba(56, 189, 248, 0.35);
    white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 280px;
  }}
  .controls {{ display: flex; align-items: center; gap: 6px; }}
  .btn {{
    background: #1e293b; color: #e2e8f0; border: 1px solid #334155;
    border-radius: 7px; padding: 5px 10px; font-size: 11.5px; font-weight: 600;
    cursor: pointer; transition: border-color 0.15s, background 0.15s;
  }}
  .btn:hover {{ border-color: #38bdf8; background: #263349; }}
  .btn.active {{ border-color: #38bdf8; color: #7dd3fc; }}
  #viewport {{
    position: relative; flex: 1; width: 100%; min-height: 360px; background: #0b1120;
  }}
  .foot {{
    display: flex; align-items: center; justify-content: space-between; gap: 10px; flex-wrap: wrap;
    padding: 8px 14px; background: #111827; border-top: 1px solid #1e293b;
    font-size: 11px; color: #94a3b8;
  }}
  .legend {{ display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }}
  .leg {{ display: inline-flex; align-items: center; gap: 5px; color: #cbd5e1; }}
  .leg i {{ width: 10px; height: 10px; border-radius: 3px; display: inline-block; }}
  .grad-b {{
    width: 72px; height: 8px; border-radius: 4px; display: inline-block;
    background: linear-gradient(90deg, #3b82f6, #f8fafc, #ef4444);
  }}
  .grad-spec {{
    width: 72px; height: 8px; border-radius: 4px; display: inline-block;
    background: linear-gradient(90deg, #3b82f6, #22c55e, #eab308, #ef4444);
  }}
  .status-msg {{ padding: 24px; font-size: 13px; color: #94a3b8; }}
</style>
</head>
<body>
<div class="shell">
  <div class="topbar">
    <div class="meta-left">
      <span style="font-size:13px;font-weight:700;color:#f1f5f9;">3D Molecular Viewport</span>
      <span class="badge" title="{source_html}">{source_html}</span>
    </div>
    <div class="controls">
      <button class="btn" id="btn-reset" type="button">Reset View</button>
      <button class="btn" id="btn-spin" type="button">Spin</button>
      <button class="btn" id="btn-png" type="button">Save PNG</button>
    </div>
  </div>
  <div id="viewport"><div class="status-msg">Loading 3Dmol.js molecular canvas…</div></div>
  <div class="foot">
    {legend_html}
    <span>Drag: rotate · Scroll: zoom · Right-drag: pan</span>
  </div>
</div>

<script>
(function() {{
  const cfg = {cfg_json};
  const root = document.getElementById("viewport");

  function decodeBase64Utf8(b64) {{
    const binary = atob(b64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
    return new TextDecoder("utf-8").decode(bytes);
  }}

  function plddtColor(atom) {{
    const raw = Number(atom.b || 0);
    const val = cfg.plddtScale100 ? raw * 100.0 : raw;
    if (val >= 90) return 0x0053d6;
    if (val >= 70) return 0x65cbf3;
    if (val >= 50) return 0xffdb13;
    return 0xff7d45;
  }}

  function ssColor(atom) {{
    if (atom.ss === "h") return 0xec4899;
    if (atom.ss === "s") return 0xf59e0b;
    return 0x38bdf8;
  }}

  function makeStyleObj(styleName, colorMode) {{
    const s = {{}};
    const colorSpec = {{}};
    if (colorMode === "Secondary Structure") colorSpec.colorfunc = ssColor;
    else if (colorMode === "By Chain") colorSpec.colorscheme = "chain";
    else if (colorMode === "N->C Spectrum") colorSpec.color = "spectrum";
    else if (colorMode === "Predicted pLDDT") colorSpec.colorfunc = plddtColor;
    else colorSpec.colorscheme = {{ prop: "b", gradient: "bwr", min: cfg.bMin, max: cfg.bMax }};

    if (styleName === "Cartoon") {{
      s.cartoon = Object.assign({{}}, colorSpec);
    }} else if (styleName === "Stick") {{
      s.stick = Object.assign({{ radius: 0.16 }}, colorSpec);
    }} else if (styleName === "Ball & Stick") {{
      s.stick = Object.assign({{ radius: 0.13 }}, colorSpec);
      s.sphere = Object.assign({{ scale: 0.23 }}, colorSpec);
    }} else {{
      s.sphere = Object.assign({{ scale: 0.34 }}, colorSpec);
    }}
    return s;
  }}

  function initViewer() {{
    if (!window.$3Dmol) {{
      root.innerHTML = '<div class="status-msg">Could not load 3Dmol.js from CDN. Check internet connectivity and refresh.</div>';
      return;
    }}
    root.innerHTML = "";
    const viewer = $3Dmol.createViewer(root, {{ backgroundColor: "#0b1120", antialias: true }});
    const text = decodeBase64Utf8(cfg.b64);
    viewer.addModel(text, cfg.fmt === "cif" ? "cif" : "pdb");

    // Base polymer selector
    const polySel = cfg.chain && cfg.chain !== "All"
      ? {{ chain: cfg.chain, hetflag: false }}
      : {{ hetflag: false }};

    // Dim non-selected chains gently if a specific chain is focused
    if (cfg.chain && cfg.chain !== "All") {{
      viewer.setStyle({{ hetflag: false }}, {{ cartoon: {{ color: "#334155", opacity: 0.25 }} }});
    }}

    viewer.setStyle(polySel, makeStyleObj(cfg.style, cfg.color));

    // Hetero non-water ligands
    if (cfg.showLigands) {{
      const ligSel = cfg.chain && cfg.chain !== "All"
        ? {{ hetflag: true, chain: cfg.chain, not: {{ resn: ["HOH", "WAT", "DOD"] }} }}
        : {{ hetflag: true, not: {{ resn: ["HOH", "WAT", "DOD"] }} }};
      viewer.setStyle(ligSel, {{ stick: {{ radius: 0.22, colorscheme: "Jmol" }}, sphere: {{ scale: 0.28, colorscheme: "Jmol" }} }});
    }}

    // Translucent molecular surface
    if (cfg.showSurface) {{
      viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: cfg.surfaceOpacity, color: "#1e293b" }}, polySel);
    }}

    // Highlight selected residues
    let zoomTarget = polySel;
    if (cfg.highlightResidues && cfg.highlightResidues.length > 0) {{
      const hiSel = {{ hetflag: false, resi: cfg.highlightResidues }};
      if (cfg.highlightChain && cfg.highlightChain !== "All") hiSel.chain = cfg.highlightChain;
      else if (cfg.chain && cfg.chain !== "All") hiSel.chain = cfg.chain;

      viewer.addStyle(hiSel, {{
        stick: {{ radius: 0.24, color: "#facc15" }},
        sphere: {{ scale: 0.30, color: "#facc15" }}
      }});

      if (cfg.showLabels) {{
        const atoms = viewer.selectedAtoms(Object.assign({{ atom: "CA" }}, hiSel));
        atoms.forEach(a => {{
          viewer.addLabel(a.resn + a.resi, {{
            position: {{ x: a.x, y: a.y, z: a.z }},
            backgroundColor: "rgba(15,23,42,0.88)",
            borderColor: "#facc15",
            borderThickness: 1,
            fontColor: "#fef08a",
            fontSize: 11,
            inFront: true
          }});
        }});
      }}
      if (cfg.zoomToSelection) zoomTarget = hiSel;
    }}

    viewer.zoomTo(zoomTarget);
    viewer.render();

    let spinning = false;
    document.getElementById("btn-reset").addEventListener("click", () => {{
      viewer.zoomTo(zoomTarget);
      viewer.render();
    }});
    const spinBtn = document.getElementById("btn-spin");
    spinBtn.addEventListener("click", () => {{
      spinning = !spinning;
      viewer.spin(spinning ? "y" : false, 0.7);
      spinBtn.classList.toggle("active", spinning);
    }});
    document.getElementById("btn-png").addEventListener("click", () => {{
      const uri = viewer.pngURI();
      const a = document.createElement("a");
      a.href = uri;
      a.download = (cfg.filename || "structure").replace(/[^a-zA-Z0-9._-]+/g, "_") + ".png";
      a.click();
    }});
    window.addEventListener("resize", () => viewer.resize());
  }}

  if (window.$3Dmol) initViewer();
  else window.addEventListener("load", initViewer);
}})();
</script>
</body>
</html>"""
