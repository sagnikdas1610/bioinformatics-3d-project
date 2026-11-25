# BioFold Workbench — Bioinformatics 3D Protein Structure & Analytics Platform

[![Open Live App](https://img.shields.io/badge/Streamlit-Live_App-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white)](https://bioinformatics-3d-project-ai.streamlit.app/)
[![Open in Colab](https://img.shields.io/badge/Google_Colab-Prediction_Notebook-F9AB00?style=for-the-badge&logo=googlecolab&logoColor=white)](https://colab.research.google.com/github/sagnikdas1610/bioinformatics-3d-project/blob/main/AI_Structure_Generator.ipynb)
![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat&logo=python&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

> **Live Application:** [https://bioinformatics-3d-project-ai.streamlit.app/](https://bioinformatics-3d-project-ai.streamlit.app/)
> **Companion Prediction Notebook:** [AI_Structure_Generator.ipynb](https://colab.research.google.com/github/sagnikdas1610/bioinformatics-3d-project/blob/main/AI_Structure_Generator.ipynb)

---

## Overview

**BioFold Workbench** is an open-source computational structural biology and protein analytics platform built with Python and Streamlit. It brings together interactive 3D molecular visualization, primary sequence biophysics, pairwise sequence alignment, coordinate-based structural validation (Ramachandran plots, Kabsch–Sander secondary structure, residue contact maps), transparent AI structure prediction integration (AlphaFold DB & ESMFold API), and downloadable scientific PDF reports in a single research-oriented workspace.

---

## Key Features

- **Interactive 3D Molecular Explorer (`3Dmol.js`)**
  - Supports **Cartoon**, **Stick**, **Ball & Stick**, and **Sphere** representations plus translucent **molecular surface** rendering.
  - Color modes: **Secondary Structure**, **By Chain**, **N→C Spectrum**, **Experimental B-factor (Å²)**, and **Predicted pLDDT (0–100)** — strictly separated by structure provenance.
  - Chain focus selector, residue number/range highlighting with labels, hetero-ligand display, camera reset, spin toggle, and PNG viewport export.
  - Supports **PDB** (`.pdb`, `.ent`), **mmCIF** (`.cif`, `.mmcif`), bundled benchmark samples (`1UBQ`, `1CRN`), and live **RCSB PDB** fetching by 4-character ID.

- **Primary Sequence Analytics**
  - FASTA file upload, raw sequence input, and multi-FASTA record selection (never concatenating distinct entries).
  - One-click sequence synchronization directly from any loaded 3D structure chain.
  - Biopython `ProtParam` calculation of sequence length, molecular weight, theoretical pI, net charge at pH 7.0, GRAVY hydropathy, instability index, aromaticity, and molar extinction coefficients ($\varepsilon_{280}$ reduced/oxidized).
  - Interactive amino-acid composition chart colored by physicochemical class, sliding-window **Kyte–Doolittle hydropathy profile**, and **net charge vs. pH curve**.

- **Pairwise Sequence Alignment Studio**
  - **Global (Needleman–Wunsch)** and **Local (Smith–Waterman)** alignments with **BLOSUM62**, **BLOSUM80**, **BLOSUM45**, **PAM250**, and **IDENTITY** scoring matrices and configurable gap open/extend penalties.
  - Explicitly distinguishes **% Identity** (exact matches) from **% Similarity** (positive substitution-matrix pairs) and **% Gaps**, rendered in a color-coded monospace alignment viewer.

- **Advanced Structural Analysis**
  - **Model-1 & Conformer Discipline:** Separates standard polymer residues from hetero ligands and waters, deduplicates alternate conformations (`altLoc`) by highest occupancy, and evaluates Model 1 on multi-model NMR ensembles.
  - **Missing-Residue Detection:** Parses `REMARK 465` (PDB) and `_pdbx_unobs_or_zero_occ_residues` (mmCIF), plus backbone numbering/peptide-bond continuity checks.
  - **Ramachandran Plot:** Computes backbone $\phi/\psi/\omega$ dihedral angles with peptide-bond continuity checks and evaluates residues against genuine **Richardson Lab MolProbity Top8000** 6-class 2° probability distributions (*General*, *Ile/Val*, *Glycine*, *trans-Proline*, *cis-Proline*, *Pre-Proline*; Williams et al. 2018; Lovell et al. 2003; CC BY 4.0).
  - **Secondary Structure Assignment:** Pure-Python **Kabsch–Sander (1983)** electrostatic backbone H-bond DSSP assignment (with PDB/mmCIF `HELIX`/`SHEET` header comparison when available).
  - **Interactive Contact Map:** Vectorized Cα and Cβ distance/contact matrices with configurable Å cutoff and sequence separation (`|i - j|`), automatically using **Cα fallback for Glycine** in Cβ mode.

- **AI-Assisted Structure Prediction (Honest Provenance)**
  - **AlphaFold DB Lookup:** Fetches precomputed AlphaFold structures by UniProt accession and plots per-residue pLDDT confidence.
  - **Live ESMFold API Folding:** Submits validated single-chain sequences (10–400 standard amino acids) to the public ESM Atlas API (`https://api.esmatlas.com/foldSequence/v1/pdb/`), normalizing pLDDT confidence to 0–100 with zero fake fallbacks on API errors.
  - **ColabFold Links:** Direct workflow links for proteins >400 residues or multi-chain complexes.

- **Publication-Style Scientific PDF & Data Exports**
  - Generates downloadable **PDF reports** (`ReportLab` + `Matplotlib`) containing structure metadata, 3D Cα backbone projection, sequence/hydropathy charts, Ramachandran plot, secondary structure distribution, contact map, methodology, and limitations — explicitly distinguishing `[MEASURED]`, `[CALCULATED]`, `[AI-PREDICTED]`, and `[UNAVAILABLE]` data.
  - Also exports plain-text `.txt` summaries, cleaned `.fasta`, `.pdb`/`.cif`, and alignment `.txt` files.

---

## Architecture

```text
bioinformatics-3d-project/
├── app.py                        # Preserved Streamlit entry point
├── requirements.txt              # Deployment dependencies
├── .streamlit/config.toml        # Unified dark scientific theme
├── data/
│   ├── 1ubq_ubiquitin.pdb        # Experimental Ubiquitin benchmark (76 aa, 1.80 Å)
│   ├── 1crn_crambin.pdb          # Experimental Crambin benchmark (46 aa, 1.50 Å)
│   └── molprobity_top8000_rama.npz # Packed MolProbity Top8000 2-deg reference tables
├── docs/screenshots/             # Verified application screenshots
├── bio3d/                        # Pure-Python scientific core
│   ├── config.py                 # Constants, AA tables, sample metadata
│   ├── sequence.py               # FASTA parser, validation, ProtParam, hydropathy, charge curve
│   ├── alignment.py              # PairwiseAligner / Gotoh affine-gap DP aligner
│   ├── structure.py              # PDB & mmCIF parser (Model 1, altLoc, B-factor vs pLDDT, missing residues)
│   ├── structural_analysis.py    # Ramachandran, Kabsch-Sander SS, Cα/Cβ contact maps
│   ├── viewer.py                 # 3Dmol.js interactive HTML builder
│   ├── prediction.py             # RCSB PDB, AlphaFold DB, ESMFold API clients
│   ├── charts.py                 # Plotly scientific charts
│   └── reporting.py              # ReportLab + Matplotlib PDF & TXT report generator
├── ui/                           # Modular Streamlit UI components & tabs
├── tests/                        # Pytest unit & scientific validation suite
└── AI_Structure_Generator.ipynb  # Companion Google Colab ESMFold notebook
```

---

## Technology Stack

| Component | Technology |
| :--- | :--- |
| Web Interface | Streamlit |
| 3D Molecular Visualization | 3Dmol.js |
| Sequence Biophysics | Biopython (`ProtParam`) |
| Alignment & Numerical Computing | Biopython (`PairwiseAligner`) / Gotoh DP, NumPy, Pandas |
| Interactive Plots | Plotly |
| Scientific PDF Reporting | ReportLab, Matplotlib |
| External Data & Folding APIs | RCSB PDB, EMBL-EBI AlphaFold DB, Meta ESM Atlas ESMFold API |
| Testing | Pytest |

---

## Supported File Formats

- **Structures:** `.pdb`, `.ent`, `.cif`, `.mmcif` (up to 20 MiB)
- **Sequences:** `.fasta`, `.fa`, `.faa`, `.txt`, or pasted raw/FASTA text

---

## Local Installation & Execution

```bash
git clone https://github.com/sagnikdas1610/bioinformatics-3d-project.git
cd bioinformatics-3d-project
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Run the automated test suite:

```bash
python -m pytest -v
```

---

## Scientific Methods, Scope & Limitations

- **No Local Model Inference:** The Streamlit server does not host large deep-learning weights locally; AI predictions rely on external public APIs (AlphaFold DB / ESMFold) or Google Colab, which may be rate-limited or unavailable.
- **Predicted vs. Experimental Provenance:** Computational predictions (AlphaFold DB, ESMFold) are structural hypotheses rather than experimentally determined models. Predicted pLDDT confidence scores (0–100) and crystallographic B-factors ($\text{\AA}^2$) are kept strictly separate.
- **Sequence Heuristics:** Descriptors such as the Guruprasad instability index and Kyte–Doolittle GRAVY are statistical sequence-derived heuristics and do not establish experimental stability or cellular localization.
- **Stereochemical Validation:** Built-in Ramachandran (Lovell et al., 2003) and Kabsch–Sander secondary structure calculations provide rapid geometric diagnostics, not a full replacement for complete crystallographic validation suites such as MolProbity or RCSB validation reports.

---

## Screenshots

Verified captures from the running platform are stored in [`docs/screenshots/`](docs/screenshots/):

- [3D Molecular Explorer (`4HHB` Chain A)](docs/screenshots/01_3d_explorer.png)
- [3D Explorer — AlphaFold `AF-P69905-F1` pLDDT coloring](docs/screenshots/08_3d_plddt_alphafold.png)
- [Sequence Analytics (`1UBQ`)](docs/screenshots/02_sequence_analytics.png)
- [Pairwise Sequence Alignment (`HBA` vs. `HBB`)](docs/screenshots/03_pairwise_alignment.png)
- [Structural Analysis — Ramachandran Plot (MolProbity Top8000)](docs/screenshots/04_structural_ramachandran.png)
- [Structural Analysis — Cβ Contact Map](docs/screenshots/05_structural_contact_map.png)
- [AI Structure Prediction (AlphaFold DB & ESMFold API)](docs/screenshots/06_ai_prediction.png)
- [Scientific PDF Reports & Exports](docs/screenshots/07_reports_exports.png)

---

## Deployment

The live application is deployed on **Streamlit Community Cloud** at [https://bioinformatics-3d-project-ai.streamlit.app/](https://bioinformatics-3d-project-ai.streamlit.app/) using root `app.py` and `requirements.txt`.

---

## Data Attribution & Licensing

- **Ramachandran Reference Tables:** `data/molprobity_top8000_rama.npz` packs the six Top8000 2°-binned probability tables (`General`, `Ile/Val`, `Glycine`, `trans-Proline`, `cis-Proline`, `Pre-Proline`) from the **Richardson Laboratory, Duke University** ([`rlabduke/reference_data`](https://github.com/rlabduke/reference_data)), used under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license.
  - Williams CJ et al. (2018), *MolProbity: More and better reference data for improved all-atom structure validation*, *Protein Science* 27:293–315.
  - Lovell SC et al. (2003), *Structure validation by Cα geometry: φ,ψ and Cβ deviation*, *Proteins* 50:437–450.
- **Benchmark Structures:** Experimental coordinates for `1UBQ` and `1CRN` are sourced from the [RCSB Protein Data Bank](https://www.rcsb.org/).

---

## Author & License

**Sagnik Das** · GitHub: [@sagnikdas1610](https://github.com/sagnikdas1610)
Released under the [MIT License](LICENSE).
