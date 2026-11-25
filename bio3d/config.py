"""Shared constants, reference tables, and sample presets for BioFold Workbench."""

from __future__ import annotations

from pathlib import Path

APP_NAME = "BioFold Workbench"
APP_VERSION = "2.0.0"

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "data"

MAX_STRUCTURE_BYTES = 20 * 1024 * 1024  # 20 MiB
MAX_FASTA_BYTES = 2 * 1024 * 1024  # 2 MiB
MAX_ESMFOLD_RESIDUES = 400
MIN_ESMFOLD_RESIDUES = 10

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AA_SET = frozenset(STANDARD_AA)

AA3_TO_AA1: dict[str, str] = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    # Common protonation / selenomethionine mappings for polymer continuity
    "MSE": "M",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "CYX": "C",
}

STANDARD_RESIDUE_NAMES = frozenset(AA3_TO_AA1.keys())
WATER_RESIDUE_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})

AA_FULL_NAMES: dict[str, str] = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartate",
    "E": "Glutamate",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
}

AA_GROUPS: dict[str, str] = {
    "A": "Hydrophobic",
    "V": "Hydrophobic",
    "I": "Hydrophobic",
    "L": "Hydrophobic",
    "M": "Hydrophobic",
    "F": "Aromatic",
    "W": "Aromatic",
    "Y": "Aromatic",
    "S": "Polar",
    "T": "Polar",
    "N": "Polar",
    "Q": "Polar",
    "C": "Special",
    "G": "Special",
    "P": "Special",
    "K": "Basic (+)",
    "R": "Basic (+)",
    "H": "Basic (+)",
    "D": "Acidic (-)",
    "E": "Acidic (-)",
}

GROUP_COLORS: dict[str, str] = {
    "Hydrophobic": "#38bdf8",
    "Aromatic": "#a78bfa",
    "Polar": "#2dd4bf",
    "Basic (+)": "#60a5fa",
    "Acidic (-)": "#f87171",
    "Special": "#fbbf24",
}

# Kyte & Doolittle (1982), J. Mol. Biol. 157:105-132
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

BUNDLED_STRUCTURES: dict[str, dict[str, str]] = {
    "1UBQ — Human Ubiquitin (76 aa, 1.80 Å)": {
        "pdb_id": "1UBQ",
        "path": str(DATA_DIR / "1ubq_ubiquitin.pdb"),
        "description": "High-resolution X-ray crystal structure of human ubiquitin (PDB: 1UBQ).",
    },
    "1CRN — Crambin (46 aa, 1.50 Å)": {
        "pdb_id": "1CRN",
        "path": str(DATA_DIR / "1crn_crambin.pdb"),
        "description": "Ultra-high-resolution X-ray structure of plant seed protein crambin (PDB: 1CRN).",
    },
}

UBIQUITIN_SEQUENCE = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)

GFP_SEQUENCE = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPD"
    "HMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYI"
    "MADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAG"
    "ITHGMDELYK"
)

HBA_HUMAN_SEQ = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDD"
    "MPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
)

HBB_HUMAN_SEQ = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGL"
    "AHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
)

ORIGIN_LABELS: dict[str, str] = {
    "experimental": "Experimental Structure (PDB / RCSB)",
    "predicted_alphafold": "AI-Predicted Structure (AlphaFold DB)",
    "predicted_esmfold": "AI-Predicted Structure (ESMFold API)",
    "uploaded_experimental": "Uploaded Experimental Structure",
    "uploaded_predicted": "Uploaded AI-Predicted Structure",
    "uploaded_unknown": "Uploaded Structure (Origin Unspecified)",
}
