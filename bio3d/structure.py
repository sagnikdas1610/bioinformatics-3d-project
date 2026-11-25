"""PDB and mmCIF structure parsing with provenance, model-0, altLoc, and missing-residue support."""

from __future__ import annotations

from dataclasses import dataclass, field
import re
import shlex
from typing import Literal

import numpy as np

from bio3d.config import AA3_TO_AA1, MODIFIED_AA3_TO_AA1, STANDARD_RESIDUE_NAMES, WATER_RESIDUE_NAMES


@dataclass
class AtomRecord:
    """Single atomic coordinate record."""

    name: str
    coord: np.ndarray
    bfactor: float
    occupancy: float
    altloc: str
    element: str


@dataclass
class ResidueRecord:
    """Single residue from Model 0 with primary-conformer atoms."""

    chain: str
    resseq: int
    icode: str
    resname: str
    aa1: str
    is_hetero: bool
    is_water: bool
    atoms: dict[str, AtomRecord] = field(default_factory=dict)

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.chain}:{self.resname}{self.resseq}{ic}"

    def primary_bfactor(self) -> float | None:
        if "CA" in self.atoms:
            return self.atoms["CA"].bfactor
        if self.atoms:
            vals = [atom.bfactor for atom in self.atoms.values()]
            return float(sum(vals) / len(vals))
        return None


@dataclass
class ParsedStructure:
    """Comprehensive parsed representation of a PDB or mmCIF structure."""

    raw_text: str
    file_format: Literal["pdb", "cif"]
    filename: str
    structure_origin: str
    is_predicted: bool
    pdb_id: str | None
    title: str | None
    experimental_method: str | None
    resolution_angstrom: float | None
    model_count: int
    active_model_index: int
    chains: list[str]
    polymer_residues: dict[str, list[ResidueRecord]]
    polymer_residue_count: int
    polymer_atom_count: int
    hetero_ligand_count: int
    hetero_ligand_names: list[str]
    water_count: int
    chain_sequences: dict[str, str]
    missing_residues: list[dict[str, object]]
    missing_residue_source: str
    missing_residue_note: str
    header_ss_ranges: dict[str, list[tuple[int, int, str]]]
    bfactor_stats: dict[str, object] | None
    plddt_stats: dict[str, object] | None
    file_size_kib: float
    modified_residues_mapped: dict[str, list[str]] = field(default_factory=dict)
    unmapped_polymer_residues: dict[str, list[str]] = field(default_factory=dict)
    chain_extraction_notes: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def detect_structure_format(raw_text: str, filename: str = "") -> Literal["pdb", "cif"]:
    """Detect whether structure text is PDB or mmCIF."""
    lower = filename.lower().strip()
    if lower.endswith((".cif", ".mmcif")):
        return "cif"
    if lower.endswith((".pdb", ".ent")):
        return "pdb"

    for line in raw_text.splitlines()[:80]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("data_") or stripped.startswith("_atom_site."):
            return "cif"
        if stripped.startswith(("HEADER", "TITLE", "ATOM  ", "HETATM", "MODEL ")):
            return "pdb"
    return "pdb"


def _infer_origin(
    header_blob: str,
    experimental_method: str | None,
    origin_hint: str | None,
    has_ma_model: bool = False,
) -> tuple[str, bool]:
    """Classify structure origin without fabricating confidence provenance."""
    if origin_hint in {"experimental", "predicted_alphafold", "predicted_esmfold"}:
        return origin_hint, origin_hint != "experimental"

    upper = header_blob.upper()
    if (
        has_ma_model
        or "ALPHAFOLD" in upper
        or "COLABFOLD" in upper
        or "ESMFOLD" in upper
        or "ESM ATLAS" in upper
    ):
        return "uploaded_predicted", True

    if experimental_method:
        return "uploaded_experimental", False

    return "uploaded_unknown", False


def _choose_altloc(existing: AtomRecord | None, candidate: AtomRecord) -> AtomRecord:
    """Select the primary conformer atom by highest occupancy then standard altLoc order."""
    if existing is None:
        return candidate
    if candidate.occupancy > existing.occupancy + 1e-4:
        return candidate
    if abs(candidate.occupancy - existing.occupancy) <= 1e-4:
        rank = {"": 0, ".": 0, " ": 0, "A": 1, "B": 2}
        if rank.get(candidate.altloc, 9) < rank.get(existing.altloc, 9):
            return candidate
    return existing


def _parse_pdb_records(raw_text: str) -> dict[str, object]:
    """Parse PDB metadata, Model-0 coordinates, REMARK 465 missing residues, and HELIX/SHEET."""
    pdb_id: str | None = None
    title_parts: list[str] = []
    expdta_parts: list[str] = []
    resolution: float | None = None
    missing_residues: list[dict[str, object]] = []
    header_ss: dict[str, list[tuple[int, int, str]]] = {}
    header_lines_blob: list[str] = []

    model_count = 0
    in_model = False
    finished_first_model = False

    # Ordered residue tracking: key = (is_hetero_bucket, chain, resseq, icode, resname)
    residue_order: list[tuple[str, str, int, str, str]] = []
    residue_map: dict[tuple[str, str, int, str, str], ResidueRecord] = {}

    for line in raw_text.splitlines():
        rec = line[:6].strip()
        if rec in {"HEADER", "TITLE", "COMPND", "EXPDTA", "REMARK", "SOURCE"}:
            header_lines_blob.append(line)

        if rec == "HEADER" and len(line) >= 66:
            cand = line[62:66].strip()
            if len(cand) == 4 and cand.isalnum():
                pdb_id = cand.upper()
        elif rec == "TITLE":
            chunk = line[10:].strip()
            if chunk:
                title_parts.append(chunk)
        elif rec == "EXPDTA":
            chunk = line[10:].strip()
            if chunk:
                expdta_parts.append(chunk)
        elif line.startswith("REMARK   2 RESOLUTION."):
            m = re.search(r"RESOLUTION\.\s+([0-9]+(?:\.[0-9]+)?)\s+ANGSTROM", line, re.IGNORECASE)
            if m and resolution is None:
                resolution = float(m.group(1))
        elif line.startswith("REMARK 465"):
            # Standard PDB REMARK 465: col 15-18 resname, col 19 chain, col 21-26 resseq
            if len(line) >= 26:
                resname_cand = line[15:18].strip().upper()
                chain_cand = line[19:20].strip() or "A"
                seq_cand = line[21:26].strip()
                if (
                    resname_cand in STANDARD_RESIDUE_NAMES
                    and seq_cand.lstrip("-+").isdigit()
                    and "IDENTIFIED" not in line
                ):
                    model_col = line[12:14].strip()
                    if not model_col or model_col == "1":
                        missing_residues.append(
                            {
                                "chain": chain_cand,
                                "resname": resname_cand,
                                "resseq": int(seq_cand),
                                "source": "REMARK 465",
                            }
                        )
        elif rec == "HELIX" and len(line) >= 38:
            chain_id = line[19:20].strip() or "A"
            s_seq = line[21:25].strip()
            e_seq = line[33:37].strip()
            if s_seq.lstrip("-").isdigit() and e_seq.lstrip("-").isdigit():
                header_ss.setdefault(chain_id, []).append((int(s_seq), int(e_seq), "H"))
        elif rec == "SHEET" and len(line) >= 38:
            chain_id = line[21:22].strip() or "A"
            s_seq = line[22:26].strip()
            e_seq = line[33:37].strip()
            if s_seq.lstrip("-").isdigit() and e_seq.lstrip("-").isdigit():
                header_ss.setdefault(chain_id, []).append((int(s_seq), int(e_seq), "E"))
        elif rec == "MODEL":
            model_count += 1
            if model_count == 1:
                in_model = True
        elif rec == "ENDMDL":
            if in_model:
                in_model = False
                finished_first_model = True
        elif rec in {"ATOM", "HETATM"}:
            if finished_first_model:
                continue
            if len(line) < 54:
                continue
            atom_name = line[12:16].strip()
            altloc = line[16:17].strip()
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip() or "A"
            resseq_str = line[22:26].strip()
            if not resseq_str.lstrip("-+").isdigit():
                continue
            resseq = int(resseq_str)
            icode = line[26:27].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            try:
                occupancy = float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else 1.0
            except ValueError:
                occupancy = 1.0
            try:
                bfactor = float(line[60:66]) if len(line) >= 66 and line[60:66].strip() else 0.0
            except ValueError:
                bfactor = 0.0
            element = line[76:78].strip().upper() if len(line) >= 78 else ""

            is_water = resname in WATER_RESIDUE_NAMES
            if is_water:
                bucket = "water"
            elif rec == "ATOM":
                bucket = "poly" if resname in AA3_TO_AA1 else "unmapped_poly"
            elif resname in MODIFIED_AA3_TO_AA1:
                bucket = "mod_het"
            else:
                bucket = "het"
            is_polymer = bucket == "poly"
            key = (bucket, chain, resseq, icode, resname)

            if key not in residue_map:
                residue_order.append(key)
                residue_map[key] = ResidueRecord(
                    chain=chain,
                    resseq=resseq,
                    icode=icode,
                    resname=resname,
                    aa1=AA3_TO_AA1.get(resname, "X"),
                    is_hetero=not is_polymer,
                    is_water=is_water,
                )

            atom = AtomRecord(
                name=atom_name,
                coord=np.array([x, y, z], dtype=np.float64),
                bfactor=bfactor,
                occupancy=occupancy,
                altloc=altloc,
                element=element,
            )
            res_obj = residue_map[key]
            res_obj.atoms[atom_name] = _choose_altloc(res_obj.atoms.get(atom_name), atom)

    title = re.sub(r"\s+", " ", " ".join(title_parts)).strip() or None
    experimental_method = re.sub(r"\s+", " ", " ".join(expdta_parts)).strip() or None

    return {
        "pdb_id": pdb_id,
        "title": title,
        "experimental_method": experimental_method,
        "resolution": resolution,
        "model_count": max(1, model_count),
        "residue_order": residue_order,
        "residue_map": residue_map,
        "missing_residues": missing_residues,
        "missing_source": "REMARK 465" if missing_residues else "none",
        "header_ss": header_ss,
        "header_blob": "\n".join(header_lines_blob),
        "has_ma_model": False,
    }


def _tokenize_cif(raw_text: str) -> list[str]:
    """Tokenize an mmCIF file handling quoted strings and semicolon multiline values."""
    tokens: list[str] = []
    lines = raw_text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.startswith(";"):
            block = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(";"):
                block.append(lines[i])
                i += 1
            tokens.append("\n".join(block).strip())
            i += 1
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            i += 1
            continue
        # Remove trailing inline comments carefully outside quotes via shlex fallback
        try:
            parts = shlex.split(stripped, comments=True)
        except ValueError:
            parts = stripped.split()
        tokens.extend(parts)
        i += 1
    return tokens


def _parse_cif_tables(raw_text: str) -> tuple[dict[str, str], list[tuple[list[str], list[list[str]]]]]:
    """Extract scalar key-values and loop tables from mmCIF tokens."""
    tokens = _tokenize_cif(raw_text)
    scalars: dict[str, str] = {}
    loops: list[tuple[list[str], list[list[str]]]] = []
    i = 0
    n = len(tokens)
    while i < n:
        tok = tokens[i]
        if tok.lower() == "loop_":
            i += 1
            cols: list[str] = []
            while i < n and tokens[i].startswith("_"):
                cols.append(tokens[i])
                i += 1
            rows: list[list[str]] = []
            ncols = len(cols)
            if ncols > 0:
                while i + ncols <= n:
                    if tokens[i].lower() == "loop_" or tokens[i].startswith("_") or tokens[i].startswith("data_"):
                        break
                    rows.append(tokens[i : i + ncols])
                    i += ncols
                loops.append((cols, rows))
            continue
        if tok.startswith("_") and i + 1 < n and not tokens[i + 1].startswith("_"):
            scalars[tok] = tokens[i + 1]
            i += 2
            continue
        i += 1
    return scalars, loops


def _cif_float(val: str | None, default: float | None = None) -> float | None:
    if val is None or val in {".", "?"}:
        return default
    try:
        return float(val)
    except ValueError:
        return default


def _parse_cif_records(raw_text: str) -> dict[str, object]:
    """Parse mmCIF metadata, Model-0 coordinates, unobs residues, and secondary structure."""
    scalars, loops = _parse_cif_tables(raw_text)

    pdb_id = scalars.get("_entry.id")
    if pdb_id in {None, ".", "?"}:
        pdb_id = None
    title = scalars.get("_struct.title")
    if title in {None, ".", "?"}:
        title = None
    experimental_method = scalars.get("_exptl.method")
    if experimental_method in {None, ".", "?"}:
        experimental_method = None
    resolution = _cif_float(scalars.get("_refine.ls_d_res_high")) or _cif_float(
        scalars.get("_em_3d_reconstruction.resolution")
    )

    has_ma_model = any(k.startswith("_ma_model_list.") for k in scalars)
    header_blob_parts = [v for k, v in scalars.items() if any(t in k for t in ("title", "method", "software", "ma_model"))]

    residue_order: list[tuple[str, str, int, str, str]] = []
    residue_map: dict[tuple[str, str, int, str, str], ResidueRecord] = {}
    missing_residues: list[dict[str, object]] = []
    header_ss: dict[str, list[tuple[int, int, str]]] = {}
    model_ids: set[str] = set()
    first_model_id: str | None = None

    for cols, rows in loops:
        col_idx = {c: idx for idx, c in enumerate(cols)}

        if any(c.startswith("_exptl.") for c in cols) and experimental_method is None:
            m_idx = col_idx.get("_exptl.method")
            if m_idx is not None and rows:
                methods = [r[m_idx] for r in rows if r[m_idx] not in {".", "?"}]
                if methods:
                    experimental_method = "; ".join(methods)

        if any(c.startswith("_software.") or c.startswith("_ma_model_list.") for c in cols):
            if any(c.startswith("_ma_model_list.") for c in cols):
                has_ma_model = True
            for r in rows[:10]:
                header_blob_parts.extend(r)

        if any(c.startswith("_pdbx_unobs_or_zero_occ_residues.") for c in cols):
            c_chain = col_idx.get("_pdbx_unobs_or_zero_occ_residues.auth_asym_id") or col_idx.get(
                "_pdbx_unobs_or_zero_occ_residues.label_asym_id"
            )
            c_comp = col_idx.get("_pdbx_unobs_or_zero_occ_residues.auth_comp_id") or col_idx.get(
                "_pdbx_unobs_or_zero_occ_residues.label_comp_id"
            )
            c_seq = col_idx.get("_pdbx_unobs_or_zero_occ_residues.auth_seq_id") or col_idx.get(
                "_pdbx_unobs_or_zero_occ_residues.label_seq_id"
            )
            c_model = col_idx.get("_pdbx_unobs_or_zero_occ_residues.PDB_model_num")
            if c_chain is not None and c_comp is not None and c_seq is not None:
                for r in rows:
                    if c_model is not None and r[c_model] not in {"1", ".", "?"}:
                        continue
                    resname = r[c_comp].strip().upper()
                    if resname in STANDARD_RESIDUE_NAMES and r[c_seq].lstrip("-+").isdigit():
                        missing_residues.append(
                            {
                                "chain": r[c_chain].strip() or "A",
                                "resname": resname,
                                "resseq": int(r[c_seq]),
                                "source": "mmCIF _pdbx_unobs_or_zero_occ_residues",
                            }
                        )

        if any(c.startswith("_struct_conf.") for c in cols):
            c_type = col_idx.get("_struct_conf.conf_type_id")
            c_chain = col_idx.get("_struct_conf.beg_auth_asym_id") or col_idx.get("_struct_conf.beg_label_asym_id")
            c_beg = col_idx.get("_struct_conf.beg_auth_seq_id") or col_idx.get("_struct_conf.beg_label_seq_id")
            c_end = col_idx.get("_struct_conf.end_auth_seq_id") or col_idx.get("_struct_conf.end_label_seq_id")
            if c_type is not None and c_chain is not None and c_beg is not None and c_end is not None:
                for r in rows:
                    if "HELX" in r[c_type].upper() and r[c_beg].lstrip("-").isdigit() and r[c_end].lstrip("-").isdigit():
                        header_ss.setdefault(r[c_chain].strip() or "A", []).append(
                            (int(r[c_beg]), int(r[c_end]), "H")
                        )

        if any(c.startswith("_struct_sheet_range.") for c in cols):
            c_chain = col_idx.get("_struct_sheet_range.beg_auth_asym_id") or col_idx.get(
                "_struct_sheet_range.beg_label_asym_id"
            )
            c_beg = col_idx.get("_struct_sheet_range.beg_auth_seq_id") or col_idx.get(
                "_struct_sheet_range.beg_label_seq_id"
            )
            c_end = col_idx.get("_struct_sheet_range.end_auth_seq_id") or col_idx.get(
                "_struct_sheet_range.end_label_seq_id"
            )
            if c_chain is not None and c_beg is not None and c_end is not None:
                for r in rows:
                    if r[c_beg].lstrip("-").isdigit() and r[c_end].lstrip("-").isdigit():
                        header_ss.setdefault(r[c_chain].strip() or "A", []).append(
                            (int(r[c_beg]), int(r[c_end]), "E")
                        )

        if any(c.startswith("_atom_site.") for c in cols):
            c_group = col_idx.get("_atom_site.group_PDB")
            c_atom = col_idx.get("_atom_site.auth_atom_id") or col_idx.get("_atom_site.label_atom_id")
            c_alt = col_idx.get("_atom_site.label_alt_id")
            c_comp = col_idx.get("_atom_site.auth_comp_id") or col_idx.get("_atom_site.label_comp_id")
            c_chain = col_idx.get("_atom_site.auth_asym_id") or col_idx.get("_atom_site.label_asym_id")
            c_seq = col_idx.get("_atom_site.auth_seq_id") or col_idx.get("_atom_site.label_seq_id")
            c_ins = col_idx.get("_atom_site.pdbx_PDB_ins_code")
            c_x = col_idx.get("_atom_site.Cartn_x")
            c_y = col_idx.get("_atom_site.Cartn_y")
            c_z = col_idx.get("_atom_site.Cartn_z")
            c_occ = col_idx.get("_atom_site.occupancy")
            c_b = col_idx.get("_atom_site.B_iso_or_equiv")
            c_elem = col_idx.get("_atom_site.type_symbol")
            c_model = col_idx.get("_atom_site.pdbx_PDB_model_num")

            if None in (c_atom, c_comp, c_chain, c_seq, c_x, c_y, c_z):
                continue

            for r in rows:
                m_num = r[c_model] if c_model is not None else "1"
                model_ids.add(m_num)
                if first_model_id is None:
                    first_model_id = m_num
                if m_num != first_model_id:
                    continue

                rec = r[c_group].upper() if c_group is not None else "ATOM"
                resname = r[c_comp].strip().upper()
                chain = r[c_chain].strip() or "A"
                seq_str = r[c_seq].strip()
                if not seq_str.lstrip("-+").isdigit():
                    continue
                resseq = int(seq_str)
                icode = "" if c_ins is None or r[c_ins] in {".", "?"} else r[c_ins].strip()
                altloc = "" if c_alt is None or r[c_alt] in {".", "?"} else r[c_alt].strip()
                atom_name = r[c_atom].strip()
                x = _cif_float(r[c_x])
                y = _cif_float(r[c_y])
                z = _cif_float(r[c_z])
                if x is None or y is None or z is None:
                    continue
                occ = _cif_float(r[c_occ] if c_occ is not None else None, 1.0) or 1.0
                bfac = _cif_float(r[c_b] if c_b is not None else None, 0.0) or 0.0
                elem = r[c_elem].strip().upper() if c_elem is not None and r[c_elem] not in {".", "?"} else ""

                is_water = resname in WATER_RESIDUE_NAMES
                if is_water:
                    bucket = "water"
                elif rec == "ATOM":
                    bucket = "poly" if resname in AA3_TO_AA1 else "unmapped_poly"
                elif resname in MODIFIED_AA3_TO_AA1:
                    bucket = "mod_het"
                else:
                    bucket = "het"
                is_polymer = bucket == "poly"
                key = (bucket, chain, resseq, icode, resname)
                if key not in residue_map:
                    residue_order.append(key)
                    residue_map[key] = ResidueRecord(
                        chain=chain,
                        resseq=resseq,
                        icode=icode,
                        resname=resname,
                        aa1=AA3_TO_AA1.get(resname, "X"),
                        is_hetero=not is_polymer,
                        is_water=is_water,
                    )

                atom = AtomRecord(
                    name=atom_name,
                    coord=np.array([x, y, z], dtype=np.float64),
                    bfactor=bfac,
                    occupancy=occ,
                    altloc=altloc,
                    element=elem,
                )
                res_obj = residue_map[key]
                res_obj.atoms[atom_name] = _choose_altloc(res_obj.atoms.get(atom_name), atom)

    return {
        "pdb_id": pdb_id,
        "title": title,
        "experimental_method": experimental_method,
        "resolution": resolution,
        "model_count": max(1, len(model_ids)),
        "residue_order": residue_order,
        "residue_map": residue_map,
        "missing_residues": missing_residues,
        "missing_source": "mmCIF _pdbx_unobs_or_zero_occ_residues" if missing_residues else "none",
        "header_ss": header_ss,
        "header_blob": "\n".join(header_blob_parts),
        "has_ma_model": has_ma_model,
    }


def _detect_coordinate_breaks(
    polymer_residues: dict[str, list[ResidueRecord]],
) -> list[dict[str, object]]:
    """Detect backbone numbering or peptide-bond breaks when formal missing records are absent."""
    breaks: list[dict[str, object]] = []
    for chain, residues in polymer_residues.items():
        for i in range(len(residues) - 1):
            r1, r2 = residues[i], residues[i + 1]
            seq_gap = r2.resseq - r1.resseq
            peptide_dist = None
            if "C" in r1.atoms and "N" in r2.atoms:
                peptide_dist = float(np.linalg.norm(r1.atoms["C"].coord - r2.atoms["N"].coord))
            if seq_gap > 1 or (peptide_dist is not None and peptide_dist > 2.2):
                breaks.append(
                    {
                        "chain": chain,
                        "resname": f"Gap {r1.resseq}->{r2.resseq}",
                        "resseq": r1.resseq + 1 if seq_gap > 1 else r1.resseq,
                        "source": "backbone_break_check",
                    }
                )
    return breaks


def parse_structure(
    raw_text: str,
    filename: str = "structure.pdb",
    origin_hint: str | None = None,
) -> ParsedStructure:
    """Parse a PDB or mmCIF structure with Model-0, altLoc, and provenance discipline."""
    if not raw_text or not raw_text.strip():
        raise ValueError("Structure input is empty.")

    fmt = detect_structure_format(raw_text, filename)
    raw_parsed = _parse_cif_records(raw_text) if fmt == "cif" else _parse_pdb_records(raw_text)

    residue_order: list[tuple[str, str, int, str, str]] = raw_parsed["residue_order"]  # type: ignore[assignment]
    residue_map: dict[tuple[str, str, int, str, str], ResidueRecord] = raw_parsed["residue_map"]  # type: ignore[assignment]

    polymer_residues: dict[str, list[ResidueRecord]] = {}
    hetero_ligands: list[ResidueRecord] = []
    water_count = 0
    modified_residues_mapped: dict[str, list[str]] = {}
    unmapped_polymer_residues: dict[str, list[str]] = {}

    for key in residue_order:
        bucket, chain, _, _, _ = key
        res = residue_map[key]
        if bucket == "poly":
            polymer_residues.setdefault(chain, []).append(res)
            if res.resname in MODIFIED_AA3_TO_AA1:
                modified_residues_mapped.setdefault(chain, []).append(
                    f"{res.resname}->{res.aa1} ({res.resseq}{res.icode})"
                )
        elif bucket == "mod_het":
            # Promote HETATM modified amino acids only when integrated with peptide backbone atoms
            if {"N", "CA", "C"}.issubset(res.atoms) or (res.resname == "MSE" and "CA" in res.atoms):
                res.is_hetero = False
                polymer_residues.setdefault(chain, []).append(res)
                modified_residues_mapped.setdefault(chain, []).append(
                    f"{res.resname}->{res.aa1} ({res.resseq}{res.icode})"
                )
            else:
                hetero_ligands.append(res)
        elif bucket == "unmapped_poly":
            if "CA" in res.atoms or {"N", "C"}.issubset(res.atoms) or res.resname == "UNK":
                polymer_residues.setdefault(chain, []).append(res)
                unmapped_polymer_residues.setdefault(chain, []).append(
                    f"{res.resname} {res.resseq}{res.icode}"
                )
            else:
                hetero_ligands.append(res)
        elif bucket == "water":
            water_count += 1
        else:
            hetero_ligands.append(res)

    if not polymer_residues:
        raise ValueError("No valid standard protein polymer residues were found in the structure.")

    chains = list(polymer_residues.keys())
    polymer_residue_count = sum(len(v) for v in polymer_residues.values())
    polymer_atom_count = sum(len(r.atoms) for chain_res in polymer_residues.values() for r in chain_res)
    chain_sequences = {
        chain: "".join(r.aa1 for r in res_list) for chain, res_list in polymer_residues.items()
    }
    chain_extraction_notes: dict[str, str] = {}
    for ch in chains:
        unmapped = unmapped_polymer_residues.get(ch, [])
        mods = modified_residues_mapped.get(ch, [])
        res_len = len(polymer_residues[ch])
        if unmapped:
            preview = ", ".join(unmapped[:5]) + ("..." if len(unmapped) > 5 else "")
            chain_extraction_notes[ch] = (
                f"Chain {ch} contains {len(unmapped)} unmapped polymer residue(s) ({preview}); "
                "standard 20-amino-acid physicochemical properties cannot be calculated without fabricating residue identity."
            )
        elif mods:
            preview = ", ".join(mods[:5]) + ("..." if len(mods) > 5 else "")
            chain_extraction_notes[ch] = (
                f"Chain {ch}: {res_len} aa extracted with {len(mods)} modified residue(s) mapped to canonical parent ({preview})."
            )
        else:
            chain_extraction_notes[ch] = (
                f"Chain {ch}: {res_len} standard amino-acid residues extracted cleanly from Model 1 coordinates."
            )
    hetero_names = sorted({r.resname for r in hetero_ligands})

    origin, is_predicted = _infer_origin(
        str(raw_parsed["header_blob"]),
        raw_parsed["experimental_method"],  # type: ignore[arg-type]
        origin_hint,
        has_ma_model=bool(raw_parsed["has_ma_model"]),
    )

    # Strictly separate experimental B-factor vs. computational pLDDT
    per_res_scores: list[dict[str, object]] = []
    raw_b_values: list[float] = []
    for chain in chains:
        for res in polymer_residues[chain]:
            val = res.primary_bfactor()
            if val is not None:
                raw_b_values.append(val)
                per_res_scores.append(
                    {
                        "chain": chain,
                        "resseq": res.resseq,
                        "resname": res.resname,
                        "aa1": res.aa1,
                        "value": float(val),
                    }
                )

    bfactor_stats: dict[str, object] | None = None
    plddt_stats: dict[str, object] | None = None
    if raw_b_values:
        arr = np.array(raw_b_values, dtype=np.float64)
        if is_predicted:
            max_raw = float(np.max(arr))
            scale_factor = 100.0 if max_raw <= 1.05 else 1.0
            raw_scale_label = "0-1 (normalized x100)" if scale_factor == 100.0 else "0-100"
            scaled = arr * scale_factor
            for entry in per_res_scores:
                entry["value"] = float(entry["value"]) * scale_factor
            n_vals = len(scaled)
            plddt_stats = {
                "raw_scale": raw_scale_label,
                "mean": float(np.mean(scaled)),
                "median": float(np.median(scaled)),
                "min": float(np.min(scaled)),
                "max": float(np.max(scaled)),
                "pct_very_high": float(np.sum(scaled >= 90.0) / n_vals * 100.0),
                "pct_confident": float(np.sum((scaled >= 70.0) & (scaled < 90.0)) / n_vals * 100.0),
                "pct_low": float(np.sum((scaled >= 50.0) & (scaled < 70.0)) / n_vals * 100.0),
                "pct_very_low": float(np.sum(scaled < 50.0) / n_vals * 100.0),
                "per_residue": per_res_scores,
            }
        else:
            all_atom_b = [
                atom.bfactor
                for ch in chains
                for r in polymer_residues[ch]
                for atom in r.atoms.values()
            ]
            all_atom_occ = [
                atom.occupancy
                for ch in chains
                for r in polymer_residues[ch]
                for atom in r.atoms.values()
            ]
            b_all_arr = np.array(all_atom_b, dtype=np.float64) if all_atom_b else arr
            occ_arr = np.array(all_atom_occ, dtype=np.float64) if all_atom_occ else np.array([1.0])
            bfactor_stats = {
                "min": float(np.min(arr)),
                "p5": float(np.percentile(arr, 5)),
                "median": float(np.median(arr)),
                "mean": float(np.mean(arr)),
                "p95": float(np.percentile(arr, 95)),
                "max": float(np.max(arr)),
                "all_atom_count": int(len(b_all_arr)),
                "all_atom_min": float(np.min(b_all_arr)),
                "all_atom_median": float(np.median(b_all_arr)),
                "all_atom_mean": float(np.mean(b_all_arr)),
                "all_atom_max": float(np.max(b_all_arr)),
                "occupancy_min": float(np.min(occ_arr)),
                "occupancy_max": float(np.max(occ_arr)),
                "occupancy_mean": float(np.mean(occ_arr)),
                "per_residue": per_res_scores,
            }

    missing_residues: list[dict[str, object]] = raw_parsed["missing_residues"]  # type: ignore[assignment]
    missing_source = str(raw_parsed["missing_source"])
    if missing_residues:
        missing_note = (
            f"Found {len(missing_residues)} missing residue record(s) from {missing_source} (Model 1 / primary)."
        )
    else:
        breaks = _detect_coordinate_breaks(polymer_residues)
        if breaks:
            missing_residues = breaks
            missing_source = "backbone_break_check"
            missing_note = (
                f"No formal missing-residue records found; detected {len(breaks)} internal residue "
                "numbering or peptide-bond break(s) from coordinates."
            )
        else:
            missing_note = (
                "No formal missing-residue records (REMARK 465 / _pdbx_unobs_or_zero_occ_residues) "
                "and no internal backbone numbering breaks detected."
            )

    warnings: list[str] = []
    model_count = int(raw_parsed["model_count"])  # type: ignore[arg-type]
    if model_count > 1:
        warnings.append(
            f"Multi-model file ({model_count} models detected); structural analytics use Model 1 only."
        )

    return ParsedStructure(
        raw_text=raw_text,
        file_format=fmt,
        filename=filename,
        structure_origin=origin,
        is_predicted=is_predicted,
        pdb_id=raw_parsed["pdb_id"],  # type: ignore[arg-type]
        title=raw_parsed["title"],  # type: ignore[arg-type]
        experimental_method=raw_parsed["experimental_method"],  # type: ignore[arg-type]
        resolution_angstrom=raw_parsed["resolution"],  # type: ignore[arg-type]
        model_count=model_count,
        active_model_index=0,
        chains=chains,
        polymer_residues=polymer_residues,
        polymer_residue_count=polymer_residue_count,
        polymer_atom_count=polymer_atom_count,
        hetero_ligand_count=len(hetero_ligands),
        hetero_ligand_names=hetero_names,
        water_count=water_count,
        chain_sequences=chain_sequences,
        missing_residues=missing_residues,
        missing_residue_source=missing_source,
        missing_residue_note=missing_note,
        header_ss_ranges=raw_parsed["header_ss"],  # type: ignore[arg-type]
        bfactor_stats=bfactor_stats,
        plddt_stats=plddt_stats,
        file_size_kib=len(raw_text.encode("utf-8")) / 1024.0,
        modified_residues_mapped=modified_residues_mapped,
        unmapped_polymer_residues=unmapped_polymer_residues,
        chain_extraction_notes=chain_extraction_notes,
        warnings=warnings,
    )
