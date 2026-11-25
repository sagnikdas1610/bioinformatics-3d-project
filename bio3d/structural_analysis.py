"""Coordinate-based structural analysis: MolProbity Top8000 Ramachandran, Kabsch-Sander SS, and Contact Maps."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

import numpy as np

from bio3d.config import DATA_DIR
from bio3d.structure import ParsedStructure, ResidueRecord

PEPTIDE_BOND_MAX_DIST = 2.0  # Angstroms (canonical C-N ~1.33 A)
MOLPROBITY_RAMA_NPZ = DATA_DIR / "molprobity_top8000_rama.npz"
MOLPROBITY_REFERENCE = (
    "Richardson Lab MolProbity Top8000 (CC BY 4.0; Williams et al., Protein Sci. 2018; "
    "Lovell et al., Proteins 2003)"
)
RamaClass = Literal[
    "General",
    "Ile/Val",
    "Glycine",
    "trans-Proline",
    "cis-Proline",
    "Pre-Proline",
]
RAMA_CLASSES: tuple[RamaClass, ...] = (
    "General",
    "Ile/Val",
    "Glycine",
    "trans-Proline",
    "cis-Proline",
    "Pre-Proline",
)


@dataclass(frozen=True)
class RamachandranResidue:
    chain: str
    resseq: int
    icode: str
    resname: str
    aa1: str
    rama_class: RamaClass
    phi: float
    psi: float
    omega: float
    prob: float
    region: Literal["Favored", "Allowed", "Outlier"]


@dataclass(frozen=True)
class RamachandranResult:
    chain_filter: str
    reference_source: str
    residues: list[RamachandranResidue]
    total_evaluated: int
    excluded_count: int
    favored_count: int
    allowed_count: int
    outlier_count: int
    favored_pct: float
    allowed_pct: float
    outlier_pct: float
    outliers: list[RamachandranResidue]
    interpretation: str


@dataclass(frozen=True)
class SecondaryStructureResidue:
    chain: str
    resseq: int
    resname: str
    aa1: str
    dssp_code: str  # H, G, I, E, B, T, C
    state3: Literal["Helix", "Sheet", "Coil/Turn"]
    header_state3: str | None


@dataclass(frozen=True)
class SecondaryStructureResult:
    chain_filter: str
    method: str
    residues: list[SecondaryStructureResidue]
    total_residues: int
    helix_count: int
    sheet_count: int
    coil_count: int
    helix_pct: float
    sheet_pct: float
    coil_pct: float
    state6_counts: dict[str, int]
    header_available: bool
    header_helix_pct: float | None
    header_sheet_pct: float | None
    header_coil_pct: float | None


@dataclass(frozen=True)
class ContactMapResult:
    chain_id: str
    atom_mode: Literal["CA", "CB"]
    threshold_angstrom: float
    min_seq_sep: int
    residue_count: int
    truncated_from: int | None
    glycine_fallback_count: int
    missing_cb_fallback_count: int
    residue_labels: list[str]
    residue_numbers: list[int]
    distance_matrix: np.ndarray
    contact_matrix: np.ndarray
    total_contacts: int
    long_range_contacts: int
    contact_density_pct: float
    top_contacts: list[dict[str, object]]


def _calc_dihedral_deg(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Compute signed dihedral angle in degrees [-180, 180] from 4 3D points."""
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1_norm = np.linalg.norm(b1)
    if b1_norm < 1e-8:
        return 0.0
    b1 = b1 / b1_norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


@lru_cache(maxsize=1)
def _load_top8000_tables() -> dict[str, np.ndarray]:
    """Load packed Richardson Lab MolProbity Top8000 2-deg probability grids."""
    if not MOLPROBITY_RAMA_NPZ.exists():
        raise FileNotFoundError(f"MolProbity Top8000 reference file missing: {MOLPROBITY_RAMA_NPZ}")
    with np.load(MOLPROBITY_RAMA_NPZ) as data:
        return {key: data[key].astype(np.float64) for key in data.files}


def top8000_probability(phi: float, psi: float, rama_class: RamaClass) -> float:
    """Evaluate MolProbity Top8000 contour probability using periodic bilinear interpolation."""
    tables = _load_top8000_tables()
    grid = tables[rama_class]  # shape (180, 180) indexed as [psi_bin, phi_bin], centers -179..+179

    ux = ((phi + 179.0) / 2.0) % 180.0
    uy = ((psi + 179.0) / 2.0) % 180.0
    ix0 = int(np.floor(ux)) % 180
    iy0 = int(np.floor(uy)) % 180
    ix1 = (ix0 + 1) % 180
    iy1 = (iy0 + 1) % 180
    fx = ux - np.floor(ux)
    fy = uy - np.floor(uy)

    return float(
        (1.0 - fx) * (1.0 - fy) * grid[iy0, ix0]
        + fx * (1.0 - fy) * grid[iy0, ix1]
        + (1.0 - fx) * fy * grid[iy1, ix0]
        + fx * fy * grid[iy1, ix1]
    )


def classify_ramachandran(
    phi: float,
    psi: float,
    rama_class: RamaClass,
) -> tuple[Literal["Favored", "Allowed", "Outlier"], float]:
    """Classify (phi, psi) using official MolProbity Top8000 cutoffs (0.02 favored, 0.0005 allowed)."""
    prob = top8000_probability(phi, psi, rama_class)
    if prob >= 0.02:
        return "Favored", prob
    if prob >= 0.0005:
        return "Allowed", prob
    return "Outlier", prob


def rama_background_grid(
    rama_class: RamaClass,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return exact Top8000 (phi, psi, z) grid (0=Outlier, 1=Allowed, 2=Favored) at 2-deg resolution."""
    tables = _load_top8000_tables()
    grid = tables[rama_class]
    centers = np.arange(-179.0, 180.0, 2.0, dtype=np.float64)
    z = np.zeros_like(grid, dtype=np.float64)
    z[grid >= 0.0005] = 1.0
    z[grid >= 0.02] = 2.0
    return centers, centers, z


def compute_ramachandran(
    parsed: ParsedStructure,
    chain_id: str | None = None,
) -> RamachandranResult:
    """Compute backbone phi/psi dihedral angles and MolProbity Top8000 6-class validation."""
    target_chains = [chain_id] if chain_id and chain_id != "All" else parsed.chains
    results: list[RamachandranResidue] = []
    excluded = 0

    for ch in target_chains:
        residues = parsed.polymer_residues.get(ch, [])
        n = len(residues)
        for i, res in enumerate(residues):
            if i == 0 or i == n - 1:
                excluded += 1
                continue
            prev_res = residues[i - 1]
            next_res = residues[i + 1]

            required_curr = ("N", "CA", "C")
            if not all(a in res.atoms for a in required_curr):
                excluded += 1
                continue
            if "C" not in prev_res.atoms or "N" not in next_res.atoms:
                excluded += 1
                continue

            d_prev = float(np.linalg.norm(prev_res.atoms["C"].coord - res.atoms["N"].coord))
            d_next = float(np.linalg.norm(res.atoms["C"].coord - next_res.atoms["N"].coord))
            if d_prev > PEPTIDE_BOND_MAX_DIST or d_next > PEPTIDE_BOND_MAX_DIST:
                excluded += 1
                continue

            phi = _calc_dihedral_deg(
                prev_res.atoms["C"].coord,
                res.atoms["N"].coord,
                res.atoms["CA"].coord,
                res.atoms["C"].coord,
            )
            psi = _calc_dihedral_deg(
                res.atoms["N"].coord,
                res.atoms["CA"].coord,
                res.atoms["C"].coord,
                next_res.atoms["N"].coord,
            )
            omega = (
                _calc_dihedral_deg(
                    prev_res.atoms["CA"].coord,
                    prev_res.atoms["C"].coord,
                    res.atoms["N"].coord,
                    res.atoms["CA"].coord,
                )
                if "CA" in prev_res.atoms
                else 180.0
            )

            if res.resname == "GLY":
                rama_class: RamaClass = "Glycine"
            elif res.resname == "PRO":
                rama_class = "cis-Proline" if abs(omega) < 90.0 else "trans-Proline"
            elif next_res.resname == "PRO":
                rama_class = "Pre-Proline"
            elif res.resname in {"ILE", "VAL"}:
                rama_class = "Ile/Val"
            else:
                rama_class = "General"

            region, prob = classify_ramachandran(phi, psi, rama_class)
            results.append(
                RamachandranResidue(
                    chain=ch,
                    resseq=res.resseq,
                    icode=res.icode,
                    resname=res.resname,
                    aa1=res.aa1,
                    rama_class=rama_class,
                    phi=phi,
                    psi=psi,
                    omega=omega,
                    prob=prob,
                    region=region,
                )
            )

    total = len(results)
    favored = sum(1 for r in results if r.region == "Favored")
    allowed = sum(1 for r in results if r.region == "Allowed")
    outliers = [r for r in results if r.region == "Outlier"]
    outlier_cnt = len(outliers)

    fav_pct = (favored / total * 100.0) if total else 0.0
    all_pct = (allowed / total * 100.0) if total else 0.0
    out_pct = (outlier_cnt / total * 100.0) if total else 0.0

    if total == 0:
        interp = "No complete internal peptide triads were available to compute phi/psi torsion angles."
    elif fav_pct >= 98.0 and out_pct <= 0.2:
        interp = (
            f"High-quality backbone stereochemistry ({fav_pct:.1f}% favored, {out_pct:.1f}% outliers vs. MolProbity Top8000), "
            "meeting standard high-resolution validation targets (>=98% favored)."
        )
    elif fav_pct >= 95.0 and out_pct <= 1.5:
        interp = (
            f"Strong backbone stereochemistry ({fav_pct:.1f}% favored, {out_pct:.1f}% outliers vs. MolProbity Top8000), "
            "consistent with well-refined experimental or high-accuracy predicted geometry."
        )
    elif fav_pct >= 90.0:
        interp = (
            f"Acceptable backbone geometry ({fav_pct:.1f}% favored, {out_pct:.1f}% outliers vs. MolProbity Top8000); "
            "inspect outlier residues for strained loops or modeling artifacts."
        )
    else:
        interp = (
            f"Elevated backbone strain ({fav_pct:.1f}% favored, {out_pct:.1f}% outliers vs. MolProbity Top8000); "
            "residues outside allowed basins warrant stereochemical refinement or inspection."
        )

    return RamachandranResult(
        chain_filter=chain_id or "All",
        reference_source=MOLPROBITY_REFERENCE,
        residues=results,
        total_evaluated=total,
        excluded_count=excluded,
        favored_count=favored,
        allowed_count=allowed,
        outlier_count=outlier_cnt,
        favored_pct=fav_pct,
        allowed_pct=all_pct,
        outlier_pct=out_pct,
        outliers=outliers,
        interpretation=interp,
    )


def _header_state_lookup(
    header_ranges: dict[str, list[tuple[int, int, str]]],
    chain: str,
    resseq: int,
) -> str | None:
    if chain not in header_ranges:
        return None
    for beg, end, code in header_ranges[chain]:
        if beg <= resseq <= end:
            return "Helix" if code == "H" else "Sheet"
    return "Coil/Turn"


def compute_secondary_structure(
    parsed: ParsedStructure,
    chain_id: str | None = None,
) -> SecondaryStructureResult:
    """Assign secondary structure via Kabsch-Sander (1983) electrostatic H-bond geometry."""
    target_chains = [chain_id] if chain_id and chain_id != "All" else parsed.chains

    flat_residues: list[ResidueRecord] = []
    chain_bounds: list[tuple[int, int]] = []
    for ch in target_chains:
        res_list = parsed.polymer_residues.get(ch, [])
        start = len(flat_residues)
        flat_residues.extend(res_list)
        chain_bounds.append((start, len(flat_residues)))

    n_res = len(flat_residues)
    if n_res == 0:
        return SecondaryStructureResult(
            chain_filter=chain_id or "All",
            method="Kabsch-Sander (1983) backbone H-bond DSSP",
            residues=[],
            total_residues=0,
            helix_count=0,
            sheet_count=0,
            coil_count=0,
            helix_pct=0.0,
            sheet_pct=0.0,
            coil_pct=0.0,
            state6_counts={},
            header_available=False,
            header_helix_pct=None,
            header_sheet_pct=None,
            header_coil_pct=None,
        )

    bonded_prev = [False] * n_res
    for start, end in chain_bounds:
        for i in range(start + 1, end):
            r_prev, r_curr = flat_residues[i - 1], flat_residues[i]
            if "C" in r_prev.atoms and "N" in r_curr.atoms:
                dist = float(np.linalg.norm(r_prev.atoms["C"].coord - r_curr.atoms["N"].coord))
                if dist <= PEPTIDE_BOND_MAX_DIST:
                    bonded_prev[i] = True

    h_coords: list[np.ndarray | None] = [None] * n_res
    for i, res in enumerate(flat_residues):
        if not bonded_prev[i] or res.resname == "PRO":
            continue
        r_prev = flat_residues[i - 1]
        if "N" in res.atoms and "C" in r_prev.atoms and "O" in r_prev.atoms:
            co_vec = r_prev.atoms["O"].coord - r_prev.atoms["C"].coord
            co_norm = np.linalg.norm(co_vec)
            if co_norm > 1e-6:
                h_coords[i] = res.atoms["N"].coord - 1.01 * (co_vec / co_norm)

    hb = np.zeros((n_res, n_res), dtype=bool)
    for j, acc in enumerate(flat_residues):
        if "C" not in acc.atoms or "O" not in acc.atoms:
            continue
        c_pos = acc.atoms["C"].coord
        o_pos = acc.atoms["O"].coord
        for i, don in enumerate(flat_residues):
            if abs(i - j) < 3 or h_coords[i] is None or "N" not in don.atoms:
                continue
            n_pos = don.atoms["N"].coord
            r_on = float(np.linalg.norm(o_pos - n_pos))
            if r_on >= 5.5:
                continue
            h_pos = h_coords[i]
            assert h_pos is not None
            r_ch = float(np.linalg.norm(c_pos - h_pos))
            r_oh = float(np.linalg.norm(o_pos - h_pos))
            r_cn = float(np.linalg.norm(c_pos - n_pos))
            if min(r_on, r_ch, r_oh, r_cn) < 0.5:
                continue
            energy = 27.888 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn)
            if energy < -0.5:
                hb[j, i] = True

    def same_segment(a: int, b: int) -> bool:
        lo, hi = min(a, b), max(a, b)
        return all(bonded_prev[k] for k in range(lo + 1, hi + 1))

    turn3 = np.zeros(n_res, dtype=bool)
    turn4 = np.zeros(n_res, dtype=bool)
    turn5 = np.zeros(n_res, dtype=bool)
    for i in range(n_res):
        if i + 3 < n_res and same_segment(i, i + 3) and hb[i, i + 3]:
            turn3[i] = True
        if i + 4 < n_res and same_segment(i, i + 4) and hb[i, i + 4]:
            turn4[i] = True
        if i + 5 < n_res and same_segment(i, i + 5) and hb[i, i + 5]:
            turn5[i] = True

    codes = ["C"] * n_res

    for i in range(n_res - 5):
        if turn4[i] and turn4[i + 1] and same_segment(i, i + 5):
            for k in range(i + 1, i + 5):
                codes[k] = "H"

    bridge_partner_count = [0] * n_res
    for i in range(1, n_res - 1):
        if not (bonded_prev[i] and bonded_prev[i + 1]):
            continue
        for j in range(1, n_res - 1):
            if abs(i - j) < 3 or not (bonded_prev[j] and bonded_prev[j + 1]):
                continue
            parallel = (hb[i - 1, j] and hb[j, i + 1]) or (hb[j - 1, i] and hb[i, j + 1])
            antiparallel = (hb[i, j] and hb[j, i]) or (hb[i - 1, j + 1] and hb[j - 1, i + 1])
            if parallel or antiparallel:
                bridge_partner_count[i] += 1
                bridge_partner_count[j] += 1

    for i in range(n_res):
        if bridge_partner_count[i] > 0 and codes[i] != "H":
            adj_bridge = (i > 0 and bonded_prev[i] and bridge_partner_count[i - 1] > 0) or (
                i + 1 < n_res and bonded_prev[i + 1] and bridge_partner_count[i + 1] > 0
            )
            codes[i] = "E" if adj_bridge else "B"

    for i in range(n_res - 4):
        if turn3[i] and turn3[i + 1] and same_segment(i, i + 4):
            for k in range(i + 1, i + 4):
                if codes[k] == "C":
                    codes[k] = "G"
    for i in range(n_res - 6):
        if turn5[i] and turn5[i + 1] and same_segment(i, i + 6):
            for k in range(i + 1, i + 6):
                if codes[k] == "C":
                    codes[k] = "I"

    for n_turn, arr in ((3, turn3), (4, turn4), (5, turn5)):
        for i in range(n_res - n_turn):
            if arr[i]:
                for k in range(i + 1, i + n_turn):
                    if codes[k] == "C":
                        codes[k] = "T"

    header_available = any(ch in parsed.header_ss_ranges for ch in target_chains)
    ss_residues: list[SecondaryStructureResidue] = []
    state6_counts: dict[str, int] = {"H": 0, "G": 0, "I": 0, "E": 0, "B": 0, "T": 0, "C": 0}
    header_counts = {"Helix": 0, "Sheet": 0, "Coil/Turn": 0}

    for res, code in zip(flat_residues, codes):
        state6_counts[code] = state6_counts.get(code, 0) + 1
        if code in {"H", "G", "I"}:
            s3: Literal["Helix", "Sheet", "Coil/Turn"] = "Helix"
        elif code in {"E", "B"}:
            s3 = "Sheet"
        else:
            s3 = "Coil/Turn"

        h_state = _header_state_lookup(parsed.header_ss_ranges, res.chain, res.resseq) if header_available else None
        if h_state:
            header_counts[h_state] += 1

        ss_residues.append(
            SecondaryStructureResidue(
                chain=res.chain,
                resseq=res.resseq,
                resname=res.resname,
                aa1=res.aa1,
                dssp_code=code,
                state3=s3,
                header_state3=h_state,
            )
        )

    helix_cnt = sum(1 for r in ss_residues if r.state3 == "Helix")
    sheet_cnt = sum(1 for r in ss_residues if r.state3 == "Sheet")
    coil_cnt = sum(1 for r in ss_residues if r.state3 == "Coil/Turn")

    h_tot = sum(header_counts.values())
    return SecondaryStructureResult(
        chain_filter=chain_id or "All",
        method="Kabsch-Sander (1983) backbone H-bond DSSP",
        residues=ss_residues,
        total_residues=n_res,
        helix_count=helix_cnt,
        sheet_count=sheet_cnt,
        coil_count=coil_cnt,
        helix_pct=helix_cnt / n_res * 100.0,
        sheet_pct=sheet_cnt / n_res * 100.0,
        coil_pct=coil_cnt / n_res * 100.0,
        state6_counts=state6_counts,
        header_available=header_available and h_tot > 0,
        header_helix_pct=(header_counts["Helix"] / h_tot * 100.0) if h_tot else None,
        header_sheet_pct=(header_counts["Sheet"] / h_tot * 100.0) if h_tot else None,
        header_coil_pct=(header_counts["Coil/Turn"] / h_tot * 100.0) if h_tot else None,
    )


def compute_contact_map(
    parsed: ParsedStructure,
    chain_id: str,
    atom_mode: Literal["CA", "CB"] = "CA",
    threshold_angstrom: float = 8.0,
    min_seq_sep: int = 1,
    max_residues: int = 1000,
) -> ContactMapResult:
    """Compute residue-residue distance and binary contact maps for a selected chain.

    In CB mode, Glycine (which lacks a C-beta atom) and any residue missing CB fall back
    explicitly to CA, and fallback counts are reported transparently.
    """
    residues = parsed.polymer_residues.get(chain_id, [])
    if not residues:
        raise ValueError(f"Chain '{chain_id}' contains no standard polymer residues.")

    truncated_from: int | None = None
    if len(residues) > max_residues:
        truncated_from = len(residues)
        residues = residues[:max_residues]

    coords: list[np.ndarray] = []
    labels: list[str] = []
    numbers: list[int] = []
    gly_fallbacks = 0
    other_fallbacks = 0

    for res in residues:
        coord: np.ndarray | None = None
        if atom_mode == "CB":
            if "CB" in res.atoms:
                coord = res.atoms["CB"].coord
            elif "CA" in res.atoms:
                coord = res.atoms["CA"].coord
                if res.resname == "GLY":
                    gly_fallbacks += 1
                else:
                    other_fallbacks += 1
        else:
            if "CA" in res.atoms:
                coord = res.atoms["CA"].coord

        if coord is not None:
            coords.append(coord)
            labels.append(f"{res.aa1}{res.resseq}")
            numbers.append(res.resseq)

    n = len(coords)
    if n == 0:
        raise ValueError(f"No {atom_mode}/CA backbone coordinates found in chain '{chain_id}'.")

    coord_arr = np.vstack(coords)
    diff = coord_arr[:, np.newaxis, :] - coord_arr[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff * diff, axis=-1))

    idx = np.arange(n)
    sep_matrix = np.abs(idx[:, np.newaxis] - idx[np.newaxis, :])
    contact_matrix = (dist_matrix <= float(threshold_angstrom)) & (sep_matrix >= int(min_seq_sep))

    upper_mask = np.triu(contact_matrix, k=1)
    total_contacts = int(np.sum(upper_mask))
    long_range_mask = upper_mask & (sep_matrix >= 12)
    long_range_contacts = int(np.sum(long_range_mask))

    possible_pairs = int(np.sum(np.triu(sep_matrix >= int(min_seq_sep), k=1)))
    density_pct = (total_contacts / possible_pairs * 100.0) if possible_pairs > 0 else 0.0

    pair_i, pair_j = np.where(upper_mask)
    contacts_list = [
        {
            "Residue 1": labels[int(i)],
            "Residue 2": labels[int(j)],
            "Seq Separation": int(numbers[int(j)] - numbers[int(i)]),
            "Distance (Å)": round(float(dist_matrix[i, j]), 2),
        }
        for i, j in zip(pair_i, pair_j)
    ]
    contacts_list.sort(key=lambda item: (-int(item["Seq Separation"]), float(item["Distance (Å)"])))

    return ContactMapResult(
        chain_id=chain_id,
        atom_mode=atom_mode,
        threshold_angstrom=float(threshold_angstrom),
        min_seq_sep=int(min_seq_sep),
        residue_count=n,
        truncated_from=truncated_from,
        glycine_fallback_count=gly_fallbacks,
        missing_cb_fallback_count=other_fallbacks,
        residue_labels=labels,
        residue_numbers=numbers,
        distance_matrix=dist_matrix,
        contact_matrix=contact_matrix,
        total_contacts=total_contacts,
        long_range_contacts=long_range_contacts,
        contact_density_pct=density_pct,
        top_contacts=contacts_list[:50],
    )
