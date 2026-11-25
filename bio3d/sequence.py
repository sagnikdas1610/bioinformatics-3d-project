"""Protein sequence parsing, validation, and physicochemical analysis."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import re

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from bio3d.config import (
    AA_FULL_NAMES,
    AA_GROUPS,
    KYTE_DOOLITTLE,
    STANDARD_AA,
    STANDARD_AA_SET,
)


@dataclass(frozen=True)
class FastaRecord:
    """Single parsed FASTA or raw sequence entry."""

    header: str
    sequence: str


@dataclass(frozen=True)
class SequenceValidation:
    """Validation result for an amino-acid sequence."""

    is_valid: bool
    cleaned_sequence: str
    length: int
    invalid_chars: list[str]
    invalid_positions: list[tuple[int, str]]
    stripped_stop_codon: bool
    message: str


@dataclass(frozen=True)
class SequenceProfile:
    """Physicochemical descriptors derived from primary amino-acid sequence."""

    sequence: str
    length: int
    molecular_weight_da: float
    molecular_weight_kda: float
    theoretical_pi: float
    aromaticity: float
    instability_index: float
    instability_class: str
    gravy: float
    extinction_reduced: int
    extinction_oxidized: int
    net_charge_ph7: float
    helix_propensity: float
    turn_propensity: float
    sheet_propensity: float
    aa_counts: dict[str, int]
    aa_percents: dict[str, float]
    group_counts: dict[str, int]
    group_percents: dict[str, float]


def parse_fasta_records(raw_text: str) -> list[FastaRecord]:
    """Parse raw or multi-FASTA text into distinct FastaRecord entries.

    Unlike naive header stripping, multiple `>` records are preserved separately
    so multi-chain FASTA files are never silently concatenated.
    """
    if not raw_text or not raw_text.strip():
        return []

    records: list[FastaRecord] = []
    current_header: str | None = None
    chunks: list[str] = []

    for line in raw_text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if chunks or current_header is not None:
                seq = "".join(chunks).upper()
                if seq:
                    records.append(FastaRecord(header=current_header or "Sequence", sequence=seq))
            current_header = stripped[1:].strip() or f"Sequence {len(records) + 1}"
            chunks = []
        else:
            chunks.append(re.sub(r"\s+", "", stripped))

    if chunks:
        seq = "".join(chunks).upper()
        if seq:
            records.append(FastaRecord(header=current_header or "Pasted sequence", sequence=seq))

    return records


def validate_sequence(raw_sequence: str) -> SequenceValidation:
    """Validate that a cleaned sequence uses only the 20 standard amino acids."""
    cleaned = re.sub(r"\s+", "", raw_sequence or "").upper()
    stripped_stop = False
    if cleaned.endswith("*"):
        cleaned = cleaned[:-1]
        stripped_stop = True

    if not cleaned:
        return SequenceValidation(
            is_valid=False,
            cleaned_sequence="",
            length=0,
            invalid_chars=[],
            invalid_positions=[],
            stripped_stop_codon=stripped_stop,
            message="Provide a non-empty amino-acid sequence using the 20 standard 1-letter codes.",
        )

    invalid_positions = [
        (idx + 1, ch) for idx, ch in enumerate(cleaned) if ch not in STANDARD_AA_SET
    ]
    invalid_chars = sorted({ch for _, ch in invalid_positions})

    if invalid_chars:
        sample_locs = ", ".join(f"{ch}@{pos}" for pos, ch in invalid_positions[:8])
        suffix = " ..." if len(invalid_positions) > 8 else ""
        return SequenceValidation(
            is_valid=False,
            cleaned_sequence=cleaned,
            length=len(cleaned),
            invalid_chars=invalid_chars,
            invalid_positions=invalid_positions,
            stripped_stop_codon=stripped_stop,
            message=(
                f"Unsupported residue symbol(s): {', '.join(invalid_chars)} "
                f"(e.g., {sample_locs}{suffix}). Supported alphabet: {STANDARD_AA}."
            ),
        )

    note = " Trailing stop symbol (*) was removed." if stripped_stop else ""
    return SequenceValidation(
        is_valid=True,
        cleaned_sequence=cleaned,
        length=len(cleaned),
        invalid_chars=[],
        invalid_positions=[],
        stripped_stop_codon=stripped_stop,
        message=f"Valid sequence ({len(cleaned)} standard residues).{note}",
    )


def analyze_sequence(sequence: str) -> SequenceProfile:
    """Compute sequence-derived physicochemical properties via Biopython ProtParam."""
    validation = validate_sequence(sequence)
    if not validation.is_valid:
        raise ValueError(validation.message)

    seq = validation.cleaned_sequence
    analysis = ProteinAnalysis(seq)
    mw_da = float(analysis.molecular_weight())
    instability = float(analysis.instability_index())
    ext_reduced, ext_oxidized = analysis.molar_extinction_coefficient()
    helix_p, turn_p, sheet_p = analysis.secondary_structure_fraction()

    counts = Counter(seq)
    length = len(seq)
    aa_counts = {aa: counts.get(aa, 0) for aa in STANDARD_AA}
    aa_percents = {aa: (aa_counts[aa] / length) * 100.0 for aa in STANDARD_AA}

    group_counts: dict[str, int] = {}
    for aa, count in aa_counts.items():
        grp = AA_GROUPS[aa]
        group_counts[grp] = group_counts.get(grp, 0) + count
    group_percents = {grp: (cnt / length) * 100.0 for grp, cnt in group_counts.items()}

    return SequenceProfile(
        sequence=seq,
        length=length,
        molecular_weight_da=mw_da,
        molecular_weight_kda=mw_da / 1000.0,
        theoretical_pi=float(analysis.isoelectric_point()),
        aromaticity=float(analysis.aromaticity()),
        instability_index=instability,
        instability_class="Stable heuristic (< 40)" if instability < 40.0 else "Unstable heuristic (>= 40)",
        gravy=float(analysis.gravy()),
        extinction_reduced=int(ext_reduced),
        extinction_oxidized=int(ext_oxidized),
        net_charge_ph7=float(analysis.charge_at_pH(7.0)),
        helix_propensity=float(helix_p),
        turn_propensity=float(turn_p),
        sheet_propensity=float(sheet_p),
        aa_counts=aa_counts,
        aa_percents=aa_percents,
        group_counts=group_counts,
        group_percents=group_percents,
    )


def kyte_doolittle_profile(sequence: str, window: int = 9) -> tuple[list[int], list[float]]:
    """Compute sliding-window Kyte-Doolittle hydropathy scores along a sequence."""
    validation = validate_sequence(sequence)
    if not validation.is_valid:
        raise ValueError(validation.message)

    seq = validation.cleaned_sequence
    eff_window = max(1, min(int(window), len(seq)))
    values = [KYTE_DOOLITTLE[aa] for aa in seq]

    positions: list[int] = []
    scores: list[float] = []
    half = eff_window // 2
    for start in range(len(values) - eff_window + 1):
        window_slice = values[start : start + eff_window]
        positions.append(start + half + 1)
        scores.append(sum(window_slice) / eff_window)
    return positions, scores


def charge_vs_ph_curve(sequence: str) -> tuple[list[float], list[float]]:
    """Compute theoretical net protein charge across pH 1.0 to 14.0."""
    validation = validate_sequence(sequence)
    if not validation.is_valid:
        raise ValueError(validation.message)

    analysis = ProteinAnalysis(validation.cleaned_sequence)
    phs = [1.0 + 0.5 * step for step in range(27)]  # 1.0 .. 14.0
    charges = [float(analysis.charge_at_pH(ph)) for ph in phs]
    return phs, charges


def composition_rows(profile: SequenceProfile) -> list[dict[str, object]]:
    """Return structured amino-acid composition rows for display and export."""
    return [
        {
            "AA": aa,
            "Name": AA_FULL_NAMES[aa],
            "Class": AA_GROUPS[aa],
            "Count": profile.aa_counts[aa],
            "Percentage (%)": round(profile.aa_percents[aa], 2),
        }
        for aa in STANDARD_AA
    ]
