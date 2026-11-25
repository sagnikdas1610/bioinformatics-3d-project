"""Pairwise protein sequence alignment using Biopython PairwiseAligner or Gotoh DP."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import re
from typing import Literal

import Bio
import numpy as np

from bio3d.sequence import validate_sequence

SUPPORTED_MATRICES = ("BLOSUM62", "BLOSUM80", "BLOSUM45", "PAM250", "IDENTITY")
_MATRIX_DATA_DIR = Path(Bio.__file__).resolve().parent / "Align" / "substitution_matrices" / "data"


@dataclass(frozen=True)
class AlignmentResult:
    """Complete global or local pairwise alignment summary."""

    mode: Literal["global", "local"]
    matrix_name: str
    open_gap_score: float
    extend_gap_score: float
    score: float
    aligned_seq1: str
    aligned_seq2: str
    match_line: str
    alignment_length: int
    aligned_pairs: int
    identical_count: int
    similar_count: int
    gap_count: int
    identity_pct: float
    identity_ungapped_pct: float
    similarity_pct: float
    gap_pct: float
    seq1_span: tuple[int, int]
    seq2_span: tuple[int, int]
    seq1_length: int
    seq2_length: int
    engine: str
    interpretation: str
    seq1_name: str = "Sequence 1"
    seq1_accession: str | None = None
    seq2_name: str = "Sequence 2"
    seq2_accession: str | None = None


def parse_sequence_identity(header: str, fallback_name: str = "Sequence") -> tuple[str, str | None]:
    """Extract a clean sequence display name and genuine accession ID from a FASTA header."""
    raw = (header or "").strip().lstrip(">").strip()
    if not raw or raw in {"Pasted sequence", "Sequence"}:
        return fallback_name, None

    # UniProt pipe format: sp|P69905|HBA_HUMAN Hemoglobin subunit alpha ...
    m_uniprot = re.match(r"^(?:sp|tr)\|([A-Z0-9]{6,10})\|([^\s]+)(?:\s+(.*))?$", raw, re.IGNORECASE)
    if m_uniprot:
        acc = m_uniprot.group(1).upper()
        entry = m_uniprot.group(2)
        desc = (m_uniprot.group(3) or "").strip()
        short_desc = desc.split(" OS=")[0].strip() if desc else ""
        name = f"{entry} ({short_desc})" if short_desc else entry
        return name, acc

    # Prefixed PDB/UniProt/RefSeq tokens in header (e.g., 1UBQ_chain_A, AF-P69905-F1, NP_000509.1)
    first_tok = raw.split()[0]
    m_af = re.search(r"\bAF-([A-Z0-9]{6,10})-F\d+\b", raw, re.IGNORECASE)
    if m_af:
        return raw, m_af.group(1).upper()
    m_pdb = re.match(r"^([A-Z0-9]{4})(?:[_:|].*)?$", first_tok, re.IGNORECASE)
    if m_pdb and any(ch.isdigit() for ch in m_pdb.group(1)):
        return raw, m_pdb.group(1).upper()
    m_acc = re.match(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}|[NX]P_[0-9]+(?:\.[0-9]+)?)$", first_tok, re.IGNORECASE)
    if m_acc:
        return raw, m_acc.group(1).upper()

    return raw, None


def _load_plain_matrix(matrix_name: str) -> dict[tuple[str, str], float]:
    """Load an NCBI/Biopython substitution matrix directly from its text file."""
    path = _MATRIX_DATA_DIR / matrix_name
    if not path.exists():
        raise FileNotFoundError(f"Substitution matrix file not found: {path}")

    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    cols = lines[0].split()
    table: dict[tuple[str, str], float] = {}
    for row_line in lines[1:]:
        parts = row_line.split()
        row_aa = parts[0]
        vals = [float(v) for v in parts[1:]]
        for col_aa, val in zip(cols, vals):
            table[(row_aa, col_aa)] = val
    return table


def _reconstruct_gapped_strings(
    seq1: str,
    seq2: str,
    blocks1: list[tuple[int, int]],
    blocks2: list[tuple[int, int]],
    mode: Literal["global", "local"],
) -> tuple[str, str, tuple[int, int], tuple[int, int]]:
    """Reconstruct gapped alignment strings cleanly from aligned coordinate blocks."""
    if not blocks1 or not blocks2:
        return "", "", (0, 0), (0, 0)

    out1: list[str] = []
    out2: list[str] = []

    if mode == "global":
        pos1, pos2 = 0, 0
        span1_start, span2_start = 1, 1
    else:
        pos1 = int(blocks1[0][0])
        pos2 = int(blocks2[0][0])
        span1_start = pos1 + 1
        span2_start = pos2 + 1

    for (s1, e1), (s2, e2) in zip(blocks1, blocks2):
        s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
        if s1 > pos1:
            seg = seq1[pos1:s1]
            out1.append(seg)
            out2.append("-" * len(seg))
        if s2 > pos2:
            seg = seq2[pos2:s2]
            out1.append("-" * len(seg))
            out2.append(seg)
        out1.append(seq1[s1:e1])
        out2.append(seq2[s2:e2])
        pos1, pos2 = e1, e2

    if mode == "global":
        if pos1 < len(seq1):
            seg = seq1[pos1:]
            out1.append(seg)
            out2.append("-" * len(seg))
        if pos2 < len(seq2):
            seg = seq2[pos2:]
            out1.append("-" * len(seg))
            out2.append(seg)
        span1_end, span2_end = len(seq1), len(seq2)
    else:
        span1_end = int(blocks1[-1][1])
        span2_end = int(blocks2[-1][1])

    return "".join(out1), "".join(out2), (span1_start, span1_end), (span2_start, span2_end)


def _gotoh_align(
    seq1: str,
    seq2: str,
    mode: Literal["global", "local"],
    score_fn,
    open_gap: float,
    extend_gap: float,
) -> tuple[float, str, str, tuple[int, int], tuple[int, int]]:
    """Exact Gotoh (1982) affine-gap global/local alignment in pure Python/NumPy."""
    n, m = len(seq1), len(seq2)
    neg_inf = -1e15

    # State encodings: 0=M, 1=X (gap in seq2), 2=Y (gap in seq1)
    M = np.full((n + 1, m + 1), neg_inf, dtype=np.float64)
    X = np.full((n + 1, m + 1), neg_inf, dtype=np.float64)
    Y = np.full((n + 1, m + 1), neg_inf, dtype=np.float64)

    bt_m = np.full((n + 1, m + 1), -1, dtype=np.int8)
    bt_x = np.full((n + 1, m + 1), -1, dtype=np.int8)
    bt_y = np.full((n + 1, m + 1), -1, dtype=np.int8)

    if mode == "global":
        M[0, 0] = 0.0
        for i in range(1, n + 1):
            X[i, 0] = open_gap + (i - 1) * extend_gap
            bt_x[i, 0] = 0 if i == 1 else 1
        for j in range(1, m + 1):
            Y[0, j] = open_gap + (j - 1) * extend_gap
            bt_y[0, j] = 0 if j == 1 else 2
    else:
        M[0, :] = 0.0
        M[:, 0] = 0.0
        X[0, :] = 0.0
        X[:, 0] = 0.0
        Y[0, :] = 0.0
        Y[:, 0] = 0.0

    best_local_score = 0.0
    best_local_pos = (0, 0)

    for i in range(1, n + 1):
        c1 = seq1[i - 1]
        for j in range(1, m + 1):
            sub = score_fn(c1, seq2[j - 1])

            # X[i, j]: seq1[i-1] aligned to gap
            x_from_m = M[i - 1, j] + open_gap
            x_from_x = X[i - 1, j] + extend_gap
            x_from_y = Y[i - 1, j] + open_gap
            if x_from_m >= x_from_x and x_from_m >= x_from_y:
                X[i, j], bt_x[i, j] = x_from_m, 0
            elif x_from_x >= x_from_y:
                X[i, j], bt_x[i, j] = x_from_x, 1
            else:
                X[i, j], bt_x[i, j] = x_from_y, 2
            if mode == "local" and X[i, j] <= 0.0:
                X[i, j], bt_x[i, j] = 0.0, -1

            # Y[i, j]: gap aligned to seq2[j-1]
            y_from_m = M[i, j - 1] + open_gap
            y_from_y = Y[i, j - 1] + extend_gap
            y_from_x = X[i, j - 1] + open_gap
            if y_from_m >= y_from_y and y_from_m >= y_from_x:
                Y[i, j], bt_y[i, j] = y_from_m, 0
            elif y_from_y >= y_from_x:
                Y[i, j], bt_y[i, j] = y_from_y, 2
            else:
                Y[i, j], bt_y[i, j] = y_from_x, 1
            if mode == "local" and Y[i, j] <= 0.0:
                Y[i, j], bt_y[i, j] = 0.0, -1

            # M[i, j]: match/mismatch
            m_prev, x_prev, y_prev = M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1]
            if m_prev >= x_prev and m_prev >= y_prev:
                val, prev_state = m_prev + sub, 0
            elif x_prev >= y_prev:
                val, prev_state = x_prev + sub, 1
            else:
                val, prev_state = y_prev + sub, 2

            if mode == "local" and val <= 0.0:
                M[i, j], bt_m[i, j] = 0.0, -1
            else:
                M[i, j], bt_m[i, j] = val, prev_state
                if mode == "local" and val > best_local_score:
                    best_local_score = float(val)
                    best_local_pos = (i, j)

    aln1: list[str] = []
    aln2: list[str] = []

    if mode == "global":
        end_scores = (M[n, m], X[n, m], Y[n, m])
        state = int(np.argmax(end_scores))
        final_score = float(end_scores[state])
        i, j = n, m
        while i > 0 or j > 0:
            if state == 0:
                prev = int(bt_m[i, j])
                aln1.append(seq1[i - 1])
                aln2.append(seq2[j - 1])
                i -= 1
                j -= 1
                state = prev
            elif state == 1:
                prev = int(bt_x[i, j])
                aln1.append(seq1[i - 1])
                aln2.append("-")
                i -= 1
                state = prev
            else:
                prev = int(bt_y[i, j])
                aln1.append("-")
                aln2.append(seq2[j - 1])
                j -= 1
                state = prev
        span1, span2 = (1, n), (1, m)
    else:
        if best_local_score <= 0.0:
            return 0.0, "", "", (0, 0), (0, 0)
        final_score = best_local_score
        end_i, end_j = best_local_pos
        i, j = end_i, end_j
        state = 0
        while i > 0 and j > 0 and state != -1:
            if state == 0:
                if M[i, j] <= 0.0:
                    break
                prev = int(bt_m[i, j])
                aln1.append(seq1[i - 1])
                aln2.append(seq2[j - 1])
                i -= 1
                j -= 1
                state = prev
            elif state == 1:
                if X[i, j] <= 0.0:
                    break
                prev = int(bt_x[i, j])
                aln1.append(seq1[i - 1])
                aln2.append("-")
                i -= 1
                state = prev
            else:
                if Y[i, j] <= 0.0:
                    break
                prev = int(bt_y[i, j])
                aln1.append("-")
                aln2.append(seq2[j - 1])
                j -= 1
                state = prev
        span1 = (i + 1, end_i)
        span2 = (j + 1, end_j)

    aln1.reverse()
    aln2.reverse()
    return final_score, "".join(aln1), "".join(aln2), span1, span2


def _interpret_alignment(
    mode: str,
    identity_pct: float,
    similarity_pct: float,
    aligned_pairs: int,
    matrix_name: str,
) -> str:
    """Generate a careful scientific interpretation avoiding homology overclaims."""
    if aligned_pairs == 0:
        return "No aligned residue pairs were found under the selected scoring parameters."

    if identity_pct >= 50.0:
        tier = "High sequence identity (>= 50%), typically consistent with closely related homologs sharing fold and often conserved function."
    elif identity_pct >= 30.0:
        tier = "Moderate sequence identity (30-50%), generally indicative of shared structural fold, though functional divergence is possible."
    elif identity_pct >= 20.0:
        tier = "Low ('twilight zone', 20-30%) sequence identity; structural homology is possible but requires structural or profile-based confirmation."
    else:
        tier = "Very low sequence identity (< 20%); sequence similarity alone does not establish homology."

    sim_note = (
        f"Exact identity ({identity_pct:.1f}%) counts identical residues at aligned columns, whereas "
        f"similarity ({similarity_pct:.1f}%) also includes conservative substitutions scoring > 0 in {matrix_name}."
        if matrix_name != "IDENTITY"
        else f"With the IDENTITY matrix, similarity ({similarity_pct:.1f}%) equals identity ({identity_pct:.1f}%) because only exact matches score positively."
    )
    return f"[{mode.capitalize()} alignment] {tier} {sim_note}"


def run_pairwise_alignment(
    seq1: str,
    seq2: str,
    mode: Literal["global", "local"] = "global",
    matrix_name: str = "BLOSUM62",
    open_gap_score: float = -10.0,
    extend_gap_score: float = -0.5,
    seq1_name: str = "Sequence 1",
    seq1_accession: str | None = None,
    seq2_name: str = "Sequence 2",
    seq2_accession: str | None = None,
) -> AlignmentResult:
    """Run Needleman-Wunsch (global) or Smith-Waterman (local) pairwise alignment."""
    v1 = validate_sequence(seq1)
    if not v1.is_valid:
        raise ValueError(f"Sequence 1 is invalid: {v1.message}")
    v2 = validate_sequence(seq2)
    if not v2.is_valid:
        raise ValueError(f"Sequence 2 is invalid: {v2.message}")

    clean1 = v1.cleaned_sequence
    clean2 = v2.cleaned_sequence

    if matrix_name not in SUPPORTED_MATRICES:
        raise ValueError(f"Unsupported matrix '{matrix_name}'. Choose one of {SUPPORTED_MATRICES}.")
    if open_gap_score > 0 or extend_gap_score > 0:
        raise ValueError("Gap open and gap extend scores must be <= 0 (penalties).")

    sub_table = None if matrix_name == "IDENTITY" else _load_plain_matrix(matrix_name)

    try:
        from Bio.Align import PairwiseAligner, substitution_matrices

        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.open_gap_score = float(open_gap_score)
        aligner.extend_gap_score = float(extend_gap_score)
        if matrix_name == "IDENTITY":
            aligner.match_score = 1.0
            aligner.mismatch_score = -1.0
        else:
            aligner.substitution_matrix = substitution_matrices.load(matrix_name)
        best = aligner.align(clean1, clean2)[0]
        blocks1 = [(int(s), int(e)) for s, e in best.aligned[0]]
        blocks2 = [(int(s), int(e)) for s, e in best.aligned[1]]
        aligned1, aligned2, span1, span2 = _reconstruct_gapped_strings(
            clean1, clean2, blocks1, blocks2, mode
        )
        score = float(best.score)
        engine = "Biopython PairwiseAligner"
    except ImportError:
        if sub_table is None:
            def score_fn(a: str, b: str) -> float:
                return 1.0 if a == b else -1.0
        else:
            def score_fn(a: str, b: str) -> float:
                return sub_table[(a, b)]

        score, aligned1, aligned2, span1, span2 = _gotoh_align(
            clean1,
            clean2,
            mode=mode,
            score_fn=score_fn,
            open_gap=float(open_gap_score),
            extend_gap=float(extend_gap_score),
        )
        engine = "Gotoh Affine-Gap DP (NCBI matrix)"

    match_chars: list[str] = []
    identical = 0
    similar = 0
    gaps = 0
    aligned_pairs = 0

    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            gaps += 1
            match_chars.append(" ")
            continue
        aligned_pairs += 1
        if a == b:
            identical += 1
            similar += 1
            match_chars.append("|")
        else:
            is_similar = sub_table is not None and sub_table[(a, b)] > 0.0
            if is_similar:
                similar += 1
                match_chars.append(":")
            else:
                match_chars.append(".")

    aln_len = len(aligned1)
    identity_pct = (identical / aln_len * 100.0) if aln_len else 0.0
    identity_ungapped_pct = (identical / aligned_pairs * 100.0) if aligned_pairs else 0.0
    similarity_pct = (similar / aln_len * 100.0) if aln_len else 0.0
    gap_pct = (gaps / aln_len * 100.0) if aln_len else 0.0

    return AlignmentResult(
        mode=mode,
        matrix_name=matrix_name,
        open_gap_score=float(open_gap_score),
        extend_gap_score=float(extend_gap_score),
        score=score,
        aligned_seq1=aligned1,
        aligned_seq2=aligned2,
        match_line="".join(match_chars),
        alignment_length=aln_len,
        aligned_pairs=aligned_pairs,
        identical_count=identical,
        similar_count=similar,
        gap_count=gaps,
        identity_pct=identity_pct,
        identity_ungapped_pct=identity_ungapped_pct,
        similarity_pct=similarity_pct,
        gap_pct=gap_pct,
        seq1_span=span1,
        seq2_span=span2,
        seq1_length=len(clean1),
        seq2_length=len(clean2),
        engine=engine,
        interpretation=_interpret_alignment(
            mode, identity_pct, similarity_pct, aligned_pairs, matrix_name
        ),
        seq1_name=seq1_name,
        seq1_accession=seq1_accession,
        seq2_name=seq2_name,
        seq2_accession=seq2_accession,
    )


def format_alignment_blocks(result: AlignmentResult, width: int = 60) -> str:
    """Format a pairwise alignment into numbered wrapped text blocks."""
    if not result.aligned_seq1:
        return "(empty alignment)"

    lines: list[str] = []
    pos1 = result.seq1_span[0] - 1
    pos2 = result.seq2_span[0] - 1

    for start in range(0, result.alignment_length, width):
        end = min(start + width, result.alignment_length)
        s1 = result.aligned_seq1[start:end]
        ml = result.match_line[start:end]
        s2 = result.aligned_seq2[start:end]

        res1_in_block = len(s1) - s1.count("-")
        res2_in_block = len(s2) - s2.count("-")
        s1_start = pos1 + 1 if res1_in_block else pos1
        s2_start = pos2 + 1 if res2_in_block else pos2
        pos1 += res1_in_block
        pos2 += res2_in_block

        lines.append(f"Seq1 {s1_start:>5d}  {s1}  {pos1}")
        lines.append(f"            {ml}")
        lines.append(f"Seq2 {s2_start:>5d}  {s2}  {pos2}")
        lines.append("")

    return "\n".join(lines).rstrip()


def format_alignment_html(result: AlignmentResult, width: int = 60) -> str:
    """Render a color-coded monospace HTML block for interactive inspection."""
    if not result.aligned_seq1:
        return "<p>No aligned region found.</p>"

    blocks: list[str] = []
    pos1 = result.seq1_span[0] - 1
    pos2 = result.seq2_span[0] - 1

    for start in range(0, result.alignment_length, width):
        end = min(start + width, result.alignment_length)
        s1 = result.aligned_seq1[start:end]
        ml = result.match_line[start:end]
        s2 = result.aligned_seq2[start:end]

        res1_in_block = len(s1) - s1.count("-")
        res2_in_block = len(s2) - s2.count("-")
        s1_start = pos1 + 1 if res1_in_block else pos1
        s2_start = pos2 + 1 if res2_in_block else pos2
        pos1 += res1_in_block
        pos2 += res2_in_block

        styled_match: list[str] = []
        for ch in ml:
            if ch == "|":
                styled_match.append('<span style="color:#38bdf8;font-weight:700;">|</span>')
            elif ch == ":":
                styled_match.append('<span style="color:#fbbf24;font-weight:700;">:</span>')
            elif ch == ".":
                styled_match.append('<span style="color:#94a3b8;">.</span>')
            else:
                styled_match.append(" ")

        blocks.append(
            f"<div>"
            f"<span style='color:#94a3b8;'>Seq1 {s1_start:>4d}</span>  "
            f"<span style='color:#e2e8f0;'>{escape(s1)}</span>  "
            f"<span style='color:#64748b;'>{pos1}</span><br>"
            f"<span style='color:#64748b;'>          </span>  "
            f"{''.join(styled_match)}<br>"
            f"<span style='color:#94a3b8;'>Seq2 {s2_start:>4d}</span>  "
            f"<span style='color:#cbd5e1;'>{escape(s2)}</span>  "
            f"<span style='color:#64748b;'>{pos2}</span>"
            f"</div>"
        )

    return (
        "<div style='font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; "
        "font-size: 13px; line-height: 1.5; background: #0f172a; border: 1px solid #1e293b; "
        "border-radius: 10px; padding: 14px 16px; overflow-x: auto; white-space: pre;'>"
        + "<div style='height:10px;'></div>".join(blocks)
        + "</div>"
    )
