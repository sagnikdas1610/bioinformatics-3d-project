"""Tests for pairwise sequence alignment (global, local, identity vs similarity)."""

from __future__ import annotations

import pytest

from bio3d.alignment import format_alignment_blocks, parse_sequence_identity, run_pairwise_alignment
from bio3d.config import HBA_HUMAN_SEQ, HBB_HUMAN_SEQ


def test_identical_sequences_global():
    res = run_pairwise_alignment("MKTAYIAKQRQ", "MKTAYIAKQRQ", mode="global", matrix_name="BLOSUM62")
    assert res.identity_pct == pytest.approx(100.0)
    assert res.similarity_pct == pytest.approx(100.0)
    assert res.gap_count == 0
    assert res.seq1_span == (1, 11)


def test_identity_vs_similarity_distinction():
    # M=M identical; K<->R, D<->E, R<->K are non-identical conservative substitutions (>0 in BLOSUM62)
    res = run_pairwise_alignment("MKDR", "MREK", mode="global", matrix_name="BLOSUM62")
    assert res.identical_count == 1  # M
    assert res.similar_count == 4  # M, K-R (+2), D-E (+2), R-K (+2)
    assert res.identity_pct == pytest.approx(25.0)
    assert res.similarity_pct == pytest.approx(100.0)
    assert res.similarity_pct > res.identity_pct

    # Include a non-conservative substitution (A<->C = 0 in BLOSUM62)
    res2 = run_pairwise_alignment("MKDA", "MREC", mode="global", matrix_name="BLOSUM62")
    assert res2.identical_count == 1
    assert res2.similar_count == 3  # M, K-R, D-E
    assert res2.identity_pct == pytest.approx(25.0)
    assert res2.similarity_pct == pytest.approx(75.0)


def test_local_alignment_extracts_matching_domain():
    s1 = "GGGGGGGGMKTAYIAKQRQGGGGGGGG"
    s2 = "PPPPPMKTAYIAKQRQTTTTT"
    res = run_pairwise_alignment(s1, s2, mode="local", matrix_name="BLOSUM62")
    assert res.aligned_seq1 == "MKTAYIAKQRQ"
    assert res.aligned_seq2 == "MKTAYIAKQRQ"
    assert res.identity_pct == pytest.approx(100.0)
    assert res.seq1_span == (9, 19)
    assert res.seq2_span == (6, 16)


def test_hemoglobin_alpha_beta_alignment():
    res = run_pairwise_alignment(HBA_HUMAN_SEQ, HBB_HUMAN_SEQ, mode="global", matrix_name="BLOSUM62")
    # Human HBA vs HBB shares ~40-46% identity and >55% similarity
    assert 35.0 <= res.identity_pct <= 52.0
    assert res.similarity_pct > res.identity_pct
    blocks = format_alignment_blocks(res, width=60)
    assert "Seq1" in blocks and "Seq2" in blocks


def test_invalid_gap_penalty_raises():
    with pytest.raises(ValueError):
        run_pairwise_alignment("MKTA", "MKTA", open_gap_score=5.0)


def test_parse_sequence_identity_extracts_accessions():
    n1, a1 = parse_sequence_identity("sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens")
    assert a1 == "P69905"
    assert "HBA_HUMAN" in n1

    n2, a2 = parse_sequence_identity("1UBQ_chain_A")
    assert a2 == "1UBQ"
    assert n2 == "1UBQ_chain_A"

    n3, a3 = parse_sequence_identity("Pasted sequence", fallback_name="Sequence 1")
    assert n3 == "Sequence 1"
    assert a3 is None
