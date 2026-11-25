"""Tests for sequence parsing, validation, and physicochemical analysis."""

from __future__ import annotations

import pytest

from bio3d.config import UBIQUITIN_SEQUENCE
from bio3d.sequence import (
    analyze_sequence,
    charge_vs_ph_curve,
    kyte_doolittle_profile,
    parse_fasta_records,
    validate_sequence,
)


def test_parse_multi_fasta_does_not_concatenate():
    raw = ">chain_A\nMVLS\nPADK\n>chain_B\nMVHL\nTPEE\n"
    records = parse_fasta_records(raw)
    assert len(records) == 2
    assert records[0].header == "chain_A"
    assert records[0].sequence == "MVLSPADK"
    assert records[1].header == "chain_B"
    assert records[1].sequence == "MVHLTPEE"


def test_validate_sequence_flags_invalid_positions():
    res = validate_sequence("ACDXEF*")
    assert not res.is_valid
    assert res.stripped_stop_codon is True
    assert res.invalid_chars == ["X"]
    assert res.invalid_positions == [(4, "X")]


def test_ubiquitin_sequence_metrics_match_reference():
    profile = analyze_sequence(UBIQUITIN_SEQUENCE)
    assert profile.length == 76
    # Human ubiquitin (1UBQ, 76 aa) theoretical MW ~ 8564.8 Da, pI ~ 6.56
    assert 8555.0 < profile.molecular_weight_da < 8575.0
    assert 6.3 < profile.theoretical_pi < 6.8
    assert sum(profile.aa_counts.values()) == 76
    assert abs(sum(profile.aa_percents.values()) - 100.0) < 1e-6


def test_kyte_doolittle_and_charge_curve():
    positions, scores = kyte_doolittle_profile(UBIQUITIN_SEQUENCE, window=9)
    assert len(positions) == len(UBIQUITIN_SEQUENCE) - 9 + 1
    assert positions[0] == 5
    phs, charges = charge_vs_ph_curve(UBIQUITIN_SEQUENCE)
    assert len(phs) == len(charges) == 27
    assert charges[0] > 0  # acidic pH -> positive net charge
    assert charges[-1] < 0  # basic pH -> negative net charge


def test_short_sequence_hydropathy_window_clamps():
    positions, scores = kyte_doolittle_profile("ACDE", window=15)
    assert len(positions) == 1
    assert len(scores) == 1


def test_invalid_sequence_raises():
    with pytest.raises(ValueError):
        analyze_sequence("ACDE123")
