"""
test_consensus.py
Marcus Viscardi,    March 04, 2025


"""
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import pysam
import pandas as pd
from andromeda.consensus import (
    extract_umi_groups,
    extract_ref_oriented_sequence,
    extract_ref_oriented_sequences,
    pick_major_base,
    pick_major_base_and_phred,
    compute_majority_consensus,
    create_consensus_cigar,
    call_consensus_from_bam,
    consensus_df_to_bam,
    build_ref_id_to_contig_dict,
)

@pytest.fixture
def mock_bam_file(tmp_path):
    bam_file = tmp_path / "mock.bam"
    with bam_file.open("wb") as bam_file:
        pysam.AlignmentFile(bam_file, "wb", header={"SQ": [{"SN": "mock_contig", "LN": 12}]})
    return bam_file

@pytest.fixture
def mock_ref_fasta(tmp_path):
    ref_fasta = tmp_path / "mock_ref.fasta"
    ref_fasta.write_text(">mock_contig\nACTGACTGACTG")
    return ref_fasta

def test_extract_umi_groups(mock_bam_file):
    with patch("pysam.AlignmentFile") as mock_alignment_file:
        mock_bam = MagicMock()
        mock_bam.mapped = 10
        mock_read = MagicMock()
        mock_read.has_tag.return_value = True
        mock_read.get_tag.return_value = "mock_umi"
        mock_bam.__iter__.return_value = [mock_read] * 10
        mock_alignment_file.return_value.__enter__.return_value = mock_bam

        umi_groups = extract_umi_groups(mock_bam_file, min_group_size=2)
        assert "mock_umi" in umi_groups
        assert len(umi_groups["mock_umi"]) == 10

def test_extract_ref_oriented_sequence():
    ref_seq = "ACTGACTGACTG"
    mock_read = MagicMock(pysam.AlignedSegment())
    mock_read.query_sequence = "ACTG"
    mock_read.query_qualities = [30, 30, 30, 30]
    mock_read.get_aligned_pairs.return_value = [(0, 0), (1, 1), (2, 2), (3, 3)]
    mock_read.reference_start = 0
    mock_read.reference_end = 3

    result = extract_ref_oriented_sequence(mock_read, ref_seq)
    assert result["query_sequence"] == ["A", "C", "T", "G"]
    assert result["ref_sequence"] == ["A", "C", "T", "G"]

def test_extract_ref_oriented_sequences():
    ref_seq = "ACTGACTGACTG"
    mock_read = MagicMock(pysam.AlignedSegment())
    mock_read.query_sequence = "ACTG"
    mock_read.query_qualities = [30, 30, 30, 30]
    mock_read.get_aligned_pairs.return_value = [(0, 0), (1, 1), (2, 2), (3, 3)]
    mock_read.reference_start = 0
    mock_read.reference_end = 3

    result = extract_ref_oriented_sequences([mock_read], ref_seq)
    assert len(result) == 1
    assert result[0]["query_sequence"] == "ACTG         "

def test_pick_major_base():
    column = ["A", "A", "C", "A"]
    major_base, major_fraction = pick_major_base(column)
    assert major_base == "A"
    assert major_fraction == 0.75

def test_pick_major_base_and_phred():
    column = [("A", "I"), ("A", "I"), ("C", "I"), ("A", "I")]
    major_base, avg_phred, major_fraction = pick_major_base_and_phred(column)
    assert major_base == "A"
    assert major_fraction == 0.75
    # Note the change below, this is due to how we decided to calc the average, it's hurt by not being 100% majority!
    assert avg_phred == "H"  # Not I!!

def test_compute_majority_consensus():
    ref_seq = "ACTGACTGACTG"
    mock_read = MagicMock(pysam.AlignedSegment())
    mock_read.query_sequence = "ACTG"
    mock_read.query_qualities = [30, 30, 30, 30]  # This is equal to "????" in phred
    mock_read.get_aligned_pairs.return_value = [(0, 0), (1, 1), (2, 2), (3, 3)]
    mock_read.reference_start = 0
    mock_read.reference_end = 3

    consensus_seq, consensus_phred, match_str, stats = compute_majority_consensus([mock_read], ref_seq)
    assert consensus_seq == "ACTG        "
    assert consensus_phred == "????        "
    assert match_str == "MMMM        "

def test_create_consensus_cigar():
    match_str = "MMMM"
    cigar, start_pos = create_consensus_cigar(match_str)
    assert cigar == "4M"
    assert start_pos == 0

def test_build_ref_id_to_contig_dict(mock_bam_file):
    with patch("pysam.AlignmentFile") as mock_alignment_file:
        mock_bam = MagicMock()
        mock_bam.references = ["mock_contig"]
        mock_alignment_file.return_value.__enter__.return_value = mock_bam

        ref_dict = build_ref_id_to_contig_dict(mock_bam_file)
        assert ref_dict == {0: "mock_contig"}