"""
test_alignment_tools.py
Marcus Viscardi,    March 04, 2025

Tests for the alignment tools in the Andromeda package.
"""
from unittest.mock import MagicMock

import pytest
from pathlib import Path
import pysam
from andromeda.alignment_tools import bam_to_tagged_bam, bam_to_df, extract_ref_and_query_region
from andromeda.io_utils import load_reference_contigs_to_dict, load_umi_positions


@pytest.fixture
def mock_files_dir():
    return Path(__file__).parent / "mock_files"


@pytest.fixture
def mock_bam_file(mock_files_dir):
    return mock_files_dir / "mock_mapped_reads.bam"


@pytest.fixture
def mock_ref_fasta(mock_files_dir):
    return mock_files_dir / "mock_ref.fasta"


@pytest.fixture
def mock_umi_positions(mock_files_dir):
    return mock_files_dir / "mock_umi_positions.csv"


def test_bam_to_tagged_bam_valid_input(mock_bam_file, mock_ref_fasta, mock_umi_positions, tmp_path):
    ref_seq_dict = load_reference_contigs_to_dict(mock_ref_fasta)
    umi_dict = load_umi_positions(mock_umi_positions)
    save_path = tmp_path / "tagged.bam"

    result_path = bam_to_tagged_bam(
        bam_file_path=mock_bam_file,
        ref_seq_dict=ref_seq_dict,
        umi_dict=umi_dict,
        save_path=save_path,
        max_del_in_umi=0,
        max_ins_in_umi=0,
        max_iupac_mismatches=0,
        restrict_to_length=True,
        bam_tag_dict=None,
        subset_count=25,
    )

    assert result_path.exists(), "Tagged BAM file was not created."


def test_bam_to_df_valid_input(mock_bam_file):
    df = bam_to_df(mock_bam_file)
    assert not df.empty, "DataFrame is empty."


def test_extract_ref_and_query_region():
    ref_seq = "ACTGACTGACTG"
    region_start = 0
    region_end = 3
    mock_entry = MagicMock(pysam.AlignedSegment())
    mock_entry.query_sequence = "ACTG"
    mock_entry.query_qualities = [30, 30, 30, 30]
    mock_entry.get_aligned_pairs.return_value = [(0, 0), (1, 1), (2, 2), (3, 3)]
    mock_entry.get_reference_positions.return_value = [0, 1, 2, 3]

    result = extract_ref_and_query_region(mock_entry, ref_seq, region_start, region_end)

    assert result["query_sequence"] == "ACTG", "Query sequence extraction failed."
    assert result["ref_sequence"] == "ACTG", "Reference sequence extraction failed."
