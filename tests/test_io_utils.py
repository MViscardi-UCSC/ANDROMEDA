"""
test_io_utils.py
Marcus Viscardi,    March 04, 2025

Tests for the io_utils functions in the Andromeda package.
"""
import pytest
from pathlib import Path
from andromeda.io_utils import load_single_reference_contig, load_reference_contigs_to_dict, load_umi_positions


def mock_fasta_file():
    return ">expected_contig_id\nexpected_sequence\n>second_contig_id\nsecond_sequence"


def mock_umi_positions_file():
    return ("ref_contig,start,end,date_time\n"
            "expected_contig_id,1,7,2025-03-28_14:50:15\n"
            "second_contig_id,1,7,2025-03-28_14:50:15\n")


def test_load_single_reference_contig(tmp_path):
    quick_fasta = tmp_path / "quick.fasta"
    quick_fasta.write_text(mock_fasta_file())
    contig, seq = load_single_reference_contig(quick_fasta)
    assert contig == "expected_contig_id"
    assert seq == "expected_sequence"


def test_load_reference_contigs_to_dict(tmp_path):
    quick_fasta = tmp_path / "quick.fasta"
    quick_fasta.write_text(mock_fasta_file())
    ref_dict = load_reference_contigs_to_dict(quick_fasta)
    assert ref_dict["expected_contig_id"] == "expected_sequence"
    assert ref_dict["second_contig_id"] == "second_sequence"


def test_load_umi_positions(tmp_path):
    quick_umi_positions = tmp_path / "quick_umi_positions.csv"
    quick_umi_positions.write_text(mock_umi_positions_file())
    umi_dict = load_umi_positions(quick_umi_positions)
    assert umi_dict["expected_contig_id"] == (0, (1, 7))
    assert umi_dict["second_contig_id"] == (0, (1, 7))
