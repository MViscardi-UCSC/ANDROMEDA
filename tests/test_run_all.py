"""
test_run_all.py
Marcus Viscardi,    March 04, 2025

Script to test out the run-all functionality of our toy data for the Andromeda package.
"""
import sys

import pytest
from pathlib import Path
from andromeda.cli import parse_args, run_all_pipeline
import pysam


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

def test_run_all_pipeline(mock_ref_fasta, mock_bam_file, mock_umi_positions, tmp_path, monkeypatch):
    # Set up arguments for the pipeline
    test_args = [
        "run-all",
        str(mock_ref_fasta),
        str(mock_bam_file),
        str(tmp_path / "output"),
        "--umi-positions",
        str(mock_umi_positions),
    ]
    monkeypatch.setattr(sys, "argv", ["andromeda"] + test_args)

    parsed_args = parse_args()
    run_all_pipeline(parsed_args)

    # Verify the output
    output_dir = tmp_path / "output"
    assert output_dir.exists()
    
    # Verify the outputs of the extraction/tagging step
    assert (output_dir / "tagging").exists()
    assert (output_dir / "tagging" / "mock_mapped_reads.tagged.bam").exists()
    assert (output_dir / "tagging" / "mock_mapped_reads.tagged.sam").exists()
    # We should check the length of the sam file to see if it's 38 lines long:
    with open(output_dir / "tagging" / "mock_mapped_reads.tagged.sam") as f:
        assert len(f.readlines()) == 36
    
    # Verify the outputs of the grouping step
    assert (output_dir / "grouped").exists()
    assert (output_dir / "grouped" / "mock_mapped_reads.tagged.grouped_0dist.bam").exists()
    assert (output_dir / "grouped" / "mock_mapped_reads.tagged.grouped_0dist.tsv").exists()
    
    # Verify the outputs of the consensus step
    assert (output_dir / "consensus").exists()
    assert (output_dir / "consensus" / "mock_mapped_reads.tagged.grouped_0dist.consensus_sequences.bam").exists()
    assert (output_dir / "consensus" / "mock_mapped_reads.tagged.grouped_0dist.consensus_sequences.sorted.bam").exists()
    assert (output_dir / "consensus" / "mock_mapped_reads.tagged.grouped_0dist_consensus_sequences.tsv").exists()
    # We should end up with 9 consensus sequences
    final_bam = pysam.AlignmentFile(output_dir /
                                    "consensus" /
                                    "mock_mapped_reads.tagged.grouped_0dist.consensus_sequences.sorted.bam")
    assert final_bam.count() == 9, "Incorrect number of consensus sequences"
    
