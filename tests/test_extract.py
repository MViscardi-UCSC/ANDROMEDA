"""
test_extract.py
Marcus Viscardi,    March 04, 2025


"""
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import pysam
import pandas as pd
from andromeda.extract import extract_umis_from_bam, save_umi_counts


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


def test_extract_umis_from_bam(mock_bam_file, mock_ref_fasta, mock_umi_positions, tmp_path):
    output_dir = tmp_path / "output"
    result_path = extract_umis_from_bam(
        bam_file=mock_bam_file,
        reference_file=mock_ref_fasta,
        umi_positions_file=mock_umi_positions,
        output_dir=output_dir,
        flanking_seq_to_capture=0,
        subset_count=-1,
        mismatch_tolerance=0,
        force=False,
        umi_tag="uM"
    )
    assert result_path.exists()


def test_save_umi_counts(mock_bam_file, tmp_path):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    with patch("andromeda.alignment_tools.bam_to_df") as mock_bam_to_df, \
         patch("pandas.Series.to_csv") as mock_to_csv, \
         patch("builtins.open", new_callable=MagicMock) as mock_open:
        mock_bam_to_df.return_value = pd.DataFrame({
            "uM": ["mock_umi"] * 10
        })
        save_umi_counts(mock_bam_file, output_dir)
        mock_bam_to_df.assert_called_once()
        mock_to_csv.assert_called_once()
        mock_open.assert_called_once()
