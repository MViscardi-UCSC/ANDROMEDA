"""
test_bam_tools.py
Marcus Viscardi,    March 04, 2025


"""
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import pysam
from andromeda.bam_tools import extract_reference_path, histogram_tag_values, scatter_tag_v_tag, barplot_contig_counts

@pytest.fixture
def mock_files_dir():
    return Path(__file__).parent / "mock_files"


@pytest.fixture
def mock_bam_file(mock_files_dir):
    return mock_files_dir / "mock_mapped_reads.bam"


def test_extract_reference_path(mock_bam_file):
    ref_path = extract_reference_path(mock_bam_file)
    assert ref_path == "nanoluc_untranslatable"
