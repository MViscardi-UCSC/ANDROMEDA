"""
bam_tools.py
Marcus Viscardi,    February 13, 2025

Place for functions that deal with BAM files specifically.
Mostly just wrappers around `pysam` functions.
"""
from pathlib import Path
from typing import List, Tuple
import pysam


def extract_reference_path(bam_file: Path) -> str:
    """
    Extract the reference path from a BAM file.

    Args:
        bam_file (Path): Path to BAM file.

    Returns:
        str: Path to the reference sequence.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        print(bam.header)
        return bam.get_reference_name(0)


if __name__ == '__main__':
    test_bam = Path("../examples/JA-NP-093/mapped_reads/250128_JANP-093_LT.AtoG.sorted.bam")
    print(extract_reference_path(test_bam))
