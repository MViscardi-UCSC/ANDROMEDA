"""
bam_tools.py
Marcus Viscardi,    February 13, 2025

Place for functions that deal with BAM files specifically.
Mostly just wrappers around `pysam` functions.
"""
from pathlib import Path
from typing import List, Tuple
import pysam

from tqdm.auto import tqdm

import seaborn as sea
import matplotlib.pyplot as plt


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


def histogram_tag_values(bam_file: Path, tag: str) -> Tuple[List[int], List[int]]:
    """
    Create a histogram of the values in a BAM tag.

    Args:
        bam_file (Path): Path to BAM file.
        tag (str): Tag to create histogram of.

    Returns:
        Tuple[List[int], List[int]]: Tuple of values and counts.
    """
    values = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total_reads = bam.mapped + bam.unmapped
        for read in tqdm(bam, total=total_reads, desc=f"Extracting {tag} tag from BAM"):
            if read.has_tag(tag):
                val = read.get_tag(tag)
                if val == 41:
                    if read.get_tag("ud") <= 1 and read.get_tag("ui") <= 1 and read.get_tag("um") <= 2:
                        print(read.get_tag("uM"),
                              "L:", read.get_tag("ul"),
                              "I:", read.get_tag("ui"),
                              "D:", read.get_tag("ud"),
                              "M:", read.get_tag("um"))
                values.append(val)
    unique_values = set(values)
    value_counts = [values.count(value) for value in unique_values]
    mean = sum(values) / len(values)
    mode = max(set(values), key=values.count)
    std_dev = (sum([(value - mean) ** 2 for value in values]) / len(values)) ** 0.5
    print(f"Mean: {mean}\nStandard Deviation: {std_dev}\nMode: {mode}")
    sea.histplot(values, binwidth=1)
    if std_dev > 10 and mean > 10:
        plt.xlim(mean - (0.5 * std_dev), mean + (0.5 * std_dev))
    plt.title(f"{tag} tag values\n[{bam_file.name}]")
    plt.show()
    return list(unique_values), value_counts

def scatter_tag_v_tag(bam_file: Path, x_tag: str, y_tag: str) -> None:
    """
    Create a scatter plot of two BAM tags.

    Args:
        bam_file (Path): Path to BAM file.
        x_tag (str): X-axis tag.
        y_tag (str): Y-axis tag.
    """
    x_values = []
    y_values = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total_reads = bam.mapped + bam.unmapped
        for read in tqdm(bam, total=total_reads, desc=f"Extracting {x_tag} and {y_tag} tags from BAM"):
            if read.has_tag(x_tag) and read.has_tag(y_tag):
                x_values.append(read.get_tag(x_tag))
                y_values.append(read.get_tag(y_tag))
    sea.histplot(x=x_values, y=y_values)
    plt.title(f"{x_tag} vs {y_tag}\n[{bam_file.name}]")
    plt.xlabel(x_tag)
    plt.ylabel(y_tag)
    plt.show()

if __name__ == '__main__':
    test_bam = Path("../examples/3RACE/umi_grouping/ont3RACE_PCR8.sorted.calmd.tagged.failed.bam")
    tags_to_test = [
        # "um",
        "ul",
        # "ud",
        # "ui",
    ]
    for tag in tags_to_test:
        histogram_tag_values(test_bam, tag)
    # scatter_tag_v_tag(test_bam, "ui", "ud")
