"""
extract.py
Marcus Viscardi,    January 31, 2025

This script should allow users to extract UMIs from BAM files based on mapped reference positions.

Currently, it works with a single contig but can extract multiple UMIs from different regions of the same contig.
TODO: Add multi-contig support.

TODO: Logging!!
"""

import argparse
from argparse import Namespace

import pysam
import pandas as pd
import pickle
import subprocess
from pathlib import Path
from Bio import SeqIO
from typing import List, Tuple, Dict
import re
from tqdm.auto import tqdm

from andromeda import alignment_tools as AlignTools
from andromeda.io_utils import load_reference_contigs_to_dict, load_umi_positions
from andromeda.logger import log

HELP_TEXT = f"Extract UMIs from mapped reads (from a BAM file) based on aligned reference positions."

KNOWN_SAM_TAGS = [
    "AM",
    "AS",
    "BC",
    "BQ",
    "BZ",
    "CB",
    "CC",
    "CG",
    "CM",
    "CO",
    "CP",
    "CQ",
    "CR",
    "CS",
    "CT",
    "CY",
    "E2",
    "FI",
    "FS",
    "FZ",
    "GC",
    "GQ",
    "GS",
    "H0",
    "H1",
    "H2",
    "HI",
    "IH",
    "LB",
    "MC",
    "MD",
    "MF",
    "MI",
    "ML",
    "MM",
    "MQ",
    "NH",
    "NM",
    "OA",
    "OC",
    "OP",
    "OQ",
    "OX",
    "PG",
    "PQ",
    "PT",
    "PU",
    "Q2",
    "QT",
    "QX",
    "R2",
    "RG",
    "RT",
    "RX",
    "S2",
    "SA",
    "SM",
    "SQ",
    "TC",
    "TS",
    "U2",
    "UQ",
    "X?",
    "Y?",
    "Z?",
]
OUR_UMI_DEFAULT_TAGS = ["ud", "ui", "um", "ul"]


def extract_umis_from_bam(
    bam_file: Path,
    reference_file: Path,
    umi_positions_file: Path,
    output_dir: Path,
    flanking_seq_to_capture: int = 0,
    subset_count: int = -1,
    mismatch_tolerance: int = 0,
    force: bool = False,
    umi_tag: str = "uM",
) -> Path:
    """Extracts UMIs from a BAM file and tags reads with UMI sequences."""
    # TODO: We are not currently handling multiple UMIs per contig, this is a smaller issue.

    assert umi_tag not in OUR_UMI_DEFAULT_TAGS, (
        f"❌ UMI tag '{umi_tag}' already in use by our added tags! "
        f"Please choose a different tag or use the default 'uM'."
    )
    assert umi_tag not in KNOWN_SAM_TAGS, (
        f"❌ UMI tag '{umi_tag}' already in use as a default SAM tag! "
        f"Please choose a different tag or use the default 'uM'."
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("📂 Loading reference sequences...")
    ref_dict = load_reference_contigs_to_dict(reference_file)

    log.info("🔍 Loading UMI positions...")
    contig_umi_positions = load_umi_positions(umi_positions_file)

    # Let's simlify the contig_umi_positions dict to only the first umi_option for each contig
    log.info("🔍 Simplifying UMI positions to only a single option for each contig...")
    simplified_contig_umi_positions = {}
    for contig, umi_positions_list in contig_umi_positions.items():
        if len(umi_positions_list) > 1:
            log.warning(f"🚨 Multiple UMI positions found for contig {contig}!")
            log.warning(f"🚨 Only the first UMI position will be used.")
        simplified_contig_umi_positions[contig] = umi_positions_list[0]
    contig_umi_positions = simplified_contig_umi_positions
    umi_pos_dict = contig_umi_positions

    if flanking_seq_to_capture > 0:
        log.info(
            f"🔍 Capturing flanking sequence of {flanking_seq_to_capture} bases around UMI positions..."
        )
        # We need to update the UMI positions to include the flanking sequence
        for contig, (index, (start, end)) in umi_pos_dict.values():
            umi_pos_dict[contig] = index(
                start - flanking_seq_to_capture, end + flanking_seq_to_capture
            )

    # We should do an "inner merge" of the two dictionaries to ensure that we only process contigs that have both
    # reference sequences and UMI positions
    umi_pos_dict = {
        contig: umi_pos_dict[contig]
        for contig in ref_dict.keys()
        if contig in umi_pos_dict
    }
    ref_dict = {contig: ref_dict[contig] for contig in umi_pos_dict.keys()}
    if not umi_pos_dict or not ref_dict:
        raise ValueError(
            "No contigs found in both the reference file and the UMI positions file!"
        )

    tag_dict = {
        "umi_sequence": umi_tag,
        "deletion_count": "ud",
        "insertion_count": "ui",
        "mismatch_count": "um",
        "umi_length": "ul",
    }

    save_path = output_dir / f"{bam_file.stem}.tagged.bam"

    tagged_bam = AlignTools.bam_to_tagged_bam(
        bam_file_path=bam_file,
        ref_seq_dict=ref_dict,
        umi_dict=umi_pos_dict,
        save_path=save_path,
        bam_tag_dict=tag_dict,
        subset_count=subset_count,
        max_del_in_umi=0,  # Keep this
        max_ins_in_umi=0,  # And keep this, b/c it allows for easier downstream analysis
        max_iupac_mismatches=mismatch_tolerance,  # This could potentially be changed
        restrict_to_length=True,  # The 0 del and 0 ins pretty much forces this already
    )
    # Sort the BAM file after tagging
    log.debug(f"🔃 Sorting tagged BAM file: {tagged_bam}")
    subprocess.run(
        ["samtools", "sort", str(tagged_bam), "-o", str(tagged_bam)], check=True
    )
    # Index the BAM file after sorting step
    log.debug(f"🔃 Indexing sorted BAM file: {tagged_bam}")
    subprocess.run(["samtools", "index", str(tagged_bam)], check=True)
    # Convert the final BAM file to SAM for easier reading
    sorted_tagged_sam = tagged_bam.with_suffix(".sam")
    log.debug(f"🔃 Converting sorted BAM file to SAM: {tagged_bam}")
    subprocess.run(
        ["samtools", "view", str(tagged_bam), "-o", str(sorted_tagged_sam)], check=True
    )
    log.info(f"📁 Finished saving and processing tagged BAM file and its derivatives.")

    log.info("✅ UMI extraction complete!")

    return save_path


def save_umi_counts(bam_file: Path, output_dir: Path, umis=None):
    """Extracts UMI sequences and saves frequency counts as TSV & Pickle."""

    log.info("📊 Generating UMI counts...")
    umi_df = AlignTools.bam_to_df(bam_file)
    umi_df["count"] = 1

    output_dir.mkdir(parents=True, exist_ok=True)
    # .tagged_u{umi_index + 1}.bam
    # Extract UMI index from BAM filename
    match = re.search(r"\.taggedu(\d+)\.bam", str(bam_file))

    if match:
        umi_index = int(match.group(1))  # Extract the number
        umi_tag_name = f"u{umi_index}"
    else:
        umi_tag_name = "uM"  # Default tag
    series = (
        umi_df[umi_tag_name]
        .value_counts()
        .groupby(level=umi_tag_name)
        .sum()
        .sort_values(ascending=False)
    )
    output_name = umi_tag_name + "_counts"

    out_tsv = output_dir / f"{output_name}.tsv"
    out_pickle = output_dir / f"{output_name}.pkl"

    log.debug(f"📝 Saving {out_tsv}...")
    series.to_csv(out_tsv, sep="\t")

    log.debug(f"📝 Saving {out_pickle}...")
    with open(out_pickle, "wb") as f:
        pickle.dump(series.to_dict(), f)

    log.info("✅ UMI counts saved!")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract UMIs from BAM files based on mapped reference positions."
    )

    parser.add_argument("ref_fasta", type=Path, help="Path to reference FASTA file.")
    parser.add_argument("mapped_bam", type=Path, help="Path to input BAM file.")
    parser.add_argument(
        "umi_positions",
        type=Path,
        help="TSV file with UMI positions per contig. You can use ref-pos-picker to generate this.",
    )
    parser.add_argument(
        "output_parent_dir",
        type=Path,
        help="Parent directory to make a new directory inside to save outputs.",
    )
    parser.add_argument(
        "--extraction-flanking",
        type=int,
        default=0,
        help="Flanking bases to capture around UMI positions [default: 0].",
    )
    parser.add_argument(
        "--extraction-subset",
        type=int,
        default=-1,
        help="Subset number of BAM reads for testing (-1 for all) [default: -1].",
    )
    parser.add_argument(
        "--extraction-do-not-confirm",
        action="store_true",
        help="Do not confirm UMI selection.",
    )
    parser.add_argument(
        "--extraction-mismatch-tolerance",
        type=int,
        default=0,
        help="Max number allowed IUPAC mismatches in UMI sequence [default: 0].",
    )
    parser.add_argument(
        "--store-umi-tag",
        type=str,
        default="uM",
        help="BAM tag to save UMIs to [default: 'uM'].",
    )

    return parser


def dependencies():
    return {
        "ref_fasta": "ref_pos_picker.ref_fasta",
        "umi_positions": "ref_pos_picker.umi_positions",
        "output_parent_dir": "ref_pos_picker.output_parent_dir",
    }


@log.catch
def pipeline_main(args):
    output_dir = args.output_parent_dir / "tagging"
    output_dir.mkdir(parents=True, exist_ok=True)
    assert args.umi_positions.exists(), (
        f"UMI positions file not found!\nLooked here: {args.umi_positions}"
    )
    tagged_bam = extract_umis_from_bam(
        bam_file=args.mapped_bam,
        reference_file=args.ref_fasta,
        umi_positions_file=args.umi_positions,
        output_dir=output_dir,
        flanking_seq_to_capture=args.extraction_flanking,
        mismatch_tolerance=args.extraction_mismatch_tolerance,
        subset_count=args.extraction_subset,
        force=args.extraction_do_not_confirm,
        umi_tag=args.store_umi_tag,
    )

    save_umi_counts(tagged_bam, output_dir)

    pass_fwd_dict = {
        "output_parent_dir": args.output_parent_dir,
        "load_umi_tag": args.store_umi_tag,
        "tagged_bam": tagged_bam,
    }
    return pass_fwd_dict


if __name__ == "__main__":
    args = parse_args().parse_args()
    pipeline_main(args)
    # TODO: Interestingly, most fails are a single nucleotide short of the expected length, due to a single deletion...
    #       Should we try to retain these species? They might be interesting to look at...
