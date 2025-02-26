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
from andromeda.io_utils import load_reference

HELP_TEXT = f"Extract UMIs from mapped reads (from a BAM file) based on aligned reference positions."


def load_umi_positions(umi_positions_file: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Loads UMI positions from a TSV file and converts them into a dictionary."""
    contig_df = pd.read_csv(umi_positions_file)
    # Columns = "ref_contig", "start", "end", "date_time"

    output_dict = {}  # Dict of umi_positions per contig
    for index, row in contig_df.set_index("ref_contig").iterrows():
        print(f"    📍 Region {index}: {row['start']}-{row['end']}")
        if index not in output_dict:
            output_dict[index] = [(0, (row["start"], row["end"]))]
        else:
            prev_index = output_dict[index][-1][0]
            output_dict[index].append((prev_index + 1, (row["start"], row["end"])))
    return output_dict


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
):
    """Extracts UMIs from a BAM file and tags reads with UMI sequences."""

    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n📂 Loading reference sequence...")
    contig, ref_seq = load_reference(reference_file)

    print("\n🔍 Loading UMI positions...")
    contig_umi_positions = load_umi_positions(umi_positions_file)

    if contig not in contig_umi_positions:
        raise KeyError(f"❌ Contig '{contig}' not found in {umi_positions_file}")

    # TODO: Ability to check multiple contigs one by one.
    umi_positions_list = contig_umi_positions[contig]
    print(f"✅ Found {len(umi_positions_list)} UMI position sets for contig {contig}")

    output_bams = []

    print("\n🛠 Processing BAM file...")
    umi_positions_with_flanking = []
    for umi_index, umi_positions in umi_positions_list:
        umi_positions_with_flanking.append((umi_index,
                                            )
                                           )
    if len(umi_positions_list) > 1:
        print(f"📝 Tagging multiple UMIs from the same contig ({contig})...")
        multiple_umis_per_contig = True
    else:
        print(f"📝 Tagging single UMI from contig {contig}...")
        multiple_umis_per_contig = False
    for umi_index, umi_positions in umi_positions_list:
        flank_adjusted_positions = (umi_positions[0] - flanking_seq_to_capture,
                                    umi_positions[1] + flanking_seq_to_capture)
        print(f"    📍 Ready to extract UMI at position {umi_index + 1}: ({umi_positions[0]}-{flanking_seq_to_capture},"
              f" {umi_positions[1]}+{flanking_seq_to_capture})")
        print(f"    🔍 Sequence: "
              f"{ref_seq[flank_adjusted_positions[0]:flank_adjusted_positions[1]+1]}")
        if not force:
            confirm_selection = input("    Confirm selection? (y/N): ").strip().lower()
            if confirm_selection != "y":
                print("    🔄 Skipping this UMI.")
                continue
            else:
                print("    ✅ UMI confirmed.")

        if multiple_umis_per_contig:
            # If multiple UMIs per contig, add a suffix to the BAM file
            suffix = f".tagged_u{umi_index + 1}.bam"
            # And changed the tags to be numbered in case you want to merge them later
            tag_dict = {
                "umi_sequence": f"u{umi_index + 1}",
                "deletion_count": f"d{umi_index + 1}",
                "insertion_count": f"i{umi_index + 1}",
                "mismatch_count": f"m{umi_index + 1}",
                "umi_length": f"l{umi_index + 1}",
            }
        else:
            suffix = ".tagged.bam"
            tag_dict = {
                "umi_sequence": umi_tag,
                "deletion_count": "ud",
                "insertion_count": "ui",
                "mismatch_count": "um",
                "umi_length": "ul",
            }
            if umi_tag != "uM":
                assert umi_tag not in [v for k, v in tag_dict.items() if k != "umi_sequence"], \
                    (f"❌ UMI tag '{umi_tag}' already in use! "
                     f"Please choose a different tag or use the default 'uM'.")

        tagged_bam = AlignTools.bam_to_tagged_bam(
            bam_file, contig, ref_seq,
            umi_positions[0], umi_positions[1],
            save_dir=output_dir,
            flanking_seq_to_capture=flanking_seq_to_capture,
            bam_tag_dict=tag_dict,
            subset_count=subset_count,
            save_suffix=suffix,
            max_del_in_umi=0,  # Keep this
            max_ins_in_umi=0,  # And keep this, b/c it allows for easier downstream analysis
            max_iupac_mismatches=mismatch_tolerance,  # This could potentially be changed
            restrict_to_length=True,  # The 0 del and 0 ins pretty much forces this already
        )
        output_bams.append(tagged_bam)

        # Index the BAM file after each tagging step
        subprocess.run(["samtools", "index", str(tagged_bam)], check=True)
        # Convert the final BAM file to SAM for easier reading
        subprocess.run(["samtools", "view", str(tagged_bam),
                        "-o", str(tagged_bam.with_suffix(".sam"))], check=True)
        print(f"    ✅ UMI {umi_index + 1} extraction complete! Saved to:")
        print(f"    📁 {tagged_bam}")

    print("\n✅ UMI extraction complete!")
    return output_bams


def save_umi_counts(bam_file: Path, output_dir: Path, umis=None):
    """Extracts UMI sequences and saves frequency counts as TSV & Pickle."""

    print("\n📊 Generating UMI counts...")
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
    series = umi_df[umi_tag_name].value_counts().groupby(level=umi_tag_name).sum().sort_values(ascending=False)
    output_name = umi_tag_name + "_counts"

    print(f"📝 Saving {output_name}...")
    series.to_csv(output_dir / f"{output_name}.tsv", sep="\t")
    with open(output_dir / f"{output_name}.pkl", "wb") as f:
        pickle.dump(series.to_dict(), f)

    print("\n✅ UMI counts saved!")


def parse_args():
    parser = argparse.ArgumentParser(description="Extract UMIs from BAM files based on mapped reference positions.")

    parser.add_argument("ref_fasta", type=Path,
                        help="Path to reference FASTA file.")
    parser.add_argument("mapped_bam", type=Path,
                        help="Path to input BAM file.")
    parser.add_argument("umi_positions", type=Path,
                        help="TSV file with UMI positions per contig. You can use ref-pos-picker to generate this.")
    parser.add_argument("output_parent_dir", type=Path,
                        help="Parent directory to make a new directory inside to save outputs.")
    parser.add_argument("--extraction-flanking", type=int, default=0,
                        help="Flanking bases to capture around UMI positions [default: 0].")
    parser.add_argument("--extraction-subset", type=int, default=-1,
                        help="Subset number of BAM reads for testing (-1 for all) [default: -1].")
    parser.add_argument("--extraction-do-not-confirm", action="store_true",
                        help="Do not confirm UMI selection.")
    parser.add_argument("--extraction-mismatch-tolerance", type=int, default=0,
                        help="Max number allowed IUPAC mismatches in UMI sequence [default: 0].")
    parser.add_argument("--store-umi-tag", type=str, default="uM",
                        help="BAM tag to save UMIs to [default: 'uM'].")

    return parser


def dependencies():
    return {
        "ref_fasta": "ref_pos_picker.ref_fasta",
        "umi_positions": "ref_pos_picker.umi_positions",
        "output_parent_dir": "ref_pos_picker.output_parent_dir",
    }


def pipeline_main(args):
    output_dir = args.output_parent_dir / "tagging"
    output_dir.mkdir(parents=True, exist_ok=True)
    assert args.umi_positions.exists(), \
        (f"UMI positions file not found!\n"
         f"Looked here: {args.umi_positions}")
    tagged_bams = extract_umis_from_bam(
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
    for bam_path in tagged_bams:
        save_umi_counts(bam_path, output_dir)
    pass_fwd_dict = {
        "output_parent_dir": args.output_parent_dir,
        "load_umi_tag": args.store_umi_tag,
    }
    if len(tagged_bams) == 1:
        pass_fwd_dict["tagged_bam"] = tagged_bams[0]
    else:
        raise NotImplementedError("Multiple tagged BAMs not yet supported.")
    return pass_fwd_dict


if __name__ == "__main__":
    args = parse_args().parse_args()
    pipeline_main(args)
    # TODO: Interestingly, most fails are a single nucleotide short of the expected length, due to a single deletion...
    #       Should we try to retain these species? They might be interesting to look at...
