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

def load_reference(reference_path: str | Path, contig: str = None):
    with open(reference_path, "r") as handle:
        records = [record for record in SeqIO.parse(handle, "fasta")]
    record_ids = [record.id for record in records]
    records_dict = {record.id: record for record in records}
    if len(records) > 1:
        if contig is None:
            print(f"Found {len(records)} records in reference file, "
                  f"because you didn't specify a contig we will use "
                  f"the first one: {record_ids[0]}")
            target_record = records[0]
        elif contig not in record_ids:
            raise ValueError(f"Contig {contig} not found in reference file, "
                             f"available contigs: {record_ids}")
        else:
            target_record = records_dict[contig]
    else:
        if contig is not None and contig not in record_ids:
            print(f"Found 1 record in reference file ({record_ids[0]}), but you specified a contig: {contig}."
                  f"We will use the record in the file.")
        elif contig is not None and contig in record_ids:
            print(f"Found 1 record in reference file (matching your request for {contig}), using it.")
        target_record = records[0]
    return str(target_record.id), str(target_record.seq)


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
              f"{ref_seq[flank_adjusted_positions[0]:flank_adjusted_positions[1]]}")
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
                "umi_sequence": "uM",
                "deletion_count": "ud",
                "insertion_count": "ui",
                "mismatch_count": "um",
                "umi_length": "ul",
            }

        tagged_bam = AlignTools.bam_to_tagged_bam(
            bam_file, contig, ref_seq,
            umi_positions[0], umi_positions[1],
            save_dir=output_dir,
            flanking_seq_to_capture=flanking_seq_to_capture,
            bam_tag_dict=tag_dict,
            subset_count=subset_count,
            save_suffix=suffix,
            max_del_in_umi=0,  # Keep this
            max_ins_in_umi=0,  # And keep this, b/c is allows for easier downstream analysis
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

    parser.add_argument("reference_fasta", type=Path,
                        help="Path to reference FASTA file.")
    parser.add_argument("bam_file", type=Path,
                        help="Path to input BAM file.")
    # parser.add_argument("umi_positions", type=Path,
    #                     help="TSV file with UMI positions per contig. You can use ref-pos-picker to generate this.")
    parser.add_argument("output_dir", type=Path,
                        help="Output directory to save results.")
    parser.add_argument("--flanking", type=int, default=0,
                        help="Flanking bases to capture around UMI positions.")
    parser.add_argument("--subset", type=int, default=-1,
                        help="Subset BAM reads for testing (-1 for all).")
    parser.add_argument("--do_not_confirm", action="store_true",
                        help="Do not confirm UMI selection.")
    parser.add_argument("-m", "--mismatch_tolerance", type=int, default=0,
                        help="Max number allowed IUPAC mismatches in UMI sequence.")

    return parser


def extract_umis_and_summarize(args):
    assert args.reference_fasta.with_suffix(".fasta.targetUMIs.csv").exists(), ("UMI positions file not found!\n"
                                                                                f"Looked here: {args.reference_fasta.with_suffix('.fasta.targetUMIs.csv')}")
    tagged_bams = extract_umis_from_bam(
        bam_file=args.bam_file,
        reference_file=args.reference_fasta,
        umi_positions_file=args.reference_fasta.with_suffix(".fasta.targetUMIs.csv"),
        output_dir=args.output_dir,
        flanking_seq_to_capture=args.flanking,
        mismatch_tolerance=args.mismatch_tolerance,
        subset_count=args.subset,
        force=args.do_not_confirm,
    )
    for bam_path in tagged_bams:
        save_umi_counts(bam_path, args.output_dir)

def main(override_args=None):
    args = override_args or parse_args().parse_args()
    extract_umis_and_summarize(args)


if __name__ == "__main__":
    print(Path.cwd())
    overiding_args = Namespace(
        bam_file=Path("../examples/3RACE/mapped_reads/ont3RACE_PCR8.sorted.calmd.bam"),
        reference_fasta=Path("../examples/3RACE/references/ref.3RACE.fasta"),
        # umi_positions=Path("../examples/3RACE/references/ref.3RACE.fasta.targetUMIs.csv"),
        output_dir=Path("../examples/3RACE/umi_tagging"),
        flanking=5,
        mismatch_tolerance=2,
        subset=-1,
        do_not_confirm=True,
    )
    main(override_args=overiding_args)
