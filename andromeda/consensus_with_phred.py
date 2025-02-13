"""
consensus_with_phred.py
Marcus Viscardi,    February 6, 2025



"""
import argparse
import pysam
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from tqdm.auto import tqdm
from typing import Dict, List, Tuple

from andromeda.alignment_tools import extract_ref_and_query_region

import andromeda.phred_tools as pT


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
    return str(target_record.seq)


def extract_umi_groups(bam_file: Path, min_group_size: int = 2) -> Dict[str, List[pysam.AlignedSegment]]:
    """
    Groups BAM reads by UMI (UG tag).

    Args:
        bam_file (Path): Path to BAM file.
        min_group_size (int): Minimum reads per UMI group.

    Returns:
        Dict[str, List[pysam.AlignedSegment]]: Dictionary {UG_tag: [reads]}
    """
    umi_groups = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in tqdm(bam, desc="📍 Extracting UMI groups", total=bam.mapped):
            if not read.has_tag("UG"):
                continue
            umi_group_id = read.get_tag("UG")
            if umi_group_id not in umi_groups:
                umi_groups[umi_group_id] = []
            umi_groups[umi_group_id].append(read)

    # Filter groups by size
    umi_groups = {k: v for k, v in umi_groups.items() if len(v) >= min_group_size}
    print(f"✅ {len(umi_groups)} UMI groups retained (min size: {min_group_size})")
    return umi_groups


def extract_ref_oriented_sequence(read_alignment: pysam.AlignedSegment,
                                  reference_seq: str) -> dict:
    """Aligns a query/read sequence and its phred scores to the reference.

    Args:
        read_alignment (pysam.AlignedSegment): A read alignment object.
        reference_seq (str): The reference sequence.

    Returns:
        dict: Dictionary with ref and query positions and sequences, all aligned to the reference.
            ref_positions: List of reference positions.
            ref_sequence: List of reference bases.
            query_positions: List of query positions (aligned to reference, from ref_start position).
            query_sequence: List of query bases.
            query_phreds: List of query phred scores.
    """
    real_ref_seq = reference_seq
    aligned_pairs = read_alignment.get_aligned_pairs(with_seq=False)  # This will just fail, whatever
    region_ref_positions = [ref for query, ref in aligned_pairs if ref is not None]
    region_ref_sequence = [real_ref_seq[i] for i in region_ref_positions]
    # Now we need to be able to pull out the same nucleotides for the actual query sequence
    region_query_positions = [query for query, ref in aligned_pairs if ref is not None]
    region_query_sequence = [read_alignment.query_sequence[i] if i else "." for i in region_query_positions]
    if read_alignment.reference_start >= 5:
        print("potentially interesting region")
    output_dict = {
        "ref_positions": region_ref_positions,
        "ref_sequence": region_ref_sequence,
        "query_positions": region_query_positions,
        "query_sequence": region_query_sequence,
    }
    return output_dict


def extract_ref_oriented_sequences(reads: List[pysam.AlignedSegment], reference_seq: str) -> List[str]:
    """Aligns reads to reference orientation."""
    # TODO: THIS IS BROKEN BECAUSE WE HAVE DELETIONS BEING RETAINED!!
    aligned_seqs = []
    for read in reads:
        query_seq_ref_oriented = extract_ref_oriented_sequence(read, reference_seq)["query_sequence"]
        cons_len_string = (f"{read.reference_start * ' '}"
                           f"{query_seq_ref_oriented}"
                           f"{(len(reference_seq) - read.reference_end) * ' '}")
        aligned_seqs.append(cons_len_string)
    return aligned_seqs


def extract_ref_oriented_sequence_with_phred(read_alignment: pysam.AlignedSegment,
                                             reference_seq: str) -> dict:
    """Aligns a query/read sequence and its phred scores to the reference.
    
    Args:
        read_alignment (pysam.AlignedSegment): A read alignment object.
        reference_seq (str): The reference sequence.
    
    Returns:
        dict: Dictionary with ref and query positions and sequences, all aligned to the reference.
            ref_positions: List of reference positions.
            ref_sequence: List of reference bases.
            query_positions: List of query positions (aligned to reference, from ref_start position).
            query_sequence: List of query bases.
            query_phreds: List of query phred scores.
    """
    real_ref_seq = reference_seq
    aligned_pairs = read_alignment.get_aligned_pairs(with_seq=False)  # This will just fail, whatever
    region_ref_positions = [ref for query, ref in aligned_pairs if ref is not None]
    region_ref_sequence = [real_ref_seq[i] for i in region_ref_positions]
    # Now we need to be able to pull out the same nucleotides for the actual query sequence
    region_query_positions = [query for query, ref in aligned_pairs if ref is not None]
    region_query_sequence = [read_alignment.query_sequence[i] if i else "." for i in region_query_positions]
    region_query_phreds = [pT.NucleotideQuality(q_score=read_alignment.query_qualities[i]).to_phred_char()
                           if i else " " for i in region_query_positions]
    if read_alignment.reference_start >= 5:
        print("potentially interesting region")
    output_dict = {
        "ref_positions": region_ref_positions,
        "ref_sequence": region_ref_sequence,
        "query_positions": region_query_positions,
        "query_sequence": region_query_sequence,
        "query_phreds": region_query_phreds,
    }
    return output_dict


def extract_ref_oriented_sequences_with_phreds(reads: List[pysam.AlignedSegment], reference_seq: str) -> List[dict]:
    """Aligns reads to reference orientation.
    
    Args:
        reads (List[pysam.AlignedSegment]): List of read alignments.
        reference_seq (str): Reference sequence.
    
    Returns:
        List[dict]: List of dictionaries with query and phred sequences aligned to the reference.
            query_positions: List of query positions (aligned to reference, from ref_start position).
            query_sequence: List of query bases.
    """
    aligned_seqs = []
    for read in reads:
        query_seq_ref_oriented = extract_ref_oriented_sequence_with_phred(read, reference_seq)
        const_len_seq = "".join([
            " " * read.reference_start,
            "".join(query_seq_ref_oriented["query_sequence"]),
            " " * (len(reference_seq) - read.reference_end),
        ])
        const_len_phreds = "".join([
            " " * read.reference_start,
            "".join(query_seq_ref_oriented["query_phreds"]),
            " " * (len(reference_seq) - read.reference_end),
        ])
        save_dict = {
            "query_sequence": const_len_seq,
            "query_phreds": const_len_phreds,
        }
        aligned_seqs.append(save_dict)
    return aligned_seqs


def pick_major_base(column: List[str]) -> Tuple[str, float]:
    base_counts = pd.Series(column).value_counts()
    total = base_counts.sum()
    try:
        major_base = base_counts.idxmax()  # can throw error if we run out of sequence
    except ValueError:
        major_base = " "
    major_fraction = base_counts.max() / total if total > 0 else 0
    return major_base, major_fraction


def pick_major_base_and_phred(column: List[Tuple[str, str]]) -> Tuple[str, str]:
    base_counts = pd.Series([base[0] for base in column]).value_counts()
    total = base_counts.sum()
    try:
        major_base = base_counts.idxmax()  # can throw error if we run out of sequence
    except ValueError:
        major_base = " "
    major_fraction = base_counts.max() / total if total > 0 else 0
    major_phred_scores = [base[1] for base in column if base[0] == major_base]
    avg_phred = pT.avg_phred_char_strings(major_phred_scores)
    avg_phred *= (1 / major_fraction)  # TODO: Try different methods here!!
    return major_base, avg_phred, major_fraction


def compute_majority_consensus(reads: List[pysam.AlignedSegment], reference_seq: str) -> Tuple[str, Dict]:
    """
    Computes a majority-rule consensus sequence for a UMI group.

    Args:
        reads (List[pysam.AlignedSegment]): Reads in the UMI group.
        reference_seq (str): The reference sequence.

    Returns:
        Tuple[str, Dict]: Consensus sequence and stats dictionary.
    """
    # TODO: work on this function to add phred averaging
    aligned_seqs = extract_ref_oriented_sequences(reads, reference_seq)

    output_seq = []
    output_phreds = []
    output_key = []
    stats = {
        "matches": 0,
        "mismatches": 0,
        "gaps": 0,
        "low_confidence": 0,
        "mid_confidence": 0,
        "group_size": len(reads)}

    for i, ref_base in enumerate(reference_seq.upper()):
        column = [seq[i].upper() for seq in aligned_seqs if i < len(seq)]

        major_base, major_fraction = pick_major_base(column)

        if major_base == ref_base:
            stats["matches"] += 1
            output_key.append("M")
        elif major_base == " ":
            stats["gaps"] += 1
            output_key.append(" ")
        else:
            stats["mismatches"] += 1
            output_key.append("X")

        if major_fraction < 0.5:
            stats["low_confidence"] += 1
        elif major_fraction < 0.85:
            stats["mid_confidence"] += 1

        output_seq.append(major_base)

    return "".join(output_seq), "".join(output_phreds), "".join(output_key), stats


def compute_majority_consensus_with_phred(reads: List[pysam.AlignedSegment],
                                          reference_seq: str) -> Tuple[str, str, Dict]:
    """
    Computes a majority-rule consensus sequence for a UMI group.
    Args:
        reads (List[pysam.AlignedSegment]): Reads in the UMI group.
        reference_seq (str): The reference sequence.
    Returns:
        Tuple[str, str, str, Dict]: Consensus sequence and stats dictionary.
            index 0: Consensus sequence
            index 1: Consensus phred scores
            index 2: Consensus key (M, X, or space)
            index 3: Stats dictionary
    """  # TODO: work on this function to add phred averaging
    aligned_seqs_and_phreds = extract_ref_oriented_sequences_with_phreds(reads, reference_seq)
    aligned_seqs = [seq["query_sequence"] for seq in aligned_seqs_and_phreds]
    aligned_phreds = [seq["query_phreds"] for seq in aligned_seqs_and_phreds]

    output_seq = []
    output_phreds = []
    output_key = []
    stats = {"matches": 0,
             "mismatches": 0,
             "gaps": 0,
             "low_confidence": 0,
             "mid_confidence": 0,
             "group_size": len(reads)}

    for i, ref_base in enumerate(reference_seq.upper()):
        seq_char_column = [seq[i].upper() for seq in aligned_seqs if i < len(seq)]
        phred_char_column = [phred[i] for phred in aligned_phreds if i < len(phred)]
        column = list(zip(seq_char_column, phred_char_column))

        major_base, avg_phred, major_fraction = pick_major_base_and_phred(column)

        if major_base == ref_base:
            stats["matches"] += 1
            output_key.append("M")
        elif major_base == " ":
            stats["gaps"] += 1
            output_key.append(" ")
        else:
            stats["mismatches"] += 1
            output_key.append("X")

        if major_fraction < 0.5:
            stats["low_confidence"] += 1
        elif major_fraction < 0.85:
            stats["mid_confidence"] += 1

        output_seq.append(major_base)
        output_phreds.append(avg_phred)

    return "".join(output_seq), "".join(output_phreds), "".join(output_key), stats


def call_consensus_from_bam(bam_file: Path, reference_fasta: Path, output_dir: Path, min_group_size: int = 2):
    """
    Main function to extract UMIs, compute consensus sequences, and save results.

    Args:
        bam_file (Path): Input BAM file with UMI groups.
        reference_fasta (Path): Reference genome.
        output_dir (Path): Directory for output.
        min_group_size (int): Minimum reads per UMI group.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print("📂 Loading reference...")
    reference_seq = load_reference(reference_fasta)

    print("📍 Extracting UMI groups...")
    umi_groups = extract_umi_groups(bam_file, min_group_size)

    results = []
    for umi_id, reads in tqdm(umi_groups.items(), desc="🔬 Calling consensus"):
        consensus_seq, consensus_phred, match_str, stats = compute_majority_consensus_with_phred(reads, reference_seq)
        results.append({"UMI": umi_id, "Consensus": consensus_seq, **stats})

    df = pd.DataFrame(results)
    df.to_csv(output_dir / "consensus_sequences.tsv", sep="\t", index=False)
    print(f"✅ Saved consensus sequences to {output_dir}/consensus_sequences.tsv")

    # plot_consensus_stats(df, output_dir)
    return df


def plot_consensus_stats(df: pd.DataFrame, output_dir: Path):
    """Generates a plot showing mismatch and confidence rates."""

    # TODO: I'd rather plot our old thing we had where we showed the average accuracy of the consensus sequences
    #       and the number of UMI groups at each membership number.

    plt.figure(figsize=(6, 4))
    sns.histplot(df["mismatches"], bins=20, kde=True, color="red", alpha=0.6, label="Mismatches")
    sns.histplot(df["low_confidence"], bins=20, kde=True, color="purple", alpha=0.6, label="Low Confidence")
    plt.xlabel("Number of Sites")
    plt.ylabel("Frequency")
    plt.legend()
    plt.title("Consensus Call Statistics")
    plt.grid(True)

    output_plot = output_dir / "consensus_quality.png"
    plt.savefig(output_plot, dpi=300)
    print(f"📊 Saved plot: {output_plot}")
    plt.close()


def call_consensus_and_plot(args: argparse.Namespace):
    call_consensus_from_bam(
        bam_file=args.bam_file,
        reference_fasta=args.reference_fasta,
        output_dir=args.output_dir,
        min_group_size=args.min_group_size,
    )  # TODO: Fix the save path to have a better name (more specific)
    print("📊 Plotting consensus quality...")
    df = pd.read_csv(args.output_dir / "consensus_sequences.tsv", sep="\t")
    if args.plot:
        plot_consensus_stats(df, args.output_dir)


def parse_args():
    parser = argparse.ArgumentParser(description="Compute consensus sequences for UMI groups in BAM files.")

    parser.add_argument("bam_file", type=Path, help="Path to input BAM file (must have UMIs tagged).")
    parser.add_argument("reference_fasta", type=Path, help="Path to reference FASTA file.")
    parser.add_argument("output_dir", type=Path, help="Directory to save consensus sequences.")
    parser.add_argument("--min-group-size", type=int, default=2, help="Minimum reads per UMI group.")
    parser.add_argument("--plot", action="store_true", help="Plot consensus quality.")

    return parser


def main():
    args = parse_args().parse_args()
    df = call_consensus_from_bam(
        bam_file=args.bam_file,
        reference_fasta=args.reference_fasta,
        output_dir=args.output_dir,
        min_group_size=args.min_group_size,
    )
    if args.plot:
        plot_consensus_stats(df, args.output_dir)


if __name__ == "__main__":
    main()
