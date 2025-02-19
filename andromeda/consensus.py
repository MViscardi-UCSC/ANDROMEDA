"""
consensus.py
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
                                  reference_seq: str, extract_phred=True) -> dict:
    """Aligns a query/read sequence and its phred scores to the reference.
    
    Args:
        read_alignment (pysam.AlignedSegment): A read alignment object.
        reference_seq (str): The reference sequence.
        extract_phred (bool): If True, extract phred scores too
    Returns:
        dict: Dictionary with ref and query positions and sequences, all aligned to the reference.
            ref_positions: List of reference positions.
            ref_sequence: List of reference bases.
            query_positions: List of query positions (aligned to reference, from ref_start position).
            query_sequence: List of query bases.
            query_phreds: (optional) List of query phred scores.
    """
    real_ref_seq = reference_seq
    aligned_pairs = read_alignment.get_aligned_pairs(with_seq=False)  # This will just fail, whatever
    region_ref_positions = [ref for query, ref in aligned_pairs if ref is not None]
    region_ref_sequence = [real_ref_seq[i] for i in region_ref_positions]
    # Now we need to be able to pull out the same nucleotides for the actual query sequence
    region_query_positions = [query for query, ref in aligned_pairs if ref is not None]
    region_query_sequence = [read_alignment.query_sequence[i] if i else "." for i in region_query_positions]
    if extract_phred:
        region_query_phreds = [pT.NucleotideQuality(q_score=read_alignment.query_qualities[i]).to_phred_char()
                               if i else " " for i in region_query_positions]
    output_dict = {
        "ref_positions": region_ref_positions,
        "ref_sequence": region_ref_sequence,
        "query_positions": region_query_positions,
        "query_sequence": region_query_sequence,
    }
    if extract_phred:
        output_dict["query_phreds"] = region_query_phreds
    return output_dict


def extract_ref_oriented_sequences(reads: List[pysam.AlignedSegment], reference_seq: str,
                                   extract_phreds=True) -> List[dict]:
    """Aligns reads to reference orientation.
    
    Args:
        reads (List[pysam.AlignedSegment]): List of read alignments.
        reference_seq (str): Reference sequence.
        extract_phreds (bool): If True, extract phred scores too.
    Returns:
        List[dict]: List of dictionaries with query and phred sequences aligned to the reference.
            query_positions: List of query positions (aligned to reference, from ref_start position).
            query_sequence: List of query bases.
    """
    aligned_seqs = []
    for read in reads:
        query_seq_ref_oriented = extract_ref_oriented_sequence(read, reference_seq, extract_phred=extract_phreds)
        const_len_seq = "".join([
            " " * read.reference_start,
            "".join(query_seq_ref_oriented["query_sequence"]),
            " " * (len(reference_seq) - read.reference_end),
        ])
        if extract_phreds:
            const_len_phreds = "".join([
                " " * read.reference_start,
                "".join(query_seq_ref_oriented["query_phreds"]),
                " " * (len(reference_seq) - read.reference_end),
            ])
        save_dict = {
            "query_sequence": const_len_seq,
        }
        if extract_phreds:
            save_dict["query_phreds"] = const_len_phreds
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
    avg_phred = pT.avg_phred_char_strings_to_obj(major_phred_scores)
    # if avg_phred.to_phred_char() != " " and avg_phred.to_prob_error() > 0.15:
    #     print("Woah, high error rate!")
    if avg_phred.to_phred_char() != " ":
        avg_phred *= (1 / major_fraction)  # TODO: Try different methods here!!
        avg_phred = avg_phred.to_phred_char()
    else:
        avg_phred = " "
    return major_base, avg_phred, major_fraction


def compute_majority_consensus(reads: List[pysam.AlignedSegment],
                               reference_seq: str, calc_avg_phreds=True) -> Tuple[str, str, str, Dict]:
    """
    Computes a majority-rule consensus sequence for a UMI group.
    Args:
        reads (List[pysam.AlignedSegment]): Reads in the UMI group.
        reference_seq (str): The reference sequence.
        calc_avg_phreds (bool): If True, calculate average phred scores for each base.
    Returns:
        Tuple[str, str, str, Dict]: Consensus sequence and stats dictionary.
            index 0: Consensus sequence
            index 1: Consensus phred scores
            index 2: Consensus key (M, X, or space)
            index 3: Stats dictionary
    """
    aligned_seqs_and_phreds = extract_ref_oriented_sequences(reads, reference_seq, extract_phreds=calc_avg_phreds)
    aligned_seqs = [seq["query_sequence"] for seq in aligned_seqs_and_phreds]
    if calc_avg_phreds:
        aligned_phreds = [seq["query_phreds"] for seq in aligned_seqs_and_phreds]
        assert len(aligned_seqs) == len(aligned_phreds), "Mismatched lengths of sequences and phreds!"
    
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
        seq_char_column = [seq[i].upper() for seq in aligned_seqs]
        if calc_avg_phreds:
            phred_char_column = [phred[i] for phred in aligned_phreds]
            column = list(zip(seq_char_column, phred_char_column))
            major_base, avg_phred, major_fraction = pick_major_base_and_phred(column)
        else:
            major_base, major_fraction = pick_major_base(seq_char_column)
            avg_phred = "!"

        if major_base == ref_base:
            stats["matches"] += 1
            output_key.append("M")
        elif major_base == " ":
            stats["gaps"] += 1
            output_key.append(" ")
        elif major_base == ".":
            stats["gaps"] += 1
            output_key.append("D")
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


def create_consensus_cigar(match_str, collapse_X_to_M=False) -> Tuple[str, int]:
    """
    Creates a CIGAR string for a consensus sequence and extracts the start position.
    :param match_str: A string of cigar-like characters (M, X, D, or space)
    :param collapse_X_to_M: If True, X characters will be converted to M characters
    :return: Tuple of CIGAR string and start position
    """
    # The match_str is a long list of values (M, X, D, or space) that we can use to build a CIGAR string
    if collapse_X_to_M:
        match_str = match_str.replace("X", "M")
    cigared_str = []
    current_char_count = 0
    current_char = None
    for char in match_str:
        if char == current_char:
            current_char_count += 1
        else:
            if current_char is not None and current_char != " ":
                cigared_str.append(f"{current_char_count}{current_char}")
            current_char = char
            current_char_count = 1
    if current_char is not None and current_char_count > 0 and current_char != " ":
        cigared_str.append(f"{current_char_count}{current_char}")
    return "".join(cigared_str), match_str.index("M")


def call_consensus_from_bam(bam_file: Path,
                            reference_fasta: Path,
                            output_dir: Path,
                            mismatches_in_cigar: bool = False,
                            calc_avg_phreds: bool = True,
                            min_group_size: int = 2) -> pd.DataFrame:
    """
    Main function to extract UMIs, compute consensus sequences, and save results.

    Args:
        bam_file (Path): Input BAM file with UMI groups.
        reference_fasta (Path): Reference genome.
        output_dir (Path): Directory for output.
        mismatches_in_cigar (bool): If True, mismatches (X) will be included in the CIGAR string.
        calc_avg_phreds (bool): If True, average phred scores will be calculated for each base.
        min_group_size (int): Minimum reads per UMI group.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print("📂 Loading reference...")
    reference_seq = load_reference(reference_fasta)

    print("📍 Extracting UMI groups...")
    umi_groups = extract_umi_groups(bam_file, min_group_size)
    
    if len(umi_groups) == 0:
        print("")
        raise ValueError("❌ No UMI groups found matching provided cutoffs! Exiting.")

    results = []
    consensus_iterator = tqdm(umi_groups.items(), desc="🔬 Calling consensus")
    for umi_id, reads in consensus_iterator:
        # consensus_iterator.set_description(f"🔬 Calling consensus (uID: {umi_id}, #: {len(reads)})")
        consensus_iterator.set_postfix_str(f"Current uID: {umi_id:0>5}; Members: {len(reads):>4}")
        consensus_seq, consensus_phred, match_str, stats = compute_majority_consensus(reads, reference_seq,
                                                                                      calc_avg_phreds=calc_avg_phreds)
        consensus_cigar, consensus_start = create_consensus_cigar(match_str,
                                                                  collapse_X_to_M=mismatches_in_cigar)
        results.append({
            "UMI": umi_id,
            "Consensus": consensus_seq,
            "CIGAR": consensus_cigar,
            "start_pos": consensus_start,
            **stats})

    df = pd.DataFrame(results)
    output_file_name = bam_file.stem + "_consensus_sequences.tsv"
    df.to_csv(output_dir / output_file_name, sep="\t", index=False)
    print(f"✅ Saved consensus sequences to {output_dir}/{output_file_name}")

    # plot_consensus_stats(df, output_dir)
    return df


def concensus_df_to_bam(df: pd.DataFrame, output_bam: Path, reference_fasta: Path) -> Path:
    """
    Converts a DataFrame of consensus sequences to a BAM file.
    Args:
        df (pd.DataFrame): DataFrame with columns: UMI, Consensus, CIGAR, start_pos.
        output_bam (Path): Path to save the output BAM file.
        reference_fasta (Path): Path to reference FASTA file.
    """
    with pysam.FastaFile(reference_fasta) as ref:
        header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": ref.get_reference_length(ref.references[0]), "SN": ref.references[0]}]}
        with pysam.AlignmentFile(output_bam, "wb", header=header) as outf:
            for _, row in df.iterrows():
                read = pysam.AlignedSegment()
                read.query_name = row["UMI"]
                read.query_sequence = row["Consensus"]
                read.cigarstring = row["CIGAR"]
                read.reference_start = row["start_pos"]
                read.reference_id = 0
                read.flag = 0
                outf.write(read)
    print(f"✅ Saved BAM file: {output_bam}")
    return output_bam


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
    )
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
        mismatches_in_cigar=True,
    )
    output_bam_name = args.bam_file.stem + "_consensus_sequences.bam"
    bam_path = concensus_df_to_bam(df, args.output_dir / output_bam_name, args.reference_fasta)
    if args.plot:
        plot_consensus_stats(df, args.output_dir)


if __name__ == "__main__":
    main()
