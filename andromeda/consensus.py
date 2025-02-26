"""
consensus.py
Marcus Viscardi,    February 6, 2025


Complute/collapse consensus sequences for UMI groups in BAM files.
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
import sys
import andromeda.phred_tools as pT
from andromeda.io_utils import load_reference

HELP_TEXT = f"Compute/collapse consensus sequences for UMI groups in BAM files."


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
        try:
            region_query_phreds = [pT.NucleotideQuality(q_score=read_alignment.query_qualities[i]).to_phred_char()
                                   if i else " " for i in region_query_positions]
        except TypeError:
            # This is a case where the query quality is None because it wasn't stored in the BAM file?
            # We'll just give every read a "perfect" quality score for now...
            # TODO: Expand this functionality a bit more to handle this case better
            region_query_phreds = ["A" if i else " " for i in region_query_positions]
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


def pick_major_base_and_phred(column: List[Tuple[str, str]]) -> Tuple[str, str, float]:
    base_counts = pd.Series([base[0] for base in column]).value_counts()
    total = base_counts.sum()
    try:
        major_base = base_counts.idxmax()  # can throw error if we run out of sequence
    except ValueError:
        major_base = " "
    major_fraction = base_counts.max() / total if total > 0 else 0
    if major_base == " " or major_base == ".":
        # Pull rip cord if we have a gap or ambiguous base, don't bother with phred scores
        return major_base, " ", major_fraction

    major_phred_scores = [base[1] for base in column if base[0] == major_base]
    avg_phred = pT.avg_phred_char_strings_to_obj(major_phred_scores)
    if avg_phred.to_phred_char() != " ":
        avg_phred *= (1 / major_fraction)  # TODO: Try different methods here!!
        avg_phred_char = avg_phred.to_phred_char()
    else:
        avg_phred_char = " "
    if avg_phred.to_prob_error() > 1:  # This should help with clamping the phred scores
        avg_phred_char = "!"
    return major_base, avg_phred_char, major_fraction


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
            # stats["gaps"] += 1  # Dropped b/c this is just "soft-clipping"
            output_key.append(" ")
        elif major_base == "." and "M" not in output_key:
            # This is just a weird edge case where we have a gap at the start of the sequence
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
                            output_parent_dir: Path,
                            mismatches_in_cigar: bool = True,
                            calc_avg_phreds: bool = True,
                            min_group_size: int = 2) -> Tuple[pd.DataFrame, Path]:
    """
    Main function to extract UMIs, compute consensus sequences, and save results.

    Args:
        bam_file (Path): Input BAM file with UMI groups.
        reference_fasta (Path): Reference genome.
        output_parent_dir (Path): Directory for output.
        mismatches_in_cigar (bool): If True, mismatches (X) will be included in the CIGAR string.
        calc_avg_phreds (bool): If True, average phred scores will be calculated for each base.
        min_group_size (int): Minimum reads per UMI group.
    """
    output_dir = Path(output_parent_dir) / "consensus"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("📂 Loading reference...")
    contig, reference_seq = load_reference(reference_fasta)

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
                                                                  collapse_X_to_M=not mismatches_in_cigar)
        consensus_seq_trimmed = consensus_seq.replace(" ", "").replace(".", "")
        # We probably need to do the same as above with the phred scores...
        consensus_phred_trimmed = consensus_phred.replace(" ", "")  # TODO: Did this work?
        assert len(consensus_seq_trimmed) == len(consensus_phred_trimmed), "Mismatched lengths of sequence and phred!"
        results.append({
            "UMI": umi_id,
            "Consensus": consensus_seq_trimmed,
            "Phred": consensus_phred_trimmed,
            "CIGAR": consensus_cigar,
            "start_pos": consensus_start,
            **stats})

    df = pd.DataFrame(results)
    output_file_name = bam_file.stem + "_consensus_sequences.tsv"
    output_file_path = output_dir / output_file_name
    df.to_csv(output_file_path, sep="\t", index=False)
    print(f"✅ Saved consensus sequences to {output_file_path}")

    # plot_consensus_stats(df, output_dir)
    return df, output_file_path


def consensus_df_to_bam(df: pd.DataFrame, output_bam: Path, reference_fasta: Path, template_bam: Path) -> Path:
    """
    Converts a DataFrame of consensus sequences to a BAM file.
    Args:
        df (pd.DataFrame): DataFrame with columns: UMI, Consensus, CIGAR, start_pos.
        output_bam (Path): Path to save the output BAM file.
        reference_fasta (Path): Path to reference FASTA file.
        template_bam (Path): Path to template BAM file. (This is used to copy the header.)
    """
    with (pysam.FastaFile(reference_fasta) as ref,
          pysam.AlignmentFile(template_bam, "rb") as template):
        header = template.header.to_dict()
        # pprint(header)
        previous_program_pn = header['PG'][-1]['PN']
        add_to_header = {"ID": "ANDROMEDA-consensus",
                         "PN": "ANDROMEDA",
                         "PP": previous_program_pn,  # This needs to be the PN of the previous program!!
                         "VN": "0.0.01",
                          "CL": " ".join(sys.argv)}
        header['PG'].append(add_to_header)
        # pprint(header)
        with pysam.AlignmentFile(output_bam, "wb",
                                 header=header,
                                 # template=template,
                                 ) as outf:
            for _, row in df.iterrows():
                assert len(row["Consensus"]) == len(row["Phred"]), "Mismatched lengths of sequence and phred!"
                read = pysam.AlignedSegment()
                read.query_name = f"UMIID{row['UMI']:0>6}"
                read.query_sequence = row["Consensus"]
                read.flag = 0
                read.cigarstring = row["CIGAR"]
                read.reference_start = row["start_pos"]
                read.reference_id = 0
                read.mapping_quality = 20  # Arbitrary for now
                read.template_length = len(read.query_sequence)
                read.query_qualities = pT.PhredString(phred_string=row["Phred"]).to_q_scores()
                read.set_tag("ug", row["group_size"], "i")
                outf.write(read)
    print(f"✅ Saved BAM file: {output_bam}")
    pysam.sort("-o", str(output_bam.with_suffix(".sorted.bam")), str(output_bam))
    print(f"✅ Sorted BAM file: {output_bam.with_suffix('.sorted.bam')}")
    pysam.index(str(output_bam.with_suffix('.sorted.bam')))
    print(f"✅ Indexed BAM file: {output_bam.with_suffix('.sorted.bam')}.bai")
    output_sam = output_bam.with_suffix('.sorted.bam').with_suffix(".sam")
    output_sam.touch()
    pysam.view("-ho", str(output_sam), str(output_bam.with_suffix('.sorted.bam')),
               save_stdout=str(output_sam))
    print(f"✅ Saved SAM file: {output_bam.with_suffix('.sam')}")
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
    df, tsv_path = call_consensus_from_bam(
        bam_file=args.grouped_bam,
        reference_fasta=args.ref_fasta,
        output_parent_dir=args.output_parent_dir,
        min_group_size=args.min_group_size,
        mismatches_in_cigar=True,
    )
    output_bam_name = args.grouped_bam.stem + ".consensus_sequences.bam"
    bam_path = consensus_df_to_bam(df, args.output_parent_dir / "consensus" / output_bam_name,
                                   args.ref_fasta,
                                   args.grouped_bam)
    if args.consensus_plot:
        plot_consensus_stats(df, args.output_dir)
    return tsv_path, bam_path


def parse_args():
    parser = argparse.ArgumentParser(description="Compute consensus sequences for UMI groups in BAM files.")

    parser.add_argument("ref_fasta", type=Path,
                        help="Path to reference FASTA file.")
    parser.add_argument("grouped_bam", type=Path,
                        help="Path to input BAM file (must have UMIs grouped).")
    parser.add_argument("output_parent_dir", type=Path,
                        help="Parent directory to make a new directory inside to save outputs.")
    parser.add_argument("--min-group-size", type=int, default=2,
                        help="Minimum reads per UMI group [default: 2].")
    parser.add_argument("--consensus-plot", action="store_true",
                        help="Plot consensus quality.")

    return parser


def dependencies():
    return {
        "grouped_bam": "umi_group.grouped_bam",
        "ref_fasta": "ref_pos_picker.ref_fasta",
        "output_parent_dir": "ref_pos_picker.output_dir",
    }


def pipeline_main(args: argparse.Namespace):
    call_consensus_and_plot(args)


if __name__ == "__main__":
    args = parse_args().parse_args()
    pipeline_main(args)
