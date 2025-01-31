"""
ref_pos_picker.py
Marcus Viscardi,    January 31, 2025

This script should let users pick the reference positions they are interested in for UMI extraction
by giving them a visual representation of the reference sequence and allowing them to pick the start
and end positions around ambiguous bases.

General plan:
Inputs: Reference sequence and (optionally) a contig/chr name (we can just assume the first contig if not given)
Outputs: Start and end positions for the UMI extraction

1. Load the reference sequence and specific contig if given
2. Identify ambiguous bases and their positions
3. Pick a window around the ambiguous bases to display
4. Display the chosen window with the ambiguous bases highlighted and reference positions labeled (every 5?)

"""
from datetime import datetime

from Bio import SeqIO
from pathlib import Path
from typing import List, Tuple
import argparse


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


def identify_ambiguous_bases(reference_seq: str):
    ambiguous_bases = "RYKMSWBDHVN"  # IUPAC codes for ambiguous bases
    ambiguous_positions = []
    for i, base in enumerate(reference_seq):
        if base.upper() in ambiguous_bases:
            ambiguous_positions.append(i)
    return ambiguous_positions


def group_ambiguous_bases(ambiguous_positions: List[int], padding: int = 5) -> List[Tuple[int, int]]:
    """
    Groups ambiguous base positions into contiguous regions, allowing for a padding of X bases around each cluster.
    
    Thanks, ChatGPT!!
    
    :param ambiguous_positions: List of positions (integers) where ambiguous bases are found.
    :param padding: Number of additional bases to include around each grouping.
    :return: List of tuples (start, end) representing grouped ambiguous base regions.
    """
    if not ambiguous_positions:
        return []

    # Sort positions to process in order
    ambiguous_positions.sort()

    grouped_regions = []
    start = ambiguous_positions[0]
    prev_position = start

    for i in range(1, len(ambiguous_positions)):
        current = ambiguous_positions[i]

        # If the current position is more than (padding * 2) bases away, start a new group
        if current > prev_position + (padding * 2):
            # Append the previous group with padding
            grouped_regions.append((max(0, start - padding), prev_position + padding))
            start = current  # Start new group

        prev_position = current  # Update previous position

    # Append the last detected group
    grouped_regions.append((max(0, start - padding), prev_position + padding))

    return grouped_regions


def extract_grouped_sequences(reference_seq: str, padding: int = 5) -> List[Tuple[int, int, str]]:
    ambiguous_positions = identify_ambiguous_bases(reference_seq)
    grouped_regions = group_ambiguous_bases(ambiguous_positions, padding=padding)
    extracted_sequences = []
    for start, end in grouped_regions:
        extracted_sequences.append((start, end, reference_seq[start:end]))
    return extracted_sequences


def print_ambiguous_regions(grouped_seq_and_regions: List[Tuple[int, int, str]],
                            label_every=10):
    """Prints extracted ambiguous regions in a readable, compact format with position labels."""
    print("\n🔬 Identified ambiguous regions:\n")
    for i, (start, end, extracted_seq) in enumerate(grouped_seq_and_regions):
        # add a space between every base for better readability
        formatted_seq = " ".join(extracted_seq)
        formatted_seq = formatted_seq.replace(" ", "")  # get rid of this later if wanted
        # Create position labels every X bases
        position_labels = []
        tick_labels = []
        tick_icon = "|"
        for idx, pos in enumerate(range(start, end)):
            if idx % label_every == 0:
                position_labels.append(f"{pos:<{label_every}}")
                tick_labels.append(f"{tick_icon:<{label_every}}")
        print(f"📍 Region {i + 1}: {start}-{end}")
        print("".join(position_labels))
        print("".join(tick_labels))
        print(formatted_seq + "\n")


def get_user_selection(ref_seq: str) -> Tuple[int, int]:
    """Prompts the user for manual start and end selection."""
    while True:
        try:
            start = int(input("\nEnter the start coordinate: "))
            end = int(input("Enter the end coordinate: "))

            if start >= end:
                print("❌ Start must be less than End. Try again.")
                continue
            if start < 0 or end > len(ref_seq):
                print(f"❌ Coordinates must be within the reference sequence length ({len(ref_seq)}). Try again.")
                continue

            print(f"\n✅ Selected UMI Region: {start} - {end}")
            print(f"🔍 Sequence: {ref_seq[start:end]}")
            confirm = input("Confirm selection? (Y/N): ").strip().lower()

            if confirm == "y":
                return start, end
            else:
                print("🔄 Let's try again.")
        except ValueError:
            print("❌ Invalid input. Please enter numeric values.")


def run_umi_region_picker(reference_fasta_path: str, contig: str = None, padding: int = 5):
    """Runs the CLI-based UMI region selection."""
    reference_fasta_PATH = Path(reference_fasta_path)
    ref_contig, ref_seq = load_reference(reference_fasta_path, contig=contig)
    
    grouped_regions_and_seq = extract_grouped_sequences(ref_seq, padding=padding)

    if not grouped_regions_and_seq:
        print("❌ No ambiguous regions found! Exiting.")
        return

    print_ambiguous_regions(grouped_regions_and_seq)

    start, end = get_user_selection(ref_seq)
    print(f"\n🎯 Final UMI Region: ({ref_contig}, {start}, {end})")
    print(f"🔍 Sequence: {ref_seq[start:end]}")
    print("✅ Use these coordinates for UMI extraction!\n")
    
    save_request = input("Would you like to save these coordinates to a file? (y/N): ").strip().lower()
    if save_request == "y":
        save_path = reference_fasta_PATH.parent / f"{reference_fasta_PATH.name}.targetUMIs.csv"
        current_datetime = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        if save_path.exists():
            print(f"\n⚠️ File already exists: {save_path}. Appending to it.")
            with open(save_path, "a") as file:
                file.write(f"{ref_contig},{start},{end},{current_datetime}\n")
        else:
            with open(save_path, "w") as file:
                file.write("ref_contig,start,end,date_time\n")
                file.write(f"{ref_contig},{start},{end},{current_datetime}\n")
        print(f"\n📝 Coordinates saved to {save_path.name} in the reference file directory.\n"
              f"{save_path.absolute()}")
    else:
        print("🚪 Exiting without saving.")


def parse_args():
    parser = argparse.ArgumentParser(description="Select UMI region from reference sequence.")
    parser.add_argument("reference_fasta", type=str,
                        help="Path to reference FASTA file with ambiguous bases contained.")
    parser.add_argument("--padding", type=int, default=5,
                        help="Padding around ambiguous bases for display.")
    parser.add_argument("--contig", type=str, default=None,
                        help="Contig/Chromosome name to extract UMI region from, default is first in FASTA.")
    return parser


def pick_umi_regions(args):
    run_umi_region_picker(args.reference_fasta, contig=args.contig, padding=args.padding)


def main():
    args = parse_args().parse_args()
    pick_umi_regions(args)


if __name__ == "__main__":
    main()
