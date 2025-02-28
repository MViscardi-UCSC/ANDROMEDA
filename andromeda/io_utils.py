"""
io_utils.py
Marcus Viscardi,    February 26, 2025

Couple of functions for reading and writing files, as well as a function to check if a file exists.
"""
import sys
from pathlib import Path
from typing import Union, List, Dict, Any, Tuple
from Bio import SeqIO
import pandas as pd

from andromeda.logger import log


def load_single_reference_contig(reference_path: str | Path, contig: str = None):
    with open(reference_path, "r") as handle:
        records = [record for record in SeqIO.parse(handle, "fasta")]
    record_ids = [record.id for record in records]
    records_dict = {record.id: record for record in records}
    if len(records) > 1:
        if contig is None:
            log.info(f"Found {len(records)} records in reference file, "
                     f"because you didn't specify a contig we will use "
                     f"the first one: {record_ids[0]}")
            target_record = records[0]
            target_contig = record_ids[0]
        elif contig not in record_ids:
            raise ValueError(f"Contig {contig} not found in reference file, "
                             f"available contigs: {record_ids}")
        else:
            target_record = records_dict[contig]
            target_contig = contig
    else:
        if contig is not None and contig not in record_ids:
            log.info(f"Found 1 record in reference file ({record_ids[0]}), "
                     f"but you specified a contig that isn't present: {contig}.\n"
                     f"We will use the record in the file.")
        elif contig is not None and contig in record_ids:
            log.info(f"Found 1 record in reference file (matching your request for {contig}), using it.")
        target_record = records[0]
        target_contig = record_ids[0]
    return target_contig, str(target_record.seq)


def load_reference_contigs_to_dict(reference_path: str | Path) -> Dict[str, str]:
    with open(reference_path, "r") as handle:
        records = [record for record in SeqIO.parse(handle, "fasta")]
    records_dict = {}
    for record in records:
        records_dict[record.id] = str(record.seq)
        log.trace(f"Loaded contig {record.id} with length {len(record.seq)}")
    log.debug(f"Loaded {len(records)} contigs ({', '.join(records_dict.keys())}) "
              f"from reference file: {reference_path})")
    return records_dict


def load_umi_positions(umi_positions_file: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Loads UMI positions from a TSV file and converts them into a dictionary."""
    contig_df = pd.read_csv(umi_positions_file)
    # Columns = "ref_contig", "start", "end", "date_time"

    output_dict = {}  # Dict of umi_positions per contig
    for contig, row in contig_df.set_index("ref_contig").iterrows():
        if contig not in output_dict:
            # Start a new list of UMI positions for this contig, with an index of 0
            umi_contig_index = 0
            output_dict[contig] = [(umi_contig_index, (row["start"], row["end"]))]
        else:
            # Add the new UMI position to the list for this contig, with an index of the previous position + 1
            prev_index = output_dict[contig][-1][0]
            umi_contig_index = prev_index + 1
            output_dict[contig].append((umi_contig_index, (row["start"], row["end"])))
        log.info(f"📍 UMI Region for {contig} (#{umi_contig_index:0>2}): {row['start']}-{row['end']}")
    return output_dict


if __name__ == '__main__':
    log.remove()
    log.add(sink=sys.stderr, level="TRACE")
    
    test_fasta = Path("../examples/JA-NP-093/references/prok_nanoluc.fa")
    contig, seq = load_single_reference_contig(test_fasta)
    ref_dict = load_reference_contigs_to_dict(test_fasta)