"""
io_utils.py
Marcus Viscardi,    February 26, 2025

Couple of functions for reading and writing files, as well as a function to check if a file exists.
"""

from pathlib import Path
from typing import Union, List, Dict, Any
from Bio import SeqIO

from andromeda.logger import log


def load_reference(reference_path: str | Path, contig: str = None):
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
