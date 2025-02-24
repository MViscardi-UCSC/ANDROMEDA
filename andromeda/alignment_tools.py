"""
alignment_tools.py
Marcus Viscardi,    January 31, 2025


This is all code from my old Jupyter Notebook that I used to extract UMIs from BAM files.
It can be found here:
    /data16/marcus/working/240118_COMETm_TadUMICollapsing/COMETm/mapBased/mapBased_umiExtraction.ipynb
And from my codebase in the `cs_parsing.py` file:
    /data16/marcus/working/240118_COMETm_TadUMICollapsing/COMETm/mapBased/cs_parsing.py
"""
import pysam
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict
from tqdm.auto import tqdm
import re

from andromeda.phred_tools import NucleotideQuality

IUPAC_DNA = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT"
}


def bam_to_tagged_bam(bam_file_path: Path,
                      target_chr: str,
                      ref_seq: str,
                      umi_ref_start: int,
                      umi_ref_end: int,
                      flanking_seq_to_capture=0,
                      save_dir=None,
                      mkdir=True,
                      read_id_suffix="",
                      add_read_id_suffix_from_tag=None,
                      save_suffix=".tagged.bam",
                      max_del_in_umi=0,
                      max_ins_in_umi=0,
                      max_iupac_mismatches=0,
                      restrict_to_length=True,
                      bam_tag_dict: Dict[str, str] = None,
                      subset_count=-1) -> Path:
    """
    Extract UMIs from a BAM file and write them to a new BAM file with the UMI sequence
    and deletion count stored in tags.

    :param bam_file_path: Path to the input BAM file.
    :param target_chr: Target chromosome/contig name.
    :param ref_seq: Reference sequence as a string.
    :param umi_ref_start: Start position of the UMI region in the reference sequence.
    :param umi_ref_end: End position of the UMI region in the reference sequence.
    :param flanking_seq_to_capture: Number of flanking bases to capture around the UMI region.
    :param save_dir: Directory to save the output BAM file.
    :param mkdir: Whether to create the save directory if it does not exist.
    :param read_id_suffix: Suffix to add to read IDs in the output BAM file.
    :param add_read_id_suffix_from_tag: Tag from which to add a suffix to read IDs.
    :param save_suffix: Suffix for the output BAM file name.
    :param max_del_in_umi: Maximum allowed deletions in the UMI sequence.
    :param max_ins_in_umi: Maximum allowed insertions in the UMI sequence.
    :param max_iupac_mismatches: Maximum allowed IUPAC mismatches in the UMI sequence.
    :param restrict_to_length: Whether to restrict UMIs to a specific length.
    :param bam_tag_dict: Dictionary mapping UMI attributes to BAM tags.
    :param subset_count: Number of reads to process (use -1 for all reads).
    :return: Path to the output BAM file.
    """
    assert bam_file_path.exists(), f"The provided bam file {bam_file_path} does not exist"

    adj_umi_ref_start = umi_ref_start - flanking_seq_to_capture
    adj_umi_ref_end = umi_ref_end + flanking_seq_to_capture

    bam_tag_dict_required_keys = [
        "umi_sequence",
        "deletion_count",
        "insertion_count",
        "mismatch_count",
        "umi_length"]
    if bam_tag_dict is None:
        bam_tag_dict = {
            "umi_sequence": "uM",
            "deletion_count": "uD",
            "insertion_count": "uI",
            "mismatch_count": "um",
            "umi_length": "uL",
        }
    else:
        for key in bam_tag_dict_required_keys:
            if key not in bam_tag_dict:
                raise ValueError(f"The provided bam_tag_dict is missing the required key: {key}")

    with pysam.AlignmentFile(bam_file_path, 'rb') as input_bam:
        iterator_total = input_bam.count(target_chr, umi_ref_start, umi_ref_end)
        if 0 < subset_count < iterator_total:
            iterator_total = subset_count
        bam_iterator = tqdm(enumerate(input_bam.fetch(target_chr,
                                                      umi_ref_start,
                                                      umi_ref_end)),
                            total=iterator_total,
                            desc="Extracting UMIs")
        if save_dir is None:
            output_bam_path = bam_file_path.with_suffix(save_suffix)
        else:
            if not Path(save_dir).exists() and not mkdir:
                raise FileNotFoundError(f"The save directory {save_dir} does not exist and mkdir is set to False")
            elif not Path(save_dir).exists() and mkdir:
                Path(save_dir).mkdir(exist_ok=True)
            output_bam_path = Path(save_dir) / bam_file_path.with_suffix(save_suffix).name
        failed_bam_path = output_bam_path.with_suffix(".failed.bam")
        umi_success_count = 0
        umi_drop_count = 0
        dropped_for_del, dropped_for_ins, dropped_for_mismatch = 0, 0, 0
        dropped_for_length, dropped_for_too_long, dropped_for_too_short = 0, 0, 0
        umi_output_set = set()
        with pysam.AlignmentFile(output_bam_path, "wb", header=input_bam.header) as out_bam, \
                pysam.AlignmentFile(failed_bam_path, "wb", header=input_bam.header) as failed_bam:
            for i, entry in bam_iterator:
                try:
                    entry_dict = extract_ref_and_query_region(entry, ref_seq,
                                                              adj_umi_ref_start,
                                                              adj_umi_ref_end)
                    extracted_seq: str = entry_dict["query_sequence"]
                    ins_count: int = entry_dict["ins_count"]
                    del_count: int = entry_dict["del_count"]
                    mismatch_count: int = entry_dict["mismatch_count"]
                    if restrict_to_length:
                        if len(extracted_seq) != adj_umi_ref_end - adj_umi_ref_start + 1:
                            length_cutoff_passed = False
                        else:
                            length_cutoff_passed = True
                    else:
                        length_cutoff_passed = True

                    # Add all the new tags to the entry
                    entry.set_tag(bam_tag_dict['umi_sequence'], extracted_seq, value_type='Z')
                    entry.set_tag(bam_tag_dict['deletion_count'], del_count, value_type='i')
                    entry.set_tag(bam_tag_dict['insertion_count'], ins_count, value_type='i')
                    entry.set_tag(bam_tag_dict['mismatch_count'], mismatch_count, value_type='i')
                    entry.set_tag(bam_tag_dict['umi_length'], len(extracted_seq), value_type='i')
                    entry.query_name += read_id_suffix
                    if del_count <= max_del_in_umi \
                            and ins_count <= max_ins_in_umi \
                            and mismatch_count <= max_iupac_mismatches \
                            and length_cutoff_passed:
                        umi_output_set.add(extracted_seq)
                        if add_read_id_suffix_from_tag:
                            entry.query_name += f"_{entry.get_tag(add_read_id_suffix_from_tag)}"
                        out_bam.write(entry)
                        umi_success_count += 1
                    else:
                        umi_drop_count += 1
                        if del_count > max_del_in_umi:
                            dropped_for_del += 1
                        if ins_count > max_ins_in_umi:
                            dropped_for_ins += 1
                        if mismatch_count > max_iupac_mismatches:
                            dropped_for_mismatch += 1
                        if not length_cutoff_passed:
                            dropped_for_length += 1
                            if len(extracted_seq) > adj_umi_ref_end - adj_umi_ref_start + 1:
                                dropped_for_too_long += 1
                            if len(extracted_seq) < adj_umi_ref_end - adj_umi_ref_start + 1:
                                dropped_for_too_short += 1
                        failed_bam.write(entry)

                    if i % 100 == 0:
                        bam_iterator.desc = (f"Extracting UMIs | {i:,} reads | "
                                             f"{len(umi_output_set):,} fitting UMIs | "
                                             f"{umi_drop_count:,} dropped UMIs")
                    if 0 < subset_count <= i:
                        break
                except IndexError as e:
                    print(entry.reference_start, entry.reference_end, entry.get_tag('cs'), e)
                    continue
                except ValueError as e:
                    print(entry.reference_start, entry.reference_end, entry.get_tag('cs'), e)
                    continue
                except Exception as e:
                    print(entry.reference_start, entry.reference_end, entry.get_tag('cs'), e)
                    raise e
    summary_string = (
        f"\n"
        f"Extracted {len(umi_output_set):>8,} unique UMIs from {i + 1:>8,} reads\n"
        f"Wrote     {umi_success_count:>8,} reads to {output_bam_path.name}\n"
        f"Dropped   {umi_drop_count:>8,} UMIs for having too many deletions, insertions, or mismatches\n"
        f"Breakdown (reads can fit into multiple categories): \n"
        f"\tDropped {dropped_for_del:>8,} UMIs for having too many deletions (>{max_del_in_umi})\n"
        f"\tDropped {dropped_for_ins:>8,} UMIs for having too many insertions (>{max_ins_in_umi})\n"
        f"\tDropped {dropped_for_mismatch:>8,} UMIs for having too many mismatches (>{max_iupac_mismatches})\n")
    if restrict_to_length:
        summary_string += (f"\tDropped {dropped_for_length:>8,} UMIs for not being the expected length "
                           f"({adj_umi_ref_end - adj_umi_ref_start + 1} nts)\n"
                           f"\t\tDropped {dropped_for_too_long:>8,} UMIs for being too long\n"
                           f"\t\tDropped {dropped_for_too_short:>8,} UMIs for being too short\n")
    print(summary_string)
    with open(output_bam_path.with_suffix(".summary.txt"), "w") as summary_file:
        summary_file.write(summary_string)
    return output_bam_path


def extract_ref_and_query_region(target_entry: pysam.AlignedSegment, ref_seq: str,
                                 region_start: int, region_end: int) -> dict:
    """
    Extracts reference and query regions from a given aligned segment.

    This function extracts the reference and query positions and sequences from a given aligned segment that
    fall within a specified range. It also provides the option to write the output to a file and print the output.
    
    The final outputs for a read with the following cigar string would look like this:
        2M2D3M2I3M
        
        query_positions: [0, 1, 2, None, None, 3, 4, 5, 6, 7, 8, 9, 10]
        
        query_sequence:   A  T  G     .     .  A  T  G  a  t  A  T   G
        
        ref_positions:   [0, 1, 2, 3, 4, 5, 6, 7, None, None, 8, 9, 10]
        
        ref_sequence:     A  T  G  t  g  A  T  G     .     .  A  T   G
    
    

    Args:
        target_entry (pysam.AlignedSegment): The aligned segment from which to extract the reference and query regions.
        ref_seq (str): The reference sequence.
        region_start (int): The start of the range within which to extract the reference and query regions.
        region_end (int): The end of the range within which to extract the reference and query regions.

    Returns:
        dict:
            ref_positions: The positions along the reference that the query aligned to (with deletions marked)
            ref_sequence: The sequence of the reference that the query aligned to (with deletions marked)
            query_positions: The positions along the query that aligned to the reference (with insertions marked)
            query_sequence: The sequence of the query that aligned to the reference (with insertions marked)
            ins_count: The number of insertions in the query sequence
            del_count: The number of deletions in the query sequence
            mismatch_count: The number of mismatches between the query and reference sequences
            perfect_match: Whether the query sequence perfectly matches the reference sequence
    """
    real_ref_seq = ref_seq
    region_ref_positions = extract_positions_within_range(target_entry.get_reference_positions(full_length=True),
                                                          region_start, region_end)
    region_ref_sequence = [real_ref_seq[i] if i else "." for i in region_ref_positions]
    # Now we need to be able to pull out the same nucleotides for the actual query sequence
    aligned_pairs = target_entry.get_aligned_pairs(with_seq=False)
    region_query_positions = extract_query_positions_from_ref(aligned_pairs, region_start, region_end)
    region_query_sequence = [target_entry.query_sequence[i] if i else "." for i in region_query_positions]
    
    region_query_phreds = [NucleotideQuality(q_score=target_entry.query_qualities[i]).to_phred_char()
                           if i else " " for i in region_query_positions]
    
    region_ref_sequence_matched, region_query_sequence_matched = [], []
    ins_count, del_count, was_perfect, mismatch_count = 0, 0, True, 0
    if not region_ref_sequence or not region_query_sequence:
        was_perfect = False
    for ref, query, phred in zip(region_ref_sequence, region_query_sequence, region_query_phreds):
        ref, query = ref.upper(), query.upper()
        if ref == query:
            region_ref_sequence_matched.append(ref)
            region_query_sequence_matched.append(query)
        elif ref == ".":
            region_ref_sequence_matched.append(ref.lower())
            region_query_sequence_matched.append(query.lower())
            ins_count += 1
            was_perfect = False
        elif query == ".":
            region_ref_sequence_matched.append(ref.lower())
            region_query_sequence_matched.append(query.lower())
            del_count += 1
            was_perfect = False
        elif query in IUPAC_DNA[ref]:
            region_ref_sequence_matched.append(ref)
            region_query_sequence_matched.append(query)
        else:
            region_ref_sequence_matched.append(ref.lower())
            region_query_sequence_matched.append(query.lower())
            was_perfect = False
            mismatch_count += 1
    return_dict = {'ref_positions': region_ref_positions,
                   'ref_sequence': ''.join(region_ref_sequence_matched),
                   'query_positions': region_query_positions,
                   'query_sequence': ''.join(region_query_sequence_matched),
                   'ins_count': ins_count, 'del_count': del_count, 'mismatch_count': mismatch_count,
                   'perfect_match': was_perfect}
    return return_dict


def extract_ref_oriented_sequences(reads: List[pysam.AlignedSegment], reference_seq: str) -> List[str]:
    pass


def extract_positions_within_range(positions, start, end):
    """
    Extract positions within a given range from a list of positions.
    None values, representing deletions, are also included if they fall within the range.

    Args:
        positions (list): List of positions.
        start (int): Start of the range.
        end (int): End of the range.

    Returns:
        list: Positions within the given range.
    """
    result = []
    for pos in positions:
        if pos is None:
            if result:  # if result list is not empty, append None
                result.append(pos)
        elif start <= pos <= end:
            result.append(pos)
        elif pos > end:
            break  # exit loop if position is beyond the end of the range
    return result


def extract_query_positions_from_ref(aligned_pairs, start, end):
    """
    Extract the query positions from a list of aligned pairs that fall within a given range.

    Args:
        aligned_pairs (list): List of aligned pairs.
        start (int): Start of the range.
        end (int): End of the range.

    Returns:
        list: Query positions within the given range.
    """
    result = []
    for query_pos, ref_pos in aligned_pairs:
        if ref_pos is None:
            if result:  # if result list is not empty, append None
                result.append(query_pos)
        elif start <= ref_pos <= end:
            result.append(query_pos)
        elif ref_pos > end:
            break  # exit loop if position is beyond the end of the range
    return result


def bam_to_df(bam_file_path: Path, subset_count=-1) -> pd.DataFrame:
    assert bam_file_path.exists(), f"The provided bam file {bam_file_path} does not exist"
    with pysam.AlignmentFile(bam_file_path, 'rb') as input_bam:
        iterator_total = input_bam.count()
        if 0 < subset_count < iterator_total:
            iterator_total = subset_count
        bam_iterator = tqdm(enumerate(input_bam.fetch()), total=iterator_total, desc="Extracting Reads")
        output_dict = {}
        for i, entry in bam_iterator:
            output_dict[i] = entry.to_dict()
            output_dict[i].pop("tags")
            output_dict[i].update(entry.get_tags())
            if 0 < subset_count <= i:
                break
    return pd.DataFrame(output_dict).T


def regex_parse_long_cs(cs_string):
    """
    This is going to parse the long format cs tags into something that I can try to use to walk the reference and query
    :param cs_string: 
    :return: Tuple containing the following:
        ref_seq_full: 
        query_seq_full:
        ref_seq_ref_oriented: 
        query_seq_ref_oriented: 
        collapsed_hits: 
    """
    re_pattern = (r"(?P<match>\=[ACGTN]+)|"
                  r"(?P<sub>\*[acgtn][acgtn])|"
                  r"(?P<insert>\+[acgtn]+)|"
                  r"(?P<del>-[acgtn]+)|"
                  r"(?P<intron>\~[acgtn]{2}[0-9]+[acgtn]{2})")

    collapsed_hits = []

    ref_seq_full = ""
    ref_seq_ref_oriented = ""
    query_seq_full = ""
    query_seq_ref_oriented = ""

    for re_hits in re.findall(re_pattern, cs_string):
        match, sub, insert, deletion, intron = re_hits
        assert sum(map(bool, re_hits)) == 1
        if match:
            collapsed_hits.append(match)
            match_seq = match[1:]

            ref_seq_full += match_seq
            query_seq_full += match_seq

            ref_seq_ref_oriented += match_seq
            query_seq_ref_oriented += match_seq
        elif sub:
            collapsed_hits.append(sub)
            sub_seq_ref, sub_seq_query = sub[1], sub[2]

            ref_seq_full += sub_seq_ref
            query_seq_full += sub_seq_query

            ref_seq_ref_oriented += sub_seq_ref
            query_seq_ref_oriented += sub_seq_query
        elif insert:
            collapsed_hits.append(insert)
            insert_seq = insert[1:]

            ref_seq_full += "i" * len(insert_seq)
            query_seq_full += insert_seq

            # ref_seq_ref_oriented is unchanged
            # query_seq_ref_oriented is unchanged
        elif deletion:
            collapsed_hits.append(deletion)
            del_seq = deletion[1:]

            ref_seq_full += del_seq
            query_seq_full += "d" * len(del_seq)

            ref_seq_ref_oriented += del_seq
            query_seq_ref_oriented += "d" * len(del_seq)
        elif intron:
            collapsed_hits.append(intron)
            raise NotImplementedError("I haven't implemented the intron parsing yet...")
        else:
            raise ValueError(f"Something went wrong with the regex parsing of the cs tag: {cs_string}")
    return ref_seq_full, query_seq_full, ref_seq_ref_oriented, query_seq_ref_oriented, collapsed_hits

