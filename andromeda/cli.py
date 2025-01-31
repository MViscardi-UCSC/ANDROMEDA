"""
cli.py
Marcus Viscardi,    January 31, 2025

This allows you to call the different functions from the command line from ~one place.
For example:
```bash
python -m andromeda.cli pick-umi-region -r reference.fasta -c contig_name
```
"""
import argparse
import sys

import andromeda.extract as extract
import andromeda.ref_pos_picker as ref_pos_picker
import andromeda.umi_group as umi_group


def main():
    parser = argparse.ArgumentParser(description="ANDROMEDA CLI: Tools for UMI extraction and processing.")

    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # 📍 Add the UMI Region Picker Command
    umi_picker_parser = subparsers.add_parser("pick-umi-region",
                                              help="Select UMI region from reference sequence")
    umi_picker_subparser = ref_pos_picker.parse_args()
    for action in umi_picker_subparser._actions:
        if action.dest not in {"help"}:
            umi_picker_parser._add_action(action)
    
    # 📍 Add the UMI Extraction Command
    umi_extract_parser = subparsers.add_parser("extract-umis",
                                               help="Extract UMIs from BAM files based on mapped reference positions")
    umi_extract_subparser = extract.parse_args()
    for action in umi_extract_subparser._actions:
        if action.dest not in {"help"}:
            umi_extract_parser._add_action(action)
    
    # 📍 Add the UMI Grouping Command
    umi_group_parser = subparsers.add_parser("group-umis",
                                             help="Group UMIs using `umi_tools group` with a single edit distance.")
    umi_group_subparser = umi_group.parse_args()
    for action in umi_group_subparser._actions:
        if action.dest not in {"help"}:
            umi_group_parser._add_action(action)

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Handle Commands
    if args.command == "pick-umi-region":
        ref_pos_picker.pick_umi_regions(args)
    
    elif args.command == "extract-umis":
        extract.extract_umis_and_summarize(args)
    
    elif args.command == "group-umis":
        umi_group.group_umis(args)


if __name__ == "__main__":
    main()
