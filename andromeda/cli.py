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
from importlib import import_module
from pathlib import Path

MODULES = [
    "andromeda.ref_pos_picker",
    "andromeda.extract",
    "andromeda.umi_group",
    "andromeda.consensus"
]

def get_dependencies():
    dependencies = {}
    for module_name in MODULES:
        module = import_module(module_name)
        if hasattr(module, "dependencies"):
            dependencies[module_name.split(".")[-1]] = module.dependencies()
    return dependencies

def resolve_dependencies(args, module_name, outputs):
    dependencies = get_dependencies()
    if module_name in dependencies:
        required_inputs = dependencies[module_name]
        for arg, dep in required_inputs.items():
            if not getattr(args, arg, None) and dep in outputs:
                setattr(args, arg, outputs[dep])
    return args

def parse_args():
    parser = argparse.ArgumentParser(description="ANDROMEDA CLI: Modular UMI Extraction and Processing.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    module_parsers = {}
    resolved_deps = get_dependencies()

    for module_name in MODULES:
        module = import_module(module_name)
        command_name = module_name.split(".")[-1]
        subcommand = subparsers.add_parser(command_name, help=f"Run {command_name} module")

        module_parser = module.parse_args()
        module_parsers[command_name] = module_parser
        for action in module_parser._actions:
            if action.dest not in {"help"}:
                subcommand._add_action(action)

    run_all_parser = subparsers.add_parser("run-all",
                                           help="Run all ANDROMEDA modules in sequence.")
    run_all_parser.add_argument("ref_fasta", type=Path,
                                help="Path to reference FASTA file.")
    run_all_parser.add_argument("mapped_bam", type=Path,
                                help="Path to BAM file for UMI extraction.")
    run_all_parser.add_argument("output_parent_dir", type=Path,
                                help="Parent directory for outputs.")
    run_all_parser.add_argument("--umi-positions", type=Path, default=None,
                                help="Path to UMI position TSV file.")

    for module_name, module_parser in module_parsers.items():
        group = run_all_parser.add_argument_group(f"{module_name} Options")
        dependencies = resolved_deps.get(module_name, {})

        for action in module_parser._actions:
            if action.dest not in {"ref_fasta", "output_parent_dir", "mapped_bam", "help"}:
                # Suppress required check if this argument is a dependency
                if action.dest in dependencies:
                    action.required = False
                    action.help = argparse.SUPPRESS
                group._add_action(action)

    return parser.parse_args()


def main():
    args = parse_args()
    outputs = {}

    if args.command in [module.split(".")[-1] for module in MODULES]:
        module = import_module(f"andromeda.{args.command}")
        module.pipeline_main(args)

    elif args.command == "run-all":
        
        assert args.ref_fasta.exists(), f"Reference FASTA file not found: {args.ref_fasta}"
        assert args.mapped_bam.exists(), f"Mapped BAM file not found: {args.mapped_bam}"
        assert args.output_parent_dir.exists(), f"Output parent directory not found: {args.output_parent_dir}"
        assert args.output_parent_dir.is_dir(), f"Output parent directory is not a directory: {args.output_parent_dir}"
        
        if not args.umi_positions:
            # Let's look to see if we can find the UMI position TSV file
            umi_positions = args.ref_fasta.with_suffix(".fasta.targetUMIs.csv")
            if umi_positions.exists():
                use_old_umi_pos = input(f"Found existing UMI position TSV file at {umi_positions}. "
                                        "Use this file? (y/n): ")
                if use_old_umi_pos.lower() == "y":
                    args.umi_positions = umi_positions
                else:
                    args.umi_positions = None
            
        if args.umi_positions:
            assert args.umi_positions.exists(), f"UMI position TSV file not found: {args.umi_positions}"
            MODULES.pop(0)  # Remove ref_pos_picker from the list of modules to run
            # Now we need to artificially add the outputs of ref_pos_picker to the outputs dictionary
            ref_pos_picker_outputs = {
                "ref_fasta": args.ref_fasta,
                "umi_positions": args.umi_positions,
                "output_parent_dir": args.output_parent_dir
            }
            outputs["ref_pos_picker"] = ref_pos_picker_outputs
        for module_name in MODULES:
            command_name = module_name.split(".")[-1]
            module = import_module(module_name)

            module_args = argparse.Namespace(**vars(args))
            module_args = resolve_dependencies(module_args, command_name, outputs)

            module_args.ref_fasta = args.ref_fasta
            module_args.output_parent_dir = args.output_parent_dir
            module_args.mapped_bam = args.mapped_bam

            print(f"\n[DEBUG] Calling {module_name} with arguments:")
            for key, value in vars(module_args).items():
                print(f"  {key}: {value}")

            result = module.pipeline_main(module_args)
            if result:
                outputs[command_name] = result

if __name__ == "__main__":
    main()
