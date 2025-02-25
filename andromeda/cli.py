"""
cli.py
Marcus Viscardi,    January 31, 2025

This allows you to call the different functions from the command line from ~one place.
For example:
```bash
python andromeda ref_pos_picker ref.fasta output_parent_directory
python andromeda extract ref.fasta mapped.bam output_parent_directory
python andromeda umi_group tagged.bam output_parent_directory
python andromeda consensus ref.fasta grouped.bam output_parent_directory
```
or
```bash
python andromeda run-all ref.fasta mapped.bam output_parent_directory
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

MODULE_NAMES = [module.split(".")[-1] for module in MODULES]

def get_dependencies():
    dependencies = {}
    for module_name in MODULES:
        module = import_module(module_name)
        if hasattr(module, "dependencies"):
            dependencies[module_name.split(".")[-1]] = module.dependencies()
    return dependencies

def resolve_dependencies(args, module_name, outputs):
    # TODO: Make sure this is working as expected!!!!
    dependencies = get_dependencies()
    if module_name in dependencies:
        required_inputs = dependencies[module_name]
        for arg, dep in required_inputs.items():
            already_set = bool(getattr(args, arg, None))
            outputs_source, outputs_name = dep.split(".")
            if outputs_source in outputs and outputs_name in outputs[outputs_source]:
                setattr(args, arg, outputs[outputs_source][outputs_name])
            elif not already_set and outputs_source in outputs:
                setattr(args, arg, outputs[outputs_source])
    return args


def peek_command():
    import sys
    if len(sys.argv) > 1:
        return sys.argv[1]
    return None


def parse_args():
    command_to_run = peek_command()
    
    parser = argparse.ArgumentParser(description="ANDROMEDA CLI: Modular UMI Extraction and Processing.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    module_parsers = {}
    resolved_deps = get_dependencies()

    for module_name in MODULES:
        module = import_module(module_name)
        module_help = module.HELP_TEXT if hasattr(module, "HELP_TEXT") else f"Run {module_name} module."
        command_name = module_name.split(".")[-1]
        subcommand = subparsers.add_parser(command_name,
                                           help=module_help)

        module_parser = module.parse_args()
        module_parsers[command_name] = module_parser
        for action in module_parser._actions:
            if action.dest not in {"help"}:
                subcommand._add_action(action)
    
    if command_to_run in MODULE_NAMES:
        # Don't bother with the run-all parser if we're just running a single module
        # This should help with time a bit too!
        return parser.parse_args()
    
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
                if command_to_run == "run-all" and action.dest in dependencies:
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
            current_suffix = args.ref_fasta.suffix
            umi_positions = args.ref_fasta.with_suffix(current_suffix + ".targetUMIs.csv")
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
            module_actions = [action.dest for action in module.parse_args()._actions]
            for key, value in vars(module_args).items():
                if key in module_actions or key == "command":
                    print(f"  {key}: {value}")
            print(f"\n[DEBUG] Additional arguments carried over for/from other steps:")
            for key, value in vars(module_args).items():
                if key not in module_actions:
                    print(f"  {key}: {value}")

            result = module.pipeline_main(module_args)
            
            print(f"\n[DEBUG] {module_name} returned:")
            print(result)
            
            if result:
                outputs[command_name] = result

if __name__ == "__main__":
    main()
