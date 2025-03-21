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
import sys, os
from importlib import import_module
from importlib.resources import Package
from pathlib import Path
import tomllib
import textwrap
import importlib.metadata as import_metadata
import shutil

from andromeda.logger import log
from andromeda.big_text import big_text, big_text_side_char, big_text_closer

MODULES = [
    "andromeda.ref_pos_picker",
    "andromeda.extract",
    "andromeda.umi_group",
    "andromeda.consensus",
]

MODULE_NAMES = [module.split(".")[-1] for module in MODULES]

PACKAGE_DIR = Path(__file__).resolve().parent


def get_project_metadata():
    """Retrieve package metadata dynamically, whether installed or run locally."""
    package_name = PACKAGE_DIR.name

    try:
        # Try to get metadata when installed
        metadata = import_metadata.metadata(package_name)
        return {
            "name": metadata["Name"],
            "version": metadata["Version"],
            "description": metadata["Summary"],
            "authors": metadata["Author"].split(", "),
            "license": metadata["License"],
            "requires-python": metadata["Requires-Python"],
            "source": metadata["Project-URL"].split(", ")[1],
        }
    except import_metadata.PackageNotFoundError:
        pass  # Fall back to searching for pyproject.toml

    # Search for pyproject.toml in parent directories
    current_path = Path(__file__).resolve()
    while current_path != current_path.root:
        pyproject_path = current_path / "pyproject.toml"
        if pyproject_path.exists():
            with pyproject_path.open("rb") as f:
                pyproject_data = tomllib.load(f)
                project_data = pyproject_data.get("project", {})
                return {
                    "name": project_data.get("name"),
                    "version": project_data.get("version"),
                    "description": project_data.get("description"),
                    "authors": project_data.get("authors", []),
                    "license": project_data.get("license"),
                    "requires-python": project_data.get("requires-python"),
                    "source": project_data.get("urls", {}).get("source"),
                }
        current_path = current_path.parent

    raise FileNotFoundError("Could not find package metadata or pyproject.toml!")


def is_uv_run():
    return os.getenv("UV", None) is not None


def get_terminal_width():
    """Returns the width of the terminal or a default if unavailable."""
    return shutil.get_terminal_size().columns


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


def get_version():
    return get_project_metadata()['version']


def global_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--version",
        action="version",
        version=f"ANDROMEDA v{get_version()}",
        help="Show the version number and exit.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
        help="Set the log level for the stdout logger [default: INFO].",
    )
    parser.add_argument(
        "--log-file-level",
        default="DEBUG",
        choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
        help="Set the log level for the file logger [default: DEBUG].",
    )
    return parser


def parse_args():
    # Parse Global Arguments First
    global_args = global_parser()
    global_only_parser = argparse.ArgumentParser(parents=[global_args], add_help=False)
    global_only_args, remaining_argv = global_only_parser.parse_known_args()

    # Main Parser
    parser = argparse.ArgumentParser(
        description="ANDROMEDA CLI: Modular UMI Map-Based Extraction, Grouping, and Collapsing."
    )
    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Available commands"
    )

    module_parsers = {}
    resolved_deps = get_dependencies()

    # Add Global Arguments to Each Module's Subparser
    for module_name in MODULES:
        module = import_module(module_name)
        command_name = module_name.split(".")[-1]
        help_text = (
            module.HELP_TEXT if module.HELP_TEXT else f"Run {command_name} module"
        )
        subcommand = subparsers.add_parser(command_name, help=help_text)

        # Add global arguments explicitly for help text visibility
        global_parser().print_help = subcommand.print_help
        for action in global_parser()._actions:
            subcommand._add_action(action)

        module_parser = module.parse_args()
        module_parsers[command_name] = module_parser
        for action in module_parser._actions:
            if action.dest not in {"help"}:
                subcommand._add_action(action)

    # "Run-All" Subparser
    run_all_parser = subparsers.add_parser(
        "run-all", help="Run all ANDROMEDA modules in sequence."
    )
    run_all_parser.add_argument(
        "ref_fasta", type=Path, help="Path to reference FASTA file."
    )
    run_all_parser.add_argument(
        "mapped_bam", type=Path, help="Path to BAM file for UMI extraction."
    )
    run_all_parser.add_argument(
        "output_parent_dir", type=Path, help="Parent directory for outputs."
    )
    run_all_parser.add_argument(
        "--umi-positions",
        type=Path,
        default=None,
        help="Path to UMI position TSV file.",
    )

    # Add global arguments explicitly to run-all
    for action in global_parser()._actions:
        run_all_parser._add_action(action)

    # Add module-specific options
    for module_name, module_parser in module_parsers.items():
        group = run_all_parser.add_argument_group(f"{module_name} Options")
        dependencies = resolved_deps.get(module_name, {})

        for action in module_parser._actions:
            if action.dest not in {
                "ref_fasta",
                "output_parent_dir",
                "mapped_bam",
                "help",
            }:
                if peek_command() == "run-all" and action.dest in dependencies:
                    action.required = False
                    action.help = argparse.SUPPRESS
                group._add_action(action)

    # Parse Subcommand Arguments
    args = parser.parse_args(remaining_argv)

    # Attach Global Arguments to Parsed Args
    args.log_level = global_only_args.log_level
    return args


def run_all_pipeline(args):
    outputs = {}
    assert args.ref_fasta.exists(), f"Reference FASTA file not found: {args.ref_fasta}"
    assert args.mapped_bam.exists(), f"Mapped BAM file not found: {args.mapped_bam}"
    if not args.output_parent_dir.exists():
        assert args.output_parent_dir.parent.exists(), (
            f"Parent directory not found: {args.output_parent_dir.parent}"
        )
        log.info(f"Creating output directory: {args.output_parent_dir}")
        args.output_parent_dir.mkdir(parents=True, exist_ok=False)
    assert args.output_parent_dir.is_dir(), (
        f"Output parent directory is not a directory: {args.output_parent_dir}"
    )

    log.success("Required files found, starting pipeline.")

    if not args.umi_positions:
        # Let's look to see if we can find the UMI position TSV file
        current_suffix = args.ref_fasta.suffix
        umi_positions = args.ref_fasta.with_suffix(current_suffix + ".targetUMIs.csv")
        # First, let's look if the umi_file exists, then we can look in the output to see if we made it previously!
        if not umi_positions.exists():
            log.debug(f"Could not find UMI position TSV file with the reference path at {umi_positions}")
            log.debug(f"Checking output directory for existing UMI position TSV file from previous pipeline runs.")
            alt_positions_file = args.output_parent_dir / "references" / umi_positions.name
            if alt_positions_file.exists():
                umi_positions = alt_positions_file
        
        if umi_positions.exists():
            if not args.extraction_do_not_confirm:
                use_old_umi_pos = input(
                    f"Found existing UMI position TSV file at {umi_positions}. "
                    "Use this file? (y/n): "
                )
                if use_old_umi_pos.lower() == "y":
                    args.umi_positions = umi_positions
                    log.debug(
                        f"Using existing UMI position TSV file at {umi_positions}. The user confirmed this."
                    )
                else:
                    args.umi_positions = None
            else:
                args.umi_positions = umi_positions
                log.debug(f"Using existing UMI position TSV file at {umi_positions}.")
    else:
        assert args.umi_positions.exists(), (
            f"UMI position TSV file not found: {args.umi_positions}"
        )
        log.debug(f"Using provided UMI position TSV file at {args.umi_positions}.")

    if args.umi_positions:
        assert args.umi_positions.exists(), (
            f"UMI position TSV file not found: {args.umi_positions}"
        )
        MODULES.pop(0)  # Remove ref_pos_picker from the list of modules to run
        # Now we need to artificially add the outputs of ref_pos_picker to the outputs dictionary
        log.debug("Skipping ref_pos_picker, using provided UMI position TSV file.")
        ref_pos_picker_outputs = {
            "ref_fasta": args.ref_fasta,
            "umi_positions": args.umi_positions,
            "output_parent_dir": args.output_parent_dir,
        }
        outputs["ref_pos_picker"] = ref_pos_picker_outputs

    for module_name in MODULES:
        command_name = module_name.split(".")[-1]
        log.success(f"Starting to run {module_name}!")

        module = import_module(module_name)
        log.trace(f"Imported module: {module_name}")

        module_args = argparse.Namespace(**vars(args))
        module_args = resolve_dependencies(module_args, command_name, outputs)

        module_args.ref_fasta = args.ref_fasta
        module_args.output_parent_dir = args.output_parent_dir
        module_args.mapped_bam = args.mapped_bam

        actions_debug_str = (
            f"Argument debugging message below:"
            f"\n  Calling {module_name} with arguments:"
        )
        module_actions = [action.dest for action in module.parse_args()._actions]
        for key, value in vars(module_args).items():
            if key in module_actions or key == "command":
                actions_debug_str += f"\n    {key}: {value}"
        actions_debug_str += (
            "\n  Additional arguments carried over for/from other steps:"
        )
        for key, value in vars(module_args).items():
            if key not in module_actions:
                actions_debug_str += f"\n    {key}: {value}"
        log.debug(actions_debug_str)

        result = module.pipeline_main(module_args)

        log.success(f"{module_name} completed successfully!")

        if result:
            log.debug(f"{module_name} returned: {result}")
            outputs[command_name] = result


def print_andromeda_header(spacer_from_left=5):
    data = get_project_metadata()
    project_items = [
        f"Project: {data['name']}",
        f"Version: v{data['version']}",
        f"Description: {data['description']}",
        f"Python Required: {data['requires-python']}",
    ]
    try:
        authors_str = f"Authors: {', '.join(data['authors'])}"
        project_items.append(authors_str)
    except TypeError:
        pass
    project_items.append(f"Source: {data['source']}")
    if is_uv_run():
        project_items.append("")
        project_items.append("Running using uv (look out we got a superuser over here!)")

    wrap_width = len(big_text_closer) - 2 - spacer_from_left
    
    if get_terminal_width() < len(big_text_closer):
        # We don't have enough space to print the big text, so we'll just print the project items
        print("\nANDROMEDA")
        print("=========")
        for line in project_items:
            print(line)
        print("\n\n")
    else:
        print(big_text)
        for line in project_items:
            for wrapped_line in textwrap.wrap(line, width=wrap_width) or [""]:
                table_print_str = (
                    big_text_side_char
                    + " " * spacer_from_left
                    + f"{wrapped_line:<{wrap_width}}"
                    + big_text_side_char
                )
                print(table_print_str)
        print(big_text_closer)
    return None


@log.catch
def main():
    print_andromeda_header()
    args = parse_args()

    log.remove()
    log.add(
        sys.stderr,
        level=args.log_level,
        format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> "
        "| <level>{level: <8}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan> - "
        "<level>{message}</level>",
    )

    if args.command in [module.split(".")[-1] for module in MODULES]:
        log.add(
            f"{args.output_parent_dir}/andromeda_{args.command}" + "_{time}.log",
            level=args.log_file_level,
        )
        log.success(f"Starting ANDROMEDA {args.command} pipeline!")
        log.success(f"CLI called: {' '.join(sys.argv)}")
        log.success(f"Log level: {args.log_level}")
        log.success(f"Log file level: {args.log_file_level}")
        log.success(f"Parsed Arguments: {args}")
        module = import_module(f"andromeda.{args.command}")
        module.pipeline_main(args)
    elif args.command == "run-all":
        log.add(
            f"{args.output_parent_dir}/andromeda_run-all" + "_{time}.log",
            level=args.log_file_level,
        )
        log.success(f"Starting full run of ANDROMEDA pipeline!")
        log.info(f"CLI called: {' '.join(sys.argv)}")
        log.info(f"Log level: {args.log_level}")
        log.info(f"Log file level: {args.log_file_level}")
        log.debug(f"Parsed Arguments: {args}")
        run_all_pipeline(args)
    else:
        # I want to trigger the more extensive help call (that comes from `run-all --help`) here
        # This should help with the case where someone doesn't provide a command
        log.warning("No command provided, please provide a command to run.")
        parse_args().print_help()


if __name__ == "__main__":
    main()
