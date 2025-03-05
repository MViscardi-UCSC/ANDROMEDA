"""
test_cli.py
Marcus Viscardi,    March 04, 2025

Tests for the CLI tools in the Andromeda package.
"""

import pytest
import sys
from pathlib import Path
from andromeda.cli import parse_args, main, run_all_pipeline


def test_cli_parse_args_run_all(monkeypatch):
    # Simulate: python andromeda run-all ./ref.fa ./reads.bam ./out_dir
    test_args = [
        "andromeda",
        "run-all",
        "./some_ref.fasta",
        "./some_mapped.bam",
        "./some_output_dir",
    ]
    monkeypatch.setattr(sys, "argv", test_args)
    args = parse_args()

    assert args.command == "run-all"
    assert Path(args.ref_fasta) == Path("./some_ref.fasta")
    assert Path(args.mapped_bam) == Path("./some_mapped.bam")
    assert Path(args.output_parent_dir) == Path("./some_output_dir")


def test_cli_main_no_command(monkeypatch, capfd):
    # Simulate: python andromeda
    test_args = ["andromeda"]
    monkeypatch.setattr(sys, "argv", test_args)

    # Run main() and expect help output
    with pytest.raises(SystemExit):
        main()

    captured = capfd.readouterr()
    assert "usage: andromeda" in captured.err


def test_cli_missing_ref_fasta(monkeypatch):
    # Testing this b/c it will trip before a missing BAM or output parent dir!
    test_args = [
        "andromeda",
        "run-all",
        "./non_existent_ref.fasta",
        "./some_mapped.bam",
        "./some_output_dir",
    ]
    monkeypatch.setattr(sys, "argv", test_args)
    args = parse_args()
    with pytest.raises(AssertionError) as exc_info:
        run_all_pipeline(args)
    assert "Reference FASTA file not found" in str(exc_info.value)
