"""
config_builder.py
Marcus Viscardi,    February 19, 2025

This should be a script to build the configuration file for the ANDROMEDA pipeline.
"""

from datetime import datetime
from pathlib import Path
from pprint import pprint
from typing import Dict

from prompt_toolkit import PromptSession
from prompt_toolkit.completion import PathCompleter
from prompt_toolkit.shortcuts import CompleteStyle

from andromeda.ref_pos_picker import parse_args as ref_pos_picker_args
from andromeda.extract import parse_args as extract_args
from andromeda.umi_group import parse_args as umi_group_args
from andromeda.consensus import parse_args as consensus_args


def main():
    # session = PromptSession()
    # path_completer = PathCompleter(expanduser=True)

    arg_options = {}
    for parser in [
        ref_pos_picker_args(),
        extract_args(),
        umi_group_args(),
        consensus_args(),
    ]:
        subparser_actions = {}
        for action in parser._actions:
            if action.dest not in {"help"}:
                subparser_actions[action.dest] = action.help
        arg_options[parser.description] = subparser_actions
        print(parser.prog)
    pprint(arg_options)

    # print("Welcome to the Config Builder CLI!")
    config = {}

    # while True:
    #     try:
    #         key = session.prompt("Enter config key: ")
    #         value = session.prompt("Enter config value: ",
    #                                completer=path_completer,
    #                                complete_style=CompleteStyle.MULTI_COLUMN)
    #         config[key] = value
    #
    #         more = session.prompt("Add more entries? (y/n): ")
    #         if more.lower() != 'y':
    #             break
    #     except KeyboardInterrupt:
    #         break
    #     except EOFError:
    #         break
    #
    # print("\nGenerated Configuration:")
    # for key, value in config.items():
    #     print(f"{key}: {value}")
    #
    # # Save the configuration to a file
    # with open("config.ini", "w") as config_file:
    #     for key, value in config.items():
    #         config_file.write(f"{key}={value}\n")


if __name__ == "__main__":
    main()
