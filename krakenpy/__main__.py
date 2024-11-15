"""
This file is part of KrakenPy (https://github.com/rmcolq/krakenpy).
Copyright 2024 Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import argparse
import sys

import krakenpy
import krakenpy.subcommands

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="krakenpy",
        description="Miscellaneous fasta manipulation tools",
    )

    parser.add_argument("--version", action="version", version=krakenpy.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    # _______________________________   common  _________________________________#
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    #common.add_argument(
    #    "--log-file", dest="log_file", metavar='<filename>', required=False, default=None,
    #    help="Log file to use (otherwise uses stdout, or stderr if out-fasta to stdout)"
    #)
    # _______________________________  consensus  __________________________________#

    subparser_merge = subparsers.add_parser(
        "merge",
        parents=[common],
        help="Combines kraken reports and consensus files representing the same dataset "
             "into a single pair of files, with later files in the input list given preference "
             "over earlier files in the list",
    )
    subparser_merge.add_argument(
        '--in-assignments', dest='in_assignments', nargs='+', metavar='<filename>', required=True,
        help='A number of kraken assignment files for the same dataset ordered by preference (later=higher)'
    )
    subparser_merge.add_argument(
        '--in-reports', dest='in_reports', nargs ='+', metavar='<filename>', required=True,
        help='A number of kraken reports for the same dataset ordered by preference (later=higher)'
    )


    subparser_merge.add_argument(
        '--out-prefix', dest='out_prefix', metavar='<filename>', default="merged",
        help='Output prefix for merged result'
    )

    subparser_merge.set_defaults(func=krakenpy.subcommands.merge.run)

    # _______________________________  end  __________________________________#


    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()