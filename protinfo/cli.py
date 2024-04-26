#!/usr/bin/env python

__doc__ = """
Command line interface for the MCCE ProtInfo tool, which gathers:
 * Info about the input protein from Bio.PDB parser
 * Info from MCCE step1 run.log & debug.log when step1 can be run

Options:
 -input_pdb (required): a pdb file name or pdbid (assumed valid)
 --fetch (False if not used): If input_pdb is a pdbid and flag is used,
   the biological assembly is downloaded from rcsb.org.
"""


from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
from protinfo import info
import sys


CLI_NAME = "ProtInfo"

def pi_parser():

    p = ArgumentParser(
            prog = f"{CLI_NAME}",
            description = __doc__,
            formatter_class = RawDescriptionHelpFormatter,
        )
    p.add_argument(
        "-input_pdb",
        required = True,
        type = str,
        help = "A pdb file name (in the current directory) or a pdbid (assumed valid)."
    )
    p.add_argument(
        "--fetch",
        default = False,
        action = "store_true",
        help = "Download the biological assembly for input_pdb (if not a file)."
    )

    return p


def prot_info_cli(argv=None):
    """Cli 'main' function: produces a single pdb report."""

    cli_parser = pi_parser()
    args = cli_parser.parse_args(argv)
    info.get_single_pdb_report(args)

    return


if __name__ == "__main__":

    prot_info_cli(sys.argv)
