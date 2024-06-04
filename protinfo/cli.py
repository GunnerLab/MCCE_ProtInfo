#!/usr/bin/env python

__doc__ = """
Command line interface for the MCCE ProtInfo tool, which gathers:
 * Info about the input protein from Bio.PDB parser.
 * Info from MCCE step1 run.log & debug.log when step1 can be run.

This is the 'main' module for the cli, which calls the function
that outputs a single protein report: `get_single_pdb_report(args)`.

Options:
 1. pdb (required): a pdb file name or pdbid (assumed valid).
 2. --fetch (False if not used): If 'pdb' is a pdbid and flag is used,
    the biological assembly is downloaded from rcsb.org.
 Step1 options:
  --wet (False): Keep water molecules.
  --noter (False): Do not label terminal residues (for making ftpl).
  -d (4.0): Protein dielectric constant for delphi.
  -u (''): User selected, comma-separated KEY=var pairs from run.prm; e.g.:
           -u HOME_MCCE=/path/to/mcce_home,EXTRA=./extra.tpl.
  -e (mcce): mcce executable location.
  -h, --help  Show this help message and exit.
  --fetch     Download the biological assembly of given pdb (if not a file).

Usage:
 >ProtInfo 1fat --fetch
 >ProtInfo 1fat.pdb --noter
"""


from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
from IPython.core.formatters import format_display_data
import logging
from pathlib import Path
from protinfo import info, io_utils as iou
import sys
from typing import Union


CLI_NAME = "ProtInfo"
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)


# error msg as fstring:
ERR_FETCH_EXISTING_FILE = """
The input pdb ({}) resolves to an existing file: to OVERWRITE it with
the biological assembly from a fresh download, remove the extension
OR do not pass the --fetch flag at the command line.
"""
ERR_MISSING_FETCH_FLAG = """
The input pdb ({}) seems to be a pdbid. To download its
biological assembly, add --fetch at the command line.
"""


def args_to_str(cliname: str, args: Namespace) -> str:
    """Return cli args to string.
    Note: Using format_display_data to obtain output
    as in a notebookk where the 'func' object ref is
    in readeable form instead of uid.
    """

    return f"{cliname} args:\n{format_display_data(vars(args))[0]['text/plain']}\n"


def validate_pdb_inputs(args: Namespace) -> Union[Path, str, None]:
    """Validate args.pdb and args.fetch"""

    pdb = iou.check_pdb_arg(args.pdb)
    if isinstance(pdb, tuple):
        # error:
        logger.error(pdb[1])
        return None
    elif isinstance(pdb, str):
        if not args.fetch:
            logger.error(ERR_MISSING_FETCH_FLAG.format(pdb))
            return None
    else:
        if args.fetch:
            logger.error(ERR_FETCH_EXISTING_FILE.format(pdb))
            return None

    return pdb


def get_single_pdb_report(args: Union[Namespace, dict]):
    """Get info and save report for a single pdb.
    This function is called by the cli.
    Expected keys in args: pdb [Path, str], fetch [bool].
    Workflow:
      1. validate_pdb_inputs
      2. queries.get_rcsb_pdb(pdb):
         download bio-assembly if pdb is id
      3. collect_info in two dicts:
         first: from Bio.parser, second: from run.log parser
      4. collect_info_lines from the two dicts
      5. save_report
    """

    if isinstance(args, dict):
        args = Namespace(**args)

    pdb = validate_pdb_inputs(args)
    if pdb is None:
        logger.error("Input validation failed.")
        return

    if not isinstance(pdb, Path):
        pdb = iou.get_rcsb_pdb(pdb)
        if not isinstance(pdb, Path):
            logger.error("Could not download from rcsb.org.")
            return

    prot_d, step1_d = info.collect_info(pdb, args)
    report_lines = info.collect_info_lines(prot_d, step1_d)
    iou.save_report(report_lines, pdb_fp=pdb)

    return


def arg_valid_pdb_len(p: str) -> Union[None, str]:
    """Return None if pdb is empty str."""
    if not len(p):
        return None
    return p


def pi_parser():
    p = ArgumentParser(
        prog=f"{CLI_NAME}",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""Report issues here:
        https://github.com/GunnerLab/MCCE_ProtInfo/issues""",
    )
    p.add_argument(
        "pdb",
        type=arg_valid_pdb_len,
        help="""A pdb file name (in the current directory) or
        a pdbid (assumed valid).""",
    )
    p.add_argument(
        "--fetch",
        default=False,
        action="store_true",
        help="Download the biological assembly of given pdb (if not a file).",
    )

    s1 = p.add_argument_group("s1", "step1 options")
    # step1.py prot.pdb {wet}{noter}{d}{u}{e}
    s1.add_argument(
        "--wet",
        default=False,
        action="store_true",
        help="Keep water molecules; %(default)s.",
    )
    s1.add_argument(
        "--noter",
        default=False,
        action="store_true",
        help="Do not label terminal residues (for making ftpl); %(default)s.",
    )
    s1.add_argument(
        "-d",
        metavar="epsilon",
        type=float,
        default=4.0,
        help="protein dielectric constant for delphi; %(default)s.",
    )
    s1.add_argument(
        "-u",
        metavar="Key=Value",
        type=str,
        default="",
        help="""
        User selected, comma-separated KEY=var pairs from run.prm; e.g.:
        -u HOME_MCCE=/path/to/mcce_home,H2O_SASCUTOFF=0.05,EXTRA=./extra.tpl; default: %(default)s.
        Note: No space after a comma!""",
    )
    s1.add_argument(
        "-e",
        metavar="/path/to/mcce",
        default="mcce",
        help="mcce executable location; default: %(default)s.",
    )

    return p


def prot_info_cli(argv=None):
    """Cli 'main' function: produces a Markdown report for a single pdb."""

    cli_parser = pi_parser()
    args = cli_parser.parse_args(argv)
    if args.pdb is None:
        logger.error("No input: you must provide a pdbid or a pdb filename.")
        sys.exit("No input: you must provide a pdbid or a pdb filename.")

    logger.info(args_to_str(CLI_NAME, args))
    get_single_pdb_report(args)

    rpt_fp = Path("ProtInfo.md")
    if rpt_fp.exists():
        with open(rpt_fp) as f:
            print(f.read())

    return


if __name__ == "__main__":
    prot_info_cli(sys.argv)
