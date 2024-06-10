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
  --dry (False if not used): Keep water molecules.
  --noter (False if not used): Do not label terminal residues (for making ftpl).
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
import logging
from pathlib import Path
from protinfo import cli_opts, info, io_utils as iou
import sys
from typing import Tuple, Union


CLI_NAME = "ProtInfo"
logger = logging.getLogger(__name__)


# error msg as fstring:
ERR_FETCH_EXISTING_FILE = """
The pdb ({}) resolves to an existing file: to OVERWRITE it with the
biological assembly from a fresh download, remove the extension, OR:
do not pass the --fetch flag at the command line to use the existing file.
"""
ERR_MISSING_FETCH_FLAG = """
The input pdb ({}) seems to be a pdbid. To download its
biological assembly, add --fetch at the command line.
"""
ERR_CALL_NOT_IN_FILE_DIR = """
Call ProtInfo from where the pdb resides."""


def check_pdb_arg(input_pdb: str) -> Union[Path, str, Tuple[None, str]]:
    """Validate input_pdb str, which can be either a pdb id or a pdb file.
    The tuple output type is used to pass the error message.
    """

    pdb = Path(input_pdb).resolve()
    if not pdb.exists():
        if not pdb.suffix:
            # if no extension, assume pdbid:
            s = len(pdb.stem)
            if s != 4:
                return None, f"Invalid pdbid length: {s}; 4 expected."

            return input_pdb

        return None, f"File not found: {pdb}"

    if not pdb.parent == Path.cwd():
        return None, ERR_CALL_NOT_IN_FILE_DIR

    if pdb.suffix != ".pdb":
        return None, f"Not a valid extension: {pdb.suffix}"

    return pdb


def validate_pdb_inputs(args: Namespace) -> Union[Path, str]:
    """Validate args.pdb and args.fetch"""

    pdb = check_pdb_arg(args.pdb)
    if isinstance(pdb, tuple):
        # error:
        logger.error(pdb[1])
        sys.exit(pdb[1])

    elif isinstance(pdb, str):
        if not args.fetch:
            logger.error(ERR_MISSING_FETCH_FLAG.format(pdb))
            sys.exit(ERR_MISSING_FETCH_FLAG.format(pdb))
    else:
        if args.fetch:
            logger.error(ERR_FETCH_EXISTING_FILE.format(pdb))
            sys.exit(ERR_FETCH_EXISTING_FILE.format(pdb))

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
    if not isinstance(pdb, Path):
        pdb = iou.get_rcsb_pdb(pdb)
        if not isinstance(pdb, Path):
            logger.error("Could not download from rcsb.org.")
            sys.exit("Could not download from rcsb.org.")

    prot_d, step1_d = info.collect_info(pdb, args)

    report_lines = f"---\n{cli_opts}\n"
    report_lines += info.get_pdb_report_lines(prot_d, step1_d)
    iou.save_report(report_lines, pdb_fp=pdb)

    return


def arg_valid_pdb_len(p: str) -> Union[None, str]:
    """Return None if pdb is empty str."""
    if not len(p):
        return None
    if not p.endswith(".pdb") and len(p) != 4:
        # pdbid given for fetching -> 4 chars
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
    # step1.py prot.pdb {dry}{noter}{d}{u}{e}
    s1.add_argument(
        "-d",
        metavar="epsilon",
        type=float,
        default=4.0,
        help="protein dielectric constant for delphi; %(default)s.",
    )
    s1.add_argument(
        "-e",
        metavar="/path/to/mcce",
        default="mcce",
        help="mcce executable location; default: %(default)s.",
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
        "--dry",
        default=False,
        action="store_true",
        help="Remove water molecules; %(default)s.",
    )
    s1.add_argument(
        "--noter",
        default=False,
        action="store_true",
        help="Do not label terminal residues (for making ftpl); %(default)s.",
    )

    return p


def prot_info_cli(argv=None):
    """Cli 'main' function: produces a Markdown report for a single pdb."""

    cli_parser = pi_parser()
    args = cli_parser.parse_args(argv)
    if args.pdb is None:
        logger.error("No input: you must provide a pdbid or a pdb filename.")
        sys.exit("No input: you must provide a pdbid or a pdb filename.")

    cli_opts.all = vars(args)
    logger.info(cli_opts)

    get_single_pdb_report(args)

    rpt_fp = Path("ProtInfo.md")
    if rpt_fp.exists():
        print(rpt_fp.read_text())

    return


if __name__ == "__main__":
    prot_info_cli(sys.argv)
