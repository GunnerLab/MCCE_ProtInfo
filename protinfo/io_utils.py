#!/usr/bin/env python

"""
Module: io_utils
"""

import Bio.PDB as PDB
import gzip
import logging
from pathlib import Path
import requests
import shutil
import subprocess
import sys
from typing import Tuple, Union


logger = logging.getLogger(__name__)
# logger.setLevel(logging.WARNING)


ERR_CALL_NOT_IN_FILE_DIR = """
Call ProtInfo from where the pdb resides."""


def subprocess_run(
    cmd: str,
    capture_output: bool = True,
    check: bool = False,
    text: bool = True,
    shell: bool = True,
) -> Union[subprocess.CompletedProcess, subprocess.CalledProcessError]:
    """Wraps subprocess.run. Return CompletedProcess or err obj."""

    try:
        data = subprocess.run(cmd, capture_output=capture_output, check=check, text=text, shell=shell)
    except subprocess.CalledProcessError as e:
        data = e

    return data


def make_executable(sh_path: str) -> None:
    """Alternative to os.chmod(sh_path, stat.S_IXUSR): permission denied."""

    sh_path = Path(sh_path)
    cmd = f"chmod +x {str(sh_path)}"

    try:
        subprocess_run(cmd, capture_output=False, check=True)
    except subprocess.CalledProcessError:
        logger.exception("Error in subprocess cmd 'chmod +x'")
        raise

    return


def decompress_gz(gfp: Path) -> Path:
    """Decompressed the gzipped file given by its filepath, gfp."""

    fpo = gfp.parent.joinpath(f"{gfp.stem}")
    with gzip.open(gfp, "rb") as f_in:
        with open(fpo, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return fpo


def get_cif_protname(cif_fp: Path):
    """Extract and return the `_struct.title` value from a .cif file."""

    cmd = "grep '^_struct.title' " + str(cif_fp)
    response = subprocess_run(cmd, check=False)

    return response.stdout.split(maxsplit=1)[1]


def get_path_keys(pdb: Path) -> Union[dict, None]:
    """
    Return path keys from run.prm.record as a dict with keys:
    ["topology files","renaming file"].
    """
    d = {}
    desc = ["topology files", "renaming file"]
    runrec = pdb.parent.joinpath("step1_run", "run.prm.record")
    if not runrec.exists():
        logger.error("Not found: run.prm.record")
        return None

    cmd = "grep -E 'MCCE_HOME|RENAME_RULES' " + str(runrec)
    response = subprocess_run(cmd, check=False)
    for i, line in enumerate(response.stdout.splitlines()):
        d[desc[i]] = line[:-1].split("(")[0].strip()

    d[desc[0]] = d[desc[0]] + "/param"

    return d


def cif2pdb(cif_fp: Path) -> Path:
    """Convert a .cif file to a .pdb file.
    The saved pdb file is truncated to the maximum number of atoms
    if their number in the .cif exceeds the 99,999 pdb limit.
    The output file header will be:
        HEADER    <cif_fp.name> converted to pdb by MCCE_ProtInfo; truncated: [True | False]
        TITLE     <protname from .cif file>
    """

    # Warning: If the file contains too many atoms, it will be saved truncated;
    # if auth_residues=False, it will do so silently! Here we need that info for
    # the pdb header
    HDR = "HEADER    {} converted to pdb by MCCE_ProtInfo; truncated: {}\n"
    HDR = HDR + "TITLE     {}\n"

    protname = get_cif_protname(cif_fp)
    cif_out = f"{cif_fp.name[:4]}.pdb"

    parser = PDB.MMCIFParser(auth_residues=True, QUIET=True)
    structure = parser.get_structure("cif", cif_fp)
    N = len(list(structure.get_atoms()))
    too_large = N > 99_999

    io = PDB.PDBIO()
    io.set_structure(structure)  # coords only
    with open(cif_out, "w") as fh:
        try:
            io.save(fh)
        except Exception as e:
            if too_large:
                logger.warning(f"The number of atoms ({N:,}) exceeds the 99,999 PDB format limit: truncated file.")
            else:
                logger.warning(f"Error while saving cif file as {cif_out}: {e}")

    hdr = HDR.format(cif_fp.name, too_large, protname)
    insert_pdb_hdr(cif_out, hdr)

    return cif_out


def insert_pdb_hdr(pdb_from_cif_fp: Path, hdr: str):
    """Insert the pdb header info into a pdb file that was converted from .cif."""

    with open(pdb_from_cif_fp, "r+") as f:
        lines = f.readlines()
        lines.insert(0, hdr)
        f.seek(0)
        f.writelines(lines)

    return


def file_in_current_dir(fp: Path) -> bool:
    """Return whether the given file path is in the current directory."""

    return fp.parent == Path.cwd()


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

    if not file_in_current_dir(pdb):
        return None, ERR_CALL_NOT_IN_FILE_DIR

    if pdb.suffix != ".pdb":
        return None, f"Not a valid extension: {pdb.suffix}"

    return pdb


def save_report(
    report_lines: str,
    pdb_fp: Union[Path, None] = None,
    report_fp: Union[Path, None] = None,
):
    """Write and save the ProtInfo report.
    Args:
      report_lines (str): The lines to write.
      pdb_fp (Union[Path, None], None): The pdb filepath if the
        report is for a single pdb. Cannot be None if report_fp is.
      report_fp (Union[Path, None], None): The filepath of the output
        report if pdb_fp is None.
        Use report_fp if the report_lines are collated from a set
        of runs.
    Note:
      Neither pdb_fp or report_fp can be both None or both set:
      one of them must be None.
    """

    if pdb_fp is None and report_fp is None:
        logger.error("ValueError: pdb_fp and report_fp cannot both be None.")
        sys.exit(1)

    if pdb_fp is not None and report_fp is not None:
        logger.error(
            """
                     ValueError: pdb_fp and report_fp cannot both be valued.
                     Set pdb_fp to None if report_lines come from multiple
                     pdbs, and conversely for a single pdb run."""
        )
        sys.exit(1)

    if pdb_fp is not None:
        # (re)set report_fp
        report_fp = pdb_fp.parent.joinpath("ProtInfo.md")

    with open(report_fp, "w") as fo:
        fo.writelines(report_lines)

    return


def rcsb_download(pdb_fname: str) -> requests.Response:
    url_rscb = "https://files.rcsb.org/download/"
    return requests.get(url_rscb + pdb_fname, allow_redirects=True)


def get_rcsb_pdb(pdbid: str) -> Union[Path, Tuple[None, str]]:
    """Given a pdb id, download the pdb file containing
    the biological assembly from rcsb.org.
    The file is downloaded with a pdb extension.
    """

    pdbid = pdbid.lower()
    bionames = pdbid + ".pdb1", f"{pdbid}-assembly1.cif.gz"
    pdb_file = pdbid + ".pdb"

    content = None
    # list of bool to identify which bio assembly was saved:
    which_ba = [False, False]  # 0: pdb, 1: cif

    # try bio assemblies first:
    r0 = rcsb_download(bionames[0])
    if r0.status_code < 400:
        which_ba[0] = True
        pdb_file = bionames[0][:-1]
        content = r0.content
    else:
        logger.error(f"Error: Could not download the pdb bio assembly:{r0.reason}")

        r1 = rcsb_download(bionames[1])
        if r1.status_code < 400:
            which_ba[1] = True
            pdb_file = bionames[1]
            content = r1.content
        else:
            logger.error(f"Error: Could not download the cif bio assembly:{r1.reason}")

    if which_ba[0] == which_ba[1]:  # both False; last try: legacy pdb
        r2 = rcsb_download(pdb_file)
        if r2.status_code < 400:
            content = r2.content
        else:
            logger.error(f"Error: Could not download the pdb file:{r2.reason}")
            return None, "Error: Could neither download the bio assembly or pdb file."

    # save file:
    with open(pdb_file, "wb") as fo:
        fo.write(content)
    logger.info("Download completed.")

    if which_ba[1]:
        decomp = decompress_gz(Path(pdb_file))
        logger.info(f"{pdb_file} saved & unzipped as {decomp.name}")
        pdb_file = cif2pdb(decomp)

    return Path(pdb_file).resolve()
