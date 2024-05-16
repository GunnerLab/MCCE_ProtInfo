#!/usr/bin/env python

"""
Module: io_utils
"""

import Bio.PDB as PDB
import gzip
import logging
from pathlib import Path
import shutil
import subprocess
from typing import Union


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def subprocess_run(
    cmd: str,
    capture_output: bool = True,
    check: bool = False,
    text: bool = True,
    shell: bool = True,
) -> Union[subprocess.CompletedProcess, subprocess.CalledProcessError]:
    """Wraps subprocess.run. Return CompletedProcess or err obj."""

    try:
        data = subprocess.run(
            cmd, capture_output=capture_output, check=check, text=text, shell=shell
        )
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


def insert_pdb_hdr(pdb_from_cif_fp: Path, hdr: str):
    """Insert the pdb header info into a pdb file that was converted from .cif."""

    with open(pdb_from_cif_fp, "r+") as f:
        lines = f.readlines()
        lines.insert(0, hdr)
        f.seek(0)
        f.writelines(lines)

    return


def cif2pdb(cif_fp: Path) -> Path:
    """Convert a .cif file to a .pdb file.
    The saved pdb file is truncated to the maximum number of atoms
    if their number in the .cif exceeds the 99,999 pdb limit.
    The output file header will be:
        HEADER    <cif_fp.name> converted to pdb by MCCE_ProtInfo; truncated: [True | False]
        TITLE     <protname from .cif file>
    """

    parser = PDB.MMCIFParser(auth_residues=True, QUIET=True)
    # Warning: If the file contains too many atoms, it will be saved truncated;
    # if auth_residues=False, it will do so silently! Here we need that info for
    # the pdb header

    HDR = "HEADER    {} converted to pdb by MCCE_ProtInfo; truncated: {}\n"
    HDR = HDR + "TITLE     {}\n"

    protname = get_cif_protname(cif_fp)
    cif_out = f"{cif_fp.name[:4]}.pdb"

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
                logger.warning(
                    f"The number of atoms ({N:,}) exceeds the 99,999 PDB format limit: truncated file."
                )
            else:
                logger.warning(f"Error while saving cif file as {cif_out}: {e}")

    hdr = HDR.format(cif_fp.name, too_large, protname)
    insert_pdb_hdr(cif_out, hdr)

    return cif_out
