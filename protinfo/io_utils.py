#!/usr/bin/env python

"""
Module: io_utils
"""

import logging
from pathlib import Path
from protinfo import queries as qrs
import subprocess
from typing import Union


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
#......................................................................

def maybe_download(pdbid:str, fetch:bool) -> Path:

    if not fetch:
        msg = f"""The input_pdb parameter ({pdb}) seems to be a pdbid. To download
        its biological assembly, the fetch option must be True."""
        raise TypeError(msg)

    pdb = qrs.get_rcsb_pd(pdb)
    if pdb is None:
        raise ValueError(f"Could not download {pdb} from rcsb.org.")
   
    return pdb


def subprocess_run(cmd:str, capture_output=True, check:bool=False,
                   text=True, shell=True) -> Union[subprocess.CompletedProcess,
                                                   subprocess.CalledProcessError]:
    """Wraps subprocess.run. Return CompletedProcess or err obj."""

    try:
        data = subprocess.run(cmd,
                              capture_output=capture_output,
                              check=check,
                              text=text,
                              shell=shell
                             )
    except subprocess.CalledProcessError as e:
        data = e

    return data


def make_executable(sh_path:str) -> None:
    """Alternative to os.chmod(sh_path, stat.S_IXUSR): permission denied."""

    sh_path = Path(sh_path)
    cmd = f"chmod +x {str(sh_path)}"

    p = subprocess_run(cmd,
                       capture_output=False,
                       check=True)
    if isinstance(p, subprocess.CalledProcessError):
        logger.exception(f"Error in subprocess cmd 'chmod +x':\nException: {p}")
        raise

    return
