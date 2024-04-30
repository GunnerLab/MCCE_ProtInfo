#!/usr/bin/env python

"""
Module: io_utils
"""

import logging
from pathlib import Path
from protinfo import queries as qrs
import subprocess
from typing import Union


logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def maybe_download(pdbid: str) -> Path:

    pdb = qrs.get_rcsb_pdb(pdbid)
    if pdb is None:
        logger.error(f"Could not download {pdbid} from rcsb.org.")

    return pdb


def subprocess_run(cmd: str,
                   capture_output: bool = True,
                   check: bool = False,
                   text: bool = True,
                   shell: bool = True) -> Union[subprocess.CompletedProcess,
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


def make_executable(sh_path: str) -> None:
    """Alternative to os.chmod(sh_path, stat.S_IXUSR): permission denied."""

    sh_path = Path(sh_path)
    cmd = f"chmod +x {str(sh_path)}"

    try:
        subprocess_run(cmd,
                       capture_output=False,
                       check=True)
    except subprocess.CalledProcessError:
        logger.exception("Error in subprocess cmd 'chmod +x'")
        raise

    return
