#!/usr/bin/env python

"""
Module: io_utils

Module with generic (enough) functions related to loading and saving files.
"""

from argparse import Namespace
import logging
import pandas as pd
from pathlib import Path
import pickle
import subprocess
from typing import Union, Any


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
#......................................................................

def Pathok(pathname:str, check_fn:str="exists", raise_err=True) -> Union[Path, bool]:
    """Return path if check passed, else raise error.
    check_fn: one of 'exists', 'is_dir', 'is_file'.
    if raise_err=False, return False instead of err.
    """

    pathname = Path(pathname)
    if check_fn not in ['exists', 'is_dir', 'is_file']:
        check_fn = 'exists'

    # failure msg:
    if check_fn == 'exists':
        msg = f"Path not found: {pathname}"
    elif check_fn == 'is_dir':
        msg = f"Directory not found: {pathname}"
    elif check_fn == 'is_file':
        msg = f"Path is not a file: {pathname}"

    if not pathname.__getattribute__(check_fn)():
        if not raise_err:
            return False
        raise FileNotFoundError(msg)

    return pathname


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


def pct_completed(book_fpath:str) -> float:
    """Return the pct of runs that are completed or finished with error."""

    book_fp = Pathok(book_fpath)
    # 2 cmds: # processed; total count => 2 lines
    cmd = f"grep ' [ce]$' {book_fp} |wc -l; cat {book_fp} |wc -l"
    data = subprocess_run(cmd)

    if isinstance(data, subprocess.SubprocessError):
        logger.error(f"Error fetching pct completed.")
        return

    if not data.stdout.strip():
        logger.info("No data from book file")

    out = data.stdout.splitlines()

    return float(out[0])/float(out[1])


def make_executable(sh_path:str) -> None:
    """Alternative to os.chmod(sh_path, stat.S_IXUSR): permission denied."""

    sh_path = Pathok(sh_path)
    cmd = f"chmod +x {str(sh_path)}"

    p = subprocess_run(cmd,
                       capture_output=False,
                       check=True)
    if isinstance(p, subprocess.CalledProcessError):
        logger.exception(f"Error in subprocess cmd 'chmod +x':\nException: {p}")
        raise

    return


def tsv_to_df(fpath:str, index_col:Union[str, int]=0) -> Union[pd.DataFrame, None]:
    """Read a tab-separated file into a pandas.DataFrame.
    Return None upon failure.
    """
    fp = Pathok(fpath, raise_err=False)
    if not fp:
        logger.error(f"Not found: {fp}")
        return None
    
    return pd.read_csv(fp, index_col=index_col, sep="\t")


def df_to_tsv(df:pd.DataFrame, fpath:str, replace:bool=False) -> None:
    """Save a pandas.DataFrame as a tab separated file."""

    fp = Pathok(fpath, raise_err=False)
    if fp:
        if not replace:
            logger.error(f"File already exists: {fp}; {replace = }")
            return None
        fp.unlink()
    df.to_csv(fp, sep="\t")
    
    return


def pk_to_float(value) -> float:
    """Out of bound values become +/-8888 or 9999 (curve too sharp)
    during conversion to float.
    """
    try:
        v = float(value)
    except ValueError:
        if value.startswith("titra"): #tion curve too sharp"
            return 9999.
        oob = value[0]  # oob sign, > or <
        if oob == "<":
            v = -8888.
        else:
            v = 8888.

    return v


def get_file_header(fp:str) -> str:
    with open(fp) as f:
        for i, line in enumerate(f.readlines()):
            if i > 0:
                break
            hdr = line

    return hdr


def get_book_dirs_for_status(book_fpath:str, status:str="c") -> list:
    """Return a list of folder names from book_fp, the Q_BOOK file path,
    if their status codes match 'status', i.e. completed ('c', default),
    or errorneous ('e').
    """

    status = status.lower()
    if not status or status not in ["c", "e"]:
        logger.error("Invalid 'status'; choices are 'c' or 'e'")
        raise ValueError("Invalid 'status'; choices are 'c' or 'e'")

    book_fp = Pathok(book_fpath)
    book_dirs = []
    with open(book_fp) as book:
        for line in book:
            # select the portion preceding any appended comment
            rawtxt = line.strip().split("#")[0]
            fields = rawtxt.split()
            if len(fields) == 2:
                if fields[1].lower() == status:
                    book_dirs.append(fields[0])

    return book_dirs


def to_pickle(obj:Any, fp:str):
   pickle.dump(obj, open(fp, "wb"))
   return


def from_pickle(fp:str) -> Any:
   obj = pickle.load(open(fp, "rb"))
   return obj
