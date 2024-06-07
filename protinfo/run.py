#!/usr/bin/env python

"""
Module: run.py
Functions to launch mcce step1.py
"""


from argparse import Namespace
import logging
from protinfo.io_utils import make_executable
from pathlib import Path
import shutil
import subprocess


logger = logging.getLogger(__name__)
# logger.setLevel(logging.WARNING)


CUSTOM_S1_SH = """#!/bin/bash

step1.py prot.pdb {wet}{noter}{d}{e}{u}
"""

s1_defaults = {
    "wet": False,
    "noter": False,
    "d": 4.0,
    "e": "mcce",
    "u": "",
}


# maping of cli arg names to mcce steps args:
cli_to_mcce_opt = {"wet": "dry"}


def cli_args_to_dict(sh_args: Namespace) -> dict:
    """Only return step1 args."""

    excluded_keys = ["pdb", "fetch"]
    d_args = {k: v for k, v in vars(sh_args).items() if k not in excluded_keys}

    return d_args


def populate_sh_template(job_args: Namespace) -> str:
    """Return the custom template string filled with values."""

    d_args = cli_args_to_dict(job_args)
    # add missing keys needed by template from s1_defaults:
    d_args.update(((k, v) for k, v in s1_defaults.items() if k not in d_args))

    d_all = {}
    # note: trailing spaces needed

    # special cases:
    v = d_args.pop("wet")
    d_all["wet"] = "" if v else "--dry "
    v = d_args.pop("noter")
    d_all["noter"] = "--noter " if v else ""

    # all remaining options:
    for k in d_args:
        v = d_args.get(k, "")
        if str(v) == str(s1_defaults[k]):
            d_all[k] = ""
        else:
            d_all[k] = f"-{cli_to_mcce_opt.get(k, k)} {v} "

    body = CUSTOM_S1_SH.format(**d_all)

    return body


def get_s1_dir(pdb_fp: Path) -> Path:
    """Return the step1 output dir path for pdb_fp."""

    return pdb_fp.parent.joinpath("step1_run").resolve()


def write_script(dest_dirpath: Path, sh_str: str):
    """Write an executable bash script to run mcce step1."""

    sh_path = dest_dirpath.joinpath("s1.sh")
    sh_path.write_text(sh_str)

    make_executable(sh_path)

    return


def prep_rundir(s1_dir: Path):
    if s1_dir.exists():
        shutil.rmtree(s1_dir)
    s1_dir.mkdir()

    return


def run_step1(s1_dir: Path) -> None:
    """Run step1 in s1_dir."""

    subprocess.Popen(f"{s1_dir}/s1.sh", cwd=str(s1_dir), close_fds=True, stdout=open(f"{s1_dir}/run.log", "w"))

    return


def do_step1(pdb_fp: Path, args: Namespace):
    """Main function"""

    # output_dir:
    s1_dir = get_s1_dir(pdb_fp)
    prep_rundir(s1_dir)
    sh_str = populate_sh_template(args)
    write_script(s1_dir, sh_str)
    # setup prot.pdb as soft link:
    prot = s1_dir.joinpath("prot.pdb")
    prot.symlink_to(f"../{pdb_fp.name}")
    # launch step1:
    run_step1(s1_dir)

    return
