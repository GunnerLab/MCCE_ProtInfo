#!/usr/bin/env python

"""
Module: run.py
Functions to launch mcce step1.py
"""

import protinfo.io_utils as iou
from pathlib import Path
import shutil
import subprocess


S1_SH = """#!/bin/bash

step1.py prot.pdb
sleep 1
"""


def get_s1_dir(pdb_fp: Path) -> Path:
    """Return the step1 output dir path for pdb_fp."""

    return pdb_fp.parent.joinpath("step1_run").resolve()


def write_script(dest_dirpath: Path):
    """Write an excecutable bash script to run mcce step1."""

    sh_path = dest_dirpath.joinpath("s1.sh")
    with open(sh_path, "w") as fh:
        fh.write(S1_SH)

    # make script executable:
    iou.make_executable(sh_path)

    return


def prep_script(s1_dir: Path):

    if s1_dir.exists():
        shutil.rmtree(s1_dir)
    s1_dir.mkdir()

    write_script(s1_dir)

    return


def run_step1(s1_dir: Path) -> None:
    """Run step1 in s1_dir."""

    subprocess.Popen(f"{s1_dir}/s1.sh",
                     cwd=str(s1_dir),
                     close_fds=True,
                     stdout=open(f"{s1_dir}/run.log", "w")
                     )

    return


def do_step1(pdb_fp: Path):
    """Main function"""

    # output_dir:
    s1_dir = get_s1_dir(pdb_fp)
    prep_script(s1_dir)
    # setup prot.pdb:
    prot = s1_dir.joinpath("prot.pdb")
    prot.symlink_to(f"../{pdb_fp.name}")
    # launch step1:
    run_step1(s1_dir)

    return
