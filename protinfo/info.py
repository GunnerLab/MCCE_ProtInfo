#!/usr/bin/env python

"""
Module: info.py

Holds functions to gather information about a protein obtained
from bio_parser.py and log_parser.py:
 - Info about the input protein from Bio.PDB parser
 - Info from MCCE step1 run.log, debug.log when step1 can be run

Main steps:
  1. validate_pdb_inputs
  2. queries.get_rcsb_pdb(pdb):
     download bio-assembly if pdb is id
  3. collect_info in two dicts:
     first: from Bio.parser, second: from run.log parser
  4. collect_info_lines from the two dicts
  5. save_report
"""


import logging
from pathlib import Path
from protinfo import USER_MCCE, bio_parser, log_parser, run
from time import sleep
from typing import Tuple, Union


logger = logging.getLogger(__name__)
# logger.setLevel(logging.WARNING)


def collect_info(pdb: Path) -> Tuple[dict, Union[dict, None]]:
    """Return at least one dict holding info from Bio.PDB.PDBParser.
    The second dict is None when step1 cannot be run, otherwise
    it contains info from the parsed run.log file.
    """

    DO_STEP1 = USER_MCCE is not None
    step1_d = None

    # Info from Bio.Parser:
    prot_d = bio_parser.info_input_prot(pdb)
    pdbid = pdb.stem
    # if multimodels, no need to run step1:
    if "MultiModels" in prot_d[pdbid]["ParsedStructure"]:
        prot_d[pdb.stem]["Invalid"] = bio_parser.ERR_MULTI_MODELS
        DO_STEP1 = False

    if "Truncation" in prot_d[pdbid]["ParsedStructure"]:
        prot_d[pdb.stem]["Failed conversion"] = bio_parser.ERR_TRUNCATED_CONVERSION
        DO_STEP1 = False

    if DO_STEP1:
        run.do_step1(pdb)
        sleep(2)
        step1_d = log_parser.info_s1_log(pdb)
        # print("step1_d from collect info:\n", step1_d)

    return prot_d, step1_d


def get_pdb_report_lines(pdbid: str, prot_d: dict, s1_d: Union[dict, None]) -> str:
    """Return the formated report lines for a pdbid for the two subsections
    in a pdb report with PDBParser info in prot_d and Step1 info in s1_d).
    Args:
      pdbid (str): The pdb id.
      prot_d (dict): The dictionary of sections from bio parser.
      s1_d ([dict, None]): The dictionary of sections from step1 log parser.
    """

    if s1_d is None:
        dict_lst = [prot_d]
    else:
        dict_lst = [prot_d, s1_d]

    report = f"---\n# {pdbid}\n"
    name = dict_lst[0].get("Name")
    if name is not None:
        report = f"---\n# {pdbid} :: {name}\n"
        _ = dict_lst[0].pop("Name")

    for i, subd in enumerate(dict_lst):
        # k0: section hdrs, ParsedStructure or MCCE.Step1
        k0 = list(subd.keys())[0]  # section hdrs: bioparser or mcce

        report = report + f"## {k0}\n"

        for k in subd[k0]:
            if not subd[k0][k]:
                continue

            report = report + f"### {k}\n"
            for val in subd[k0][k]:
                if i == 0 and k == "Warnings":
                    report = report + f"  * <strong><font color='red'>{val}</font> </strong>\n"
                    d = subd[k0][k][val]
                    for w in d:
                        report = report + f"    - {w}: {d[w]}\n"

                elif i == 1 and isinstance(val, str) and val.startswith("Generic"):
                    report = report + f"  - <strong><font color='red'>{val}</font> </strong>\n"

                elif i == 1 and k == "Distance Clashes":
                    if i == 1 and isinstance(val, str) and val.startswith("Clashes"):
                        report = report + f"<details><summary>{val}</summary>\n"
                    elif i == 1 and isinstance(val, str) and val.endswith("end_clash"):
                        report = report + "</details>\n"
                    else:
                        report = report + f"  {val}\n"

                elif isinstance(val, tuple) or isinstance(val, list):
                    ter, lst = val
                    report = report + f"  * <strong>{ter} </strong> : {', '.join(lst)}\n"
                else:
                    report = report + f"  - {val}\n"

            report = report + "\n"

    report = report + "---\n"

    return report


def collect_info_lines(prot_d: dict, s1_d: Union[dict, None]) -> str:
    """Transform the info in each dict into printable lines."""

    rpt_lines = ""
    for pdb in prot_d:
        if s1_d is not None:
            rpt_lines = rpt_lines + get_pdb_report_lines(pdb, prot_d[pdb], s1_d[pdb])
        else:
            rpt_lines = rpt_lines + get_pdb_report_lines(pdb, prot_d[pdb], None)

    return rpt_lines
