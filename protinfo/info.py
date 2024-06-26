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


from argparse import Namespace
import logging
from pathlib import Path
from protinfo import USER_MCCE, bio_parser, log_parser, run
from typing import Tuple, Union


logger = logging.getLogger(__name__)


def collect_info(pdb: Path, args: Namespace) -> Tuple[dict, Union[dict, None]]:
    """Return at least one dict holding info from Bio.PDB.PDBParser.
    The second dict is None when step1 cannot be run, otherwise
    it contains info from the parsed run.log file.
    Args:
      pdb (Path): File path of validated pdb.
      args (argparse.Namespace): Cli arguments (for creating step1 script).
    Returns:
      A 2-tuple of dicts when step1 can run, else (dict1, None).
    """

    DO_STEP1 = USER_MCCE is not None
    step1_d = None

    # Info from Bio.Parser:
    prot_d = bio_parser.info_input_prot(pdb)
    # if multimodels, no need to run step1:
    if "MultiModels" in prot_d["ParsedStructure"]:
        prot_d["Invalid"] = bio_parser.ERR_MULTI_MODELS
        DO_STEP1 = False

    if "Truncation" in prot_d["ParsedStructure"]:
        prot_d["Failed conversion"] = bio_parser.ERR_TRUNCATED_CONVERSION
        DO_STEP1 = False

    if DO_STEP1:
        run.do_step1(pdb, args)
        step1_d = log_parser.info_s1_log(pdb)

    return prot_d, step1_d


def get_pdb_report_lines(prot_d: dict, s1_d: Union[dict, None]) -> str:
    """Return the formated report lines for a pdbid for the two subsections
    in a pdb report with PDBParser info in prot_d and Step1 info in s1_d).
    Args:
      prot_d (dict): The dictionary of sections from bio parser.
      s1_d ([dict, None]): The dictionary of sections from step1 log parser.
    """

    if s1_d is None:
        dict_lst = [prot_d]
    else:
        dict_lst = [prot_d, s1_d]

    name = dict_lst[0].pop("Name")
    report = f"---\n# {name}\n"

    for i, subd in enumerate(dict_lst):
        # h2: section hdrs, ParsedStructure or MCCE.Step1
        h2 = list(subd.keys())[0]  # section hdrs: bioparser or mcce

        report = report + f"## {h2}\n"

        for k in subd[h2]:
            if not subd[h2][k]:
                continue

            if k in ["Chains", "Residues", "Waters", "Buried"]:
                if k == "Buried":
                    report = report + f"### {k} {bio_parser.BURIED_THR_MSG}"
                else:
                    report = report + f"### {k}: "
            else:
                report = report + f"### {k}:\n"

            for val in subd[h2][k]:
                if i == 0 and k == "Warnings":
                    warnstr = ""
                    d = subd[h2][k][val]
                    for w in d:
                        warnstr = warnstr + f"{w} ({', '.join(str(i) for i in d[w])}); "
                    report = report + f"  <strong><font color='red'>{val}</font></strong>: {warnstr}\n"

                elif i == 1 and isinstance(val, str) and (val.startswith("Generic") or val.startswith("Unloadable")):
                    report = report + f"  - <strong><font color='red'>{val}</font></strong>:\n"

                elif i == 1 and k == "Distance Clashes":
                    if i == 1 and isinstance(val, str) and val.startswith("Clashes"):
                        report = report + f"<details><summary>{val}</summary>\n"
                    elif i == 1 and isinstance(val, str) and val.endswith("end_clash"):
                        report = report + "</details>\n"
                    else:
                        report = report + f"  {val}\n"

                elif i == 0 and k == "Buried":
                    for d in val:
                        n, ids = val[d]
                        report = report + f"  - <strong>{d}</strong>: {n}: {', '.join(x for x in ids)}\n"

                elif isinstance(val, tuple):
                    ter, lst = val
                    report = report + f"  <strong>{ter}</strong>: {', '.join(lst)}\n"

                elif isinstance(val, list):
                    ter, lst = val
                    report = report + f"  <strong>{ter}</strong>: {', '.join(x for x in lst)}\n"

                else:
                    report = report + f"  {val}\n"

            report = report + "\n"

    report = report + "---\n"

    return report
