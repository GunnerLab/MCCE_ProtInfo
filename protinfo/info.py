#!/usr/bin/env python

"""
Module: info.py
Functions to process information about a protein:
 - Info about the input protein from Bio.PDB parser
 - Info from MCCE step1 run.log, debug.log when step1 can be run

 1. collect_info
 2. collect_info_lines
 3. save_report

"""

from argparse import Namespace
import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from collections import Counter, defaultdict
from pathlib import Path
from protinfo import USER_MCCE
from protinfo.log_parser import info_s1_log, blocks_specs
from protinfo import io_utils as iou, run
import time
from typing import Tuple, Union
import warnings


# error msg as fstrings:
ERR_MULTI_MODELS = "MCCE cannot handle multi-model proteins such as {}."
ERR_FETCH_EXISTING_FILE = """
The input_pdb parameter ({}) resolves to an existing file: to OVERWRITE it
with the biological assembly from a fresh download, remove the extension.
"""


def check_pdb_arg(input_pdb:str) -> Union[Path, str]:
    """Validate input_pdb str, which can be either a pdb id or a pdb file."""

    pdb = Path(input_pdb).resolve()
    if not pdb.exists():
        # if no extension, assume pdbid:
        if not pdb.suffix:
            return input_pdb

        raise FileNotFoundError(f"Not found: {pdb}")

    if pdb.suffix != ".pdb":
        raise TypeError(f"Not a valid extension: {pdb.suffix}")

    return pdb


def validate_input(args:Namespace) -> Path:
    """Validate args.input_pdb and args.fetch """
    
    pdb = check_pdb_arg(args.input_pdb)
    if isinstance(pdb, str):
        pdb = iou.maybe_download(args.input_pdb, args.fetch)
    else:
        if args.fetch:
            raise TypeError(ERR_FETCH_EXISTING_FILE.format(pdb))
    
    return pdb


def process_warnings(w:PDBConstructionWarning) -> dict:
    """Function to collapse warnings into categories."""

    warn_d = defaultdict(list) # output dict
    chn_d = defaultdict(list)
    unrec_l = []
    neg_l = []
    miss_l = []

    # Unpack the warnings list
    for i, info in enumerate(w):
        # extract str past "WARNING: ":
        line = info.message.args[0].removeprefix("WARNING: ") #[9:]

        if line.startswith("Chain"):
            newl = line[:-1].replace(' is discontinuous at ', '').split('line')
            chn_d[newl[0]].append("Line"+newl[1])
        elif line.startswith("Ignoring"):
            newl = line.removeprefix("Ignoring unrecognized ").split('at')
            unrec_l.append((newl[0].strip().capitalize(), newl[1].strip().capitalize()))
        elif line.startswith("Negative"):
            neg_l.append(line)
        elif line.startswith("Some atoms"):
            miss_l.append(line)
    
    if chn_d:
        warn_d["Discontinuity"] = dict(chn_d)
    if unrec_l:
        warn_d["Unrecognized records"] = unrec_l
    if neg_l:
        warn_d["Negative occupancy"] = neg_l
    if miss_l:
        warn_d["Missing"] = miss_l
        
    return dict(warn_d)


def info_input_prot(pdb:Path) -> dict:
    """Return information about input_pdb from Bio.PDB.PDBParser.
    The output dict structure follow the Bio.PDB.PDBParser
    structural hierarchy: Model, Chains, Residues, Atoms, along with
    a warnings section if any.
    """

    parser = Bio.PDB.PDBParser()

    dout = defaultdict(dict)
    dinner = defaultdict(list)

    pdbid = pdb.stem

    try:
        with warnings.catch_warnings(record=True, append=True) as w:
            warnings.simplefilter("always", category=PDBConstructionWarning)
            
            structure = parser.get_structure(pdbid, pdb.name)

    except Exception as ex:
        #raise ValueError(f"Could not parse {pdb.name}:\n{ex}")
        dout[pdbid]["ParsedStructure"]['ERROR'] = ex.args
        return dict(dout)
    
    n_models = len(structure)
    if n_models > 1:
        dinner["MultiModels"].append(n_models)
    else:
        s0 = structure[0]
        chains = list(s0.get_chains())
        n_chains = len(chains)
        n_res = 0
        n_hoh = 0
        cnames = []
        for c in chains:
            cname = c.get_id()
            cnames.append(cname)

            dc = defaultdict(dict)

            res = list(c.get_residues())
            n_res += len(res)
            n_hoh = len([r for r in res if r.resname.strip() == "HOH"])

            atoms = c.get_atoms()
            altlocs = []
            c = Counter()
            for a in atoms:
                alt = a.get_altloc()
                if alt != " ":
                    r = a.get_parent()
                    rname = r.get_resname()
                    rid = r.get_id()[1]
                    idx = a.get_serial_number()
                    c[(rname, rid, a.get_name(), idx)] += 1

            altlocs = [(itm, cnt) for itm, cnt in c.items() if cnt > 1]
            if altlocs:
                dc[cname] = altlocs
                dinner["Atoms.MultipleAltLocs"].append(dict(dc))

        dinner["Chains"].append((n_chains, cnames))
        dinner["Residues"].append(n_res)
        dinner["Waters"].append(n_hoh)

    dout[pdbid]["ParsedStructure"] = dict(dinner)
    if len(w):
        warn_d = process_warnings(w)
        dout[pdbid]["ParsedStructure"]["Warnings"] = warn_d

    return dict(dout)


def get_pdb_report_lines(pdbid:str, prot_d:dict, s1_d:Union[dict, None]) -> str:
    """Return the report lines for a pdbid for the two subsections in a
    pdb report with PDBParser info in prot_d and Step1 info in s1_d).
    Note:
    The input dictionaries are assumed to be the values of the
    complete dict 'filtered' by pdbid, the top level key, e.g.:
        prot_d = prot_info_d[pdbid]
    """

    if s1_d is None:
        dict_lst = [prot_d]
    else:
        dict_lst = [prot_d, s1_d]

    report = f"---\n# {pdbid}\n"
    for i, subd in enumerate(dict_lst):
        k0 = list(subd.keys())[0] # section hdrs: bioparser or mcce step1

        report = report + f"## {k0}\n"

        for k in subd[k0]:
            if not subd[k0][k]:
                continue

            report = report + f"### {k}\n"
            for val in subd[k0][k]:
                if i == 0 and k == "Warnings":
                    report = report + f"  - <strong><font color='red'>{val}</font> </strong>\n"
                    d = subd[k0][k][val]
                    for w in d:
                        report = report + f"    - {w}: {d[w]}\n"

                elif i == 1 and isinstance(val, str) and val.startswith("Generic"):
                    report = report + f"  - <strong><font color='red'>{val}</font> </strong>\n"

                elif i == 1 and k == "Distance Clashes:":
                    report = report + f"{val}\n"

                elif (isinstance(val, tuple) or isinstance(val, list)):
                    ter, lst = val
                    report = report + f"  * <strong>{ter} </strong> : {", ".join(lst)}\n"
                else:
                    report = report + f"  - {val}\n"

            report = report + "\n"

    report = report + "---\n"

    return report


def collect_info_lines(prot_d:dict,
                       s1_d:Union[dict, None]) -> str:
    """Transform the info in each dict into printable lines."""

    rpt_lines = ""
    for pdb in prot_d:
        rpt_lines = rpt_lines + get_pdb_report_lines(pdb,
                                                 prot_d[pdb],
                                                 s1_d[pdb])

    return rpt_lines


def save_report(report_lines:str,
                pdb_fp:Union[Path, None]=None,
                report_fp:Union[Path, None]=None):
    """Write and save the ProtInfo report.
    Args:
      report_lines (str): The lines to write
      pdb_fp (Union[Path, None], None): The pdb filepath if the
        report is for a single pdb. Cannot be None if report_fp is.
      report_fp (Union[Path, None], None): The filepath of the output
        report if pdb_fp is None (pdb_fp has precedence).
    """

    if pdb_fp is None and report_fp is None:
        raise ValueError("pdb_fp and report_fp cannot both be None.")
    
    if pdb_fp is not None:
       # (re)set report_fp
       report_fp = pdb_fp.parent.joinpath("ProtInfo.md")
   
    with open(report_fp, "w") as fo:
        fo.writelines(report_lines)
    
    return


def collect_info(pdb:Path) -> Tuple[dict, Union[dict, None]]:
    """Return at least one dict holding info from Bio.PDB.PDBParser.
    The second dict is None when step1 cannot be run, otherwise
    it contains info from the parsed run.log file.
    """

    DO_STEP1 = USER_MCCE is not None

    # Info from Bio.Parser:
    prot_d = info_input_prot(pdb)
    pdbid = pdb.stem
    # if multimodels, no need to run step1:
    if "MultiModels" in prot_d[pdbid]["ParsedStructure"]:
        prot_d[pdb.stem]["Invalid"] = ERR_MULTI_MODELS.format(pdbid)
        DO_STEP1 = False
        step1_d = None

    if DO_STEP1:
        #s1_start = time.time()
        run.do_step1(pdb)
        #elapsed = time.time() - s1_start
        #print(f"step1 took {elapsed:,.2f} s ({elapsed/60:,.2f} min).")
        time.sleep(2)
        step1_d = info_s1_log(pdb)

    return prot_d, step1_d


def get_single_pdb_report(args:Union[Namespace, dict]):
    """Get info and save report for a single pdb.
    This function is called by the cli.
    Expected keys in args: input_pdb [Path, str], fetch [bool].
    Workflow:
      1. validate_input
      2. collect_info in dicts
      3. collect_info_lines
      4. save_report
    """

    if isinstance(args, dict):
        args = Namespace(**args)

    pdb = validate_input(args)
    prot_d, step1_d = collect_info(pdb)
    report_lines = collect_info_lines(prot_d, step1_d)
    save_report(report_lines, pdb_fp=pdb)

    return
