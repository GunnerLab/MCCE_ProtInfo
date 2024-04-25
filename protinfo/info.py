#!/usr/bin/env python

"""
Module: info.py
Functions to process information about a protein:
 - Info about the input protein from Bio.PDB parser
 - Info from MCCE step1 run.log, debug.log when step1 can be run
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
from typing import Union
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


def process_warnings(w:PDBConstructionWarning) -> dict:
    """Function to collapse warnings into categories."""

    warn_d = defaultdict(list)
    chn_d = defaultdict(list)
    unrec_l = []
    neg_l = []
    miss_l = []

    # Unpack the warnings list
    for i, info in enumerate(w):
        # extract str past "WARNING: ":
        line = info.message.args[0][9:]

        if line.startswith("Chain"):
            newl = line[:-1].replace(' is discontinuous at ', '').split('line')
            chn_d[newl[0]].append("Line"+newl[1])
        elif line.startswith("Ignoring"):
            newl = line.removeprefix("Ignoring unrecognized ").split('at')
            unrec_l.append((newl[0].strip().capitalize(), new_l[1].strip().capitalize()))
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

    with warnings.catch_warnings(record=True, append=True) as w:
        warnings.simplefilter("always", category=PDBConstructionWarning)
        structure = parser.get_structure(pdbid, pdb.name)

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
                dinner["SingleModel.Atoms.MultipleAltLocs"].append(dict(dc))

        dinner["SingleModel.Chains"].append((n_chains, cnames))
        dinner["SingleModel.Residues"].append(n_res)
        dinner["SingleModel.Waters"].append(n_hoh)

    dout[pdbid]["Input.ParsedStructure"] = dict(dinner)
    if len(w):
        warn_d = process_warnings(w)
        #dout[pdbid]["Input.ParsedStructure"] = dict(dinner)
        dout[pdbid]["Input.ParsedStructure"]["ParsedStructure.Warnings"] = warn_d

    return dict(dout)


def write_report(input_info_d:dict=None, s1_info_d:dict=None):
    # TODO: write_report

    if input_info_d is not None:
        print(input_info_d)

    if s1_info_d is not None:
        print(s1_info_d)

    return


def main(args):

    if isinstance(args, dict):
        args = Namespace(**args)

    pdb = check_pdb_arg(args.input_pdb)
    if isinstance(pdb, str):
        pdb = iou.maybe_download(args.input_pdb, args.fetch)
    else:
        if args.fetch:
            raise TypeError(ERR_FETCH_EXISTING_FILE.format(pdb))

    assert isinstance(pdb, Path)

    DO_STEP1 = USER_MCCE is not None

    # Info from Bio.Parser:
    input_info_d = info_input_prot(pdb)

    # if multimodels, no need to run step1:
    if "MultiModels" in input_info_d[pdb.stem]["Input.ParsedStructure"]:
        input_info_d[pdb.stem]["Input.Invalid"] = ERR_MULTI_MODELS.format(pdb.stem)
        DO_STEP1 = False
        s1_info_d = None

    if DO_STEP1:
        #s1_start = time.time()
        run.do_step1(pdb)
        #elapsed = time.time() - s1_start
        #print(f"step1 took {elapsed:,.2f} s ({elapsed/60:,.2f} min).")
        time.sleep(2)
        s1_info_d = info_s1_log(pdb)

    write_report(input_info_d, s1_info_d)

    return
