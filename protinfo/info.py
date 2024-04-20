#!/usr/bin/env python

"""
Module: info.py
Functions to process information about a protein:
 - Info about the input protein from Bio.PDB parser
 - Info from MCCE step1 run.log when step1 can be run
"""


from argparse import ArgumentParser, RawDescriptionHelpFormatter, Namespace
import Bio.PDB
from collections import Counter, defaultdict
from protinfo import USER_MCCE, run
import protinfo.queries as qry
from pathlib import Path
import time
from typing import Union
import warnings


DO_STEP1 = USER_MCCE is not None

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


def info_input_prot(pdb:Path) -> dict:
    """Return information about input_pdb from Bio.PDB.PDBParser."""

    parser = Bio.PDB.PDBParser()

    dout = defaultdict(dict)
    dinner = defaultdict(list)

    pdbid = pdb.stem

    with warnings.catch_warnings(record=True, append=True) as w:
        warnings.simplefilter("always",
                              category=Bio.PDB.PDBExceptions.PDBConstructionWarning)
        structure = parser.get_structure(pdbid, pdb.name)

    n_models = len(structure)
    if n_models > 1:
        dinner["MultiModels"].append((pdbid, n_models))
    else:
        s0 = structure[0]
        chains = list(s0.get_chains())
        n_chains = len(chains)
        n_res = 0
        cnames = []
        for c in chains:
            cname = c.get_id()
            cnames.append(cname)

            dc = defaultdict(dict)

            res = list(c.get_residues())
            n_res += len(res)

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
                dc[(pdbid,cname)] = altlocs
                dinner["SingleModel.Atoms.MultipleAltLocs"].append(dict(dc))

        dinner["SingleModel.Chains"].append((pdbid, n_chains, cnames))
        dinner["SingleModel.Residues"].append((pdb.stem, n_res))

    if len(w):
        # Unpack the warnings list
        for i, info in enumerate(w):
            msg = info.message.args
            dinner["ParsedStructure.Warnings"].append((pdbid, msg))

    dout[pdbid]["Input.ParsedStructure"] = dict(dinner)

    return dict(dout)


runlog_headers = ["   Rename residue and atom names...",
"   Identify NTR and CTR...",
"   Label backbone, sidechain and altLoc conformers...",
"   Load pdb lines into data structure...",
"   Strip free cofactors with SAS >   5%...",
"   Check missing heavy atoms and complete altLoc conformers...",
"   Find distance clash (<2.000)...",
"   Make connectivity network ...",
]

blocks_dict = dict(list((r, hdr) for (r, hdr) in enumerate(runlog_headers, start=1)))


def extract_content_between_tags(text:str, tag1:str, tag2:str="   Done"):
  """Extracts the content between two string tags in a text.
  Args:
    text: A text.
    tag1: The first tag.
    tag2: The second tag, default: "   Done".

  Returns:
    A string containing the content between the two headers.
  """

  start_pos = text.find(tag1) + len(tag1)
  end_pos = text.find(tag2, start_pos)

  return text[start_pos:end_pos]



def info_s1_log(pdb:Path) -> dict:

    pdbname = pdb.stem
    s1_dir = run.get_s1_dir(pdb)

    with open(s1_dir.joinpath("run.log")) as fp:
        text = fp.read()

    block_txt = dict()

    for k in blocks_dict:
        content = extract_content_between_tags(text, blocks_dict[k])
        block_txt[k] = [line+"\n" for line in content.splitlines() if line.strip()]

    return block_txt



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
    is_pdbid = isinstance(pdb, str)

    if is_pdbid:
        if not args.fetch:
            msg = f"input_pdb ({pdb}) is not a file. If it is a valid "
            msg = msg + "pdb id and you meant to download it, the fetch option must be True."
            raise TypeError(msg)
        else:
            pdb = qry.get_rcsb_pd(pdb)
            if pdb is None:
                raise ValueError(f"Could not download {pdb} from rcsb.org.")
            
    pdb = check_pdb_arg(args.input_pdb)
    assert isinstance(pdb, Path)

    # Info from Bio.Parser:
    input_info_d = info_input_prot(pdb)
    # if multimodels, no need to run step1:
    if "MultiModels" in input_info_d[pdb.stem]["Input.ParsedStructure"]:
        print(f"Cannot run step1 on multi-model protein {pdb.stem}")
        DO_STEP1 = False
        s1_info_d = None

    if DO_STEP1:
        s1_start = time.time()
        run.do_step1(pdb)
        elapsed = time.time() - s1_start
        print(f"step1 took {elapsed:,.2f} s ({elapsed/60:,.2f} min).")
        #time.sleep(3)
        s1_info_d = info_s1_log(pdb)

    write_report(input_info_d, s1_info_d)

    return
