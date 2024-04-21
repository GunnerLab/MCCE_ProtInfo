#!/usr/bin/env python

"""
Module: info.py
Functions to process information about a protein:
 - Info about the input protein from Bio.PDB parser
 - Info from MCCE step1 run.log, debug.log when step1 can be run
"""

from argparse import Namespace
import Bio.PDB
from collections import Counter, defaultdict
from dataclasses import dataclass
import pandas as pd
from protinfo import USER_MCCE, run
import protinfo.queries as qry
from pathlib import Path
import time
from typing import Union
import warnings


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
        dinner["MultiModels"].append(n_models)
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
                dc[cname] = altlocs
                dinner["SingleModel.Atoms.MultipleAltLocs"].append(dict(dc))

        dinner["SingleModel.Chains"].append((n_chains, cnames))
        dinner["SingleModel.Residues"].append(n_res)

    if len(w):
        # Unpack the warnings list
        for i, info in enumerate(w):
            # extract str past "WARNING: ":
            msg = info.message.args[0][9:]
            dinner["ParsedStructure.Warnings"].append(msg)

    dout[pdbid]["Input.ParsedStructure"] = dict(dinner)

    return dict(dout)


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


@dataclass
class LogHdr:
    idx: int
    hdr: str
    rpt_hdr: str
    line_start: Union[None, str] = None
    # special behaviour:
    # if list::full line, line is skipped if in list
    # if tuple::substring, line skipped if substr in line
    skip_lines: Union[None, list, tuple] = None
    debuglog: bool = False

    def has_debuglog(self, line:str, i:int=5) -> bool:
        """Set debug_log to True if line ends with 'saved in debug.log.
        i is meant to be the index key of the calling dict.
        '"""
        if i != 5:  # only known case: free cofactor stripping
            return
        
        if not self.debuglog:
            self.debuglog = line.endswith("saved in debug.log.")
        return
    

runlog_headers = [
    "   Rename residue and atom names...",
    "   Identify NTR and CTR...",
    "   Label backbone, sidechain and altLoc conformers...",
    "   Load pdb lines into data structure...",
    "   Strip free cofactors with SAS >   5%...",
    "   Check missing heavy atoms and complete altLoc conformers...",
    "   Find distance clash (<2.000)...",
    "   Make connectivity network ...",
    ]


def get_loghdr_specs(loghdrs:list) -> dict:
    """Return a dict of LogHdr classes.
    Note: Single quotes needed here for list of full lines?
    """

    all = defaultdict(dict)
    for i, hdr in enumerate(loghdrs, start=1):
        if i == 1:
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Renamed:",
                            line_start="   Renaming ",
                            )
        elif i == 2:
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Terminii:",
                            line_start="      Labeling ",
                            )
        elif i == 3:
            b3_exclude = (
                '   Creating temporary parameter file for unrecognized residue...',
                '   Trying labeling again...',
                '   Try delete this entry and run MCCE again',
                '   Error! premcce_confname()',
                ' is already loaded somewhere else.',
            )
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Labeling:",
                            line_start="      Labeling ",
                            skip_lines=b3_exclude
                            )
        elif i == 4:
            # keep as is until error found
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Load Structure:",
                            )
        elif i == 5:
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Free Cofactors:",
                            skip_lines=("free cofactors were stripped off in this round",
                                        "saved in debug.log.")
                           )
        elif i == 6:
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Missing Heavy atoms:",
                            line_start="   Missing heavy atom  ",
                            skip_lines=['   Missing heavy atoms detected.']
                           )
        elif i == 7:
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Distance Clashes:",
                            )
        elif i == 8:
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Connectivity:",
                            )
        else:
            # unknown
            all[i] = LogHdr(i, hdr,
                            rpt_hdr="Other:",
                            )

    return dict(all)


#blocks_dict = dict(list((r, hdr) for (r, hdr) in enumerate(runlog_headers, start=1)))
blocks_specs = get_loghdr_specs(runlog_headers)


class S1Log:
    """A class to parse run.log into sections, and
    process each one of them.
    """

    def __init__(self, pdb:Path) -> None:
        self.pdb = pdb
        self.pdbid = self.pdb.stem
        self.s1_dir = self.pdb.parent.joinpath("step1_run")
        self.check_debuglog_idx = [5]
        self.txt_blocks = self.get_blocks()


    def get_debuglog_species(self) -> str:

        fp =self.s1_dir.joinpath("debug.log")
        df = pd.read_csv(fp, sep=r"\s+", header=None, engine="python")
        txt = "Species and properties with assigned default values in debug.log:\n"
        for k in df[1].unique():
            txt = txt + f"  {k}: {list(df[df[1]==k][0].unique())}\n"

        return txt
    

    @staticmethod
    def process_content_block(content:list, lhdr:LogHdr) -> list:

        out = []
        skip = lhdr.skip_lines is not None
        change = lhdr.line_start is not None
        newtpl = None
        if lhdr.idx == 3:
            newtpl = []
        
        for line in content:
            if not line:
                continue

            if lhdr.idx == 3:
                if line.startswith("   Error! premcce_confname()"):
                    # add conf name:
                    newtpl.append(line.rsplit(maxsplit=1)[1])

            if lhdr.idx == 5:
                 # flag if debug.log was in line:
                if not lhdr.debuglog:
                    lhdr.has_debuglog(line)

            if skip:
                if isinstance(lhdr.skip_lines, tuple):
                    found = False
                    for t in lhdr.skip_lines:
                        found = found or t in line
                    if found:
                        continue
                else:
                    if line in lhdr.skip_lines:
                        continue
            
            if change:
                # remove common start:
                if line.startswith(lhdr.line_start):
                    line = line[len(lhdr.line_start):]

            out.append(line)
        
        # check if new tpl confs:
        if lhdr.idx == 3 and newtpl:
            out.append(f"Generic topology file created for: {newtpl}")
        
        return out
    
    def get_blocks(self):
        
        with open(self.s1_dir.joinpath("run.log")) as fp:
            text = fp.read()

        block_txt = {}
        for k in blocks_specs:
            lhdr = blocks_specs[k]
            content = extract_content_between_tags(text, lhdr.hdr).splitlines()

            if (lhdr.line_start is not None) or (lhdr.skip_lines is not None):
                content = self.process_content_block(content, lhdr)

            block_txt[k] = [line for line in content if line.strip()]
        
        if blocks_specs[5].debuglog:
            # add extra line:
            block_txt[5].append(self.get_debuglog_species())

        return block_txt


def info_s1_log(pdb:Path) -> dict:
    dout = {}
    silog = S1Log(pdb)
    #match struc in prot info dict:
    dout[pdb.stem] = {"MCCE.Info": silog.txt_blocks}

    return dout


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
            msg = f"""
            The input_pdb parameter ({pdb}) does not resolve to a file: If it is a pdbid
            and you meant to download its biological assembly, the fetch option must be True."""
            raise TypeError(msg)
        else:
            pdb = qry.get_rcsb_pd(pdb)
            if pdb is None:
                raise ValueError(f"Could not download {pdb} from rcsb.org.")
    else:
        if args.fetch:
            msg = f"""
            The input_pdb parameter ({pdb}) resolves to an existing file: to OVERWRITE it
            with the biological assembly from a fresh download, remove the extension."""
            raise TypeError(msg)

            
    pdb = check_pdb_arg(args.input_pdb)
    #assert isinstance(pdb, Path)

    # Info from Bio.Parser:
    input_info_d = info_input_prot(pdb)

    DO_STEP1 = USER_MCCE is not None

    # if multimodels, no need to run step1:
    if "MultiModels" in input_info_d[pdb.stem]["Input.ParsedStructure"]:
        #print(f"MCCE cannot handle multi-model proteins such as {pdb.stem}.")
        input_info_d[pdb.stem]["Input.Invalid"] = f"MCCE cannot handle multi-model proteins such as {pdb.stem}."

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
