#!/usr/bin/env python

"""
Module: bio_parser.py

Holds functions to process information about a protein
returned from from Bio.PDB parser

Functions:
 - info_input_prot(pdb: Path) -> dict:
 - process_warnings(w: PDBConstructionWarning) -> dict:

 TODO:
 feat: Get neighbors of HETATMs not waters.
"""


import Bio.PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from collections import Counter, defaultdict
import logging
from pathlib import Path
import warnings


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


ERR_MULTI_MODELS = "MCCE cannot handle multi-model proteins."
ERR_TRUNCATED_CONVERSION = """
This pdb was truncated during the .cif to .pdb format conversion
because the number of atoms exceeds 99,999."""


def process_warnings(w: PDBConstructionWarning) -> dict:
    """Function to collapse warnings into categories."""

    warn_d = defaultdict(list)  # output dict
    chn_d = defaultdict(list)
    unrec_l = []
    neg_l = []
    miss_l = []

    # Unpack the warnings list
    for i, info in enumerate(w):
        # extract str past "WARNING: ":
        line = info.message.args[0].removeprefix("WARNING: ")

        if line.startswith("Chain"):
            newl = line[:-1].replace(" is discontinuous at ", "").split("line")
            chn_d[newl[0]].append("Line" + newl[1])
        elif line.startswith("Ignoring"):
            newl = line.removeprefix("Ignoring unrecognized ").split("at")
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


def info_input_prot(pdb: Path) -> dict:
    """Return information about 'pdb' from Bio.PDB.PDBParser.
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
        dout[pdbid]["ParsedStructure"]["ERROR"] = ex.args
        return dict(dout)

    pdb_hdr_d = Bio.PDB.parse_pdb_header(pdb)
    protname = pdb_hdr_d.get("name")
    if protname is not None:
        dinner["Name"].append(protname.title())
    # check if header has truncation warning:
    note_hdr = pdb_hdr_d.get("head")
    if note_hdr is not None:
        if "truncated" in note_hdr and note_hdr.endswith("true"):
            dinner["Truncation"].append(f"WARNING: {note_hdr}")

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