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

from Bio.PDB import PDBParser, parse_pdb_header
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from collections import Counter, defaultdict
import logging
from pathlib import Path
import warnings


logger = logging.getLogger(__name__)


ERR_MULTI_MODELS = "MCCE cannot handle multi-model proteins."
ERR_TRUNCATED_CONVERSION = """
This pdb was truncated during the .cif to .pdb format conversion
because the number of atoms exceeds 99,999."""


BURIED_THRESH = 0.05  # mcce default; res with sasa < this are buried.
BURIED_THR_MSG = f"(using default mcce SASA threshold of {BURIED_THRESH:.0%}):\n"
# This ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7817970/
# has benchmarked a SASA threshold of 20%.


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

        if line.startswith("Chain "):
            newl = line[6:-1].replace(" is discontinuous at ", "").split("line")
            chn_d[newl[0]].append(int(newl[1]))

        elif line.startswith("Ignoring"):
            newl = line.removeprefix("Ignoring unrecognized ").split("at")
            unrec_l.append((newl[0].strip().capitalize(), newl[1].strip().capitalize()))
        elif line.startswith("Negative"):
            neg_l.append(line)
        elif line.startswith("Some atoms"):
            miss_l.append(line)

    if chn_d:
        warn_d["Chain discontinuity"] = dict(chn_d)
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

    pdbid = pdb.stem
    parser = PDBParser()
    dout = defaultdict(dict)
    dinner = defaultdict(list)

    try:
        with warnings.catch_warnings(record=True) as w:
            # catch_warnings w/o happen (3.10): same output
            # possibly due to specified filter:
            warnings.simplefilter("always", category=PDBConstructionWarning)

            structure = parser.get_structure(pdbid, pdb.name)

    except Exception as ex:
        # dout[pdbid]["ParsedStructure"] = {"ERROR": ex.args + " (Possibly not biopython related.)"}
        dout["ParsedStructure"] = {"ERROR": ex.args + " (Possibly not biopython related.)"}
        logger.error(f"ERROR: {ex.args} - (Possibly not biopython related.)")
        return dict(dout)

    pdb_hdr_d = parse_pdb_header(pdb)
    protname = pdb_hdr_d.get("name")  # == "" if no header

    # check if header has truncation warning:
    note_hdr = pdb_hdr_d.get("head")
    if note_hdr:
        if "truncated" in note_hdr and note_hdr.endswith("true"):
            dinner["Truncation"].append(f"WARNING: {note_hdr}")

    n_models = len(structure)
    if n_models > 1:
        dinner["MultiModels"].append(n_models)
    else:
        s0 = structure[0]
        sasa = ShrakeRupley()
        sasa.compute(s0, level="R")

        chains = list(s0.get_chains())
        n_chains = len(chains)
        n_res = 0
        n_hoh = 0
        n_het = 0
        cnames = []
        waters = []
        buried_wat = []
        heteros = []
        buried_het = []

        for c in chains:
            cname = c.get_id()
            cnames.append(cname)

            dc = defaultdict(dict)

            res = list(c.get_residues())
            n_res += len(res)
            heteros = [r for r in res if r.get_id()[0].startswith("H_")]
            for het in heteros:
                if het.sasa < BURIED_THRESH:
                    hetid = het.get_id()
                    buried_het.append(f"{c.id} {hetid[0][2:]} {hetid[1]}")
            n_het += len(heteros)

            waters = [r for r in res if r.resname.strip() == "HOH"]
            n_hoh += len(waters)
            for wat in waters:
                if wat.sasa < BURIED_THRESH:
                    buried_wat.append(f"{c.id} {wat.get_id()[1]}")

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

        if buried_wat:
            dinner["Buried"].append({"Waters": (len(buried_wat), buried_wat)})
        if buried_het:
            dinner["Buried"].append({"Heteros": (len(buried_het), buried_het)})

    if protname:
        # dout[pdbid]["Name"] = f"{pdbid.upper()} :: {protname.title()}"
        dout["Name"] = f"{pdbid.upper()} :: {protname.title()}"
    else:
        dout["Name"] = f"{pdbid.upper()}"

    # dout[pdbid]["ParsedStructure"] = dict(dinner)
    dout["ParsedStructure"] = dict(dinner)
    if len(w):
        warn_d = process_warnings(w)
        # dout[pdbid]["ParsedStructure"]["Warnings"] = warn_d
        dout["ParsedStructure"]["Warnings"] = warn_d

    return dict(dout)
