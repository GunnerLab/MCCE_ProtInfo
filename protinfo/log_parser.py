#!/usr/bin/env python

"""
Module: log_parser.py

Functions and classes to process information from mcce run.log.
For this app purpose, only the output from step1 is considered, which
is handled by the class RunLog1.
"""

from collections import defaultdict
from dataclasses import dataclass
import pandas as pd
from pathlib import Path
import logging
from typing import Union


logger = logging.getLogger(__name__)
# logger.setLevel(logging.WARNING)


def extract_content_between_tags(text: str, tag1: str, tag2: str = "   Done") -> str:
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


def get_pubchem_compound_link(compound_id: str) -> str:
    """Return the unvalidated link of the PubChem page for compund_id.
    subtance tab
    """
    if compound_id:
        url_fstr = "https://pubchem.ncbi.nlm.nih.gov/#query={}&tab=substance"
        return url_fstr.format(compound_id.upper())
    else:
        return ""


@dataclass
class LogHdr:
    """Dataclass to store information and line transformations for
    each 'processing block' returned in run.log after step1 has run.
    Each processing block starts with a header line (hdr) and ends
    with a common '   Done' line.

    Attributes:
        idx (int, start=1): Index of the block in order of appearance
        hdr (str): header line
        rpt_hdr (str): Corresponding header in the report
        line_start (Union[None, str]): Remove the substring
          from the start of the line
        skip_lines (Union[None, list, tuple]): List of lines to skip
          OR skip a line if substr in line when skip_lines is a tuple.
        debuglog (bool): Whether a line in the given block mentions
          'debug.log'; the debug.log file will then be parsed.
    """

    idx: int
    hdr: str
    rpt_hdr: str
    line_start: Union[None, str] = None
    skip_lines: Union[None, list, tuple] = None
    debuglog: bool = False

    def has_debuglog(self, line: str, calling_key: int = 5) -> bool:
        """Set debuglog proprerty to True if line ends with
        'saved in debug.log.'
        calling_key is meant to be the index key of the calling dict.
        The only known case where debug.log is mentioned is in 'block 5'
        of step1.py: 'free cofactor stripping' (see runlog_headers list).
        '"""

        # condition to be removed if debug.log found in other steps:
        if calling_key != 5:
            return

        if not self.debuglog:
            self.debuglog = line.endswith("saved in debug.log.")
        return


# list of run.log headers returned by mcce step1:
runlog1_headers = [
    "   Rename residue and atom names...",
    "   Identify NTR and CTR...",
    "   Label backbone, sidechain and altLoc conformers...",
    "   Load pdb lines into data structure...",
    "   Strip free cofactors with SAS >   5%...",  # 5
    "   Check missing heavy atoms and complete altLoc conformers...",
    "   Find distance clash (<2.000)...",
    "   Make connectivity network ...",
]


def get_log1_specs(loghdrs: list) -> dict:
    """Return a dict of LogHdr classes for processing step1 sections
    in run.log.
    """

    all = defaultdict(dict)
    for i, hdr in enumerate(loghdrs, start=1):
        if i == 1:
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Renamed",
                line_start="   Renaming ",
            )
        elif i == 2:
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Termini",
                line_start="      Labeling ",
            )
        elif i == 3:
            b3_exclude = (
                "   Creating temporary parameter file for unrecognized residue...",
                "   Trying labeling again...",
                "   Try delete this entry and run MCCE again",
                "   Error! premcce_confname()",
                " is already loaded somewhere else.",
            )
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Labeling",
                line_start="      Labeling ",
                skip_lines=b3_exclude,
            )
        elif i == 4:
            # keep as is until error found
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Load Structure",
            )
        elif i == 5:
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Free Cofactors",
                skip_lines=(
                    "free cofactors were stripped off in this round",
                    "saved in debug.log.",
                ),
            )
        elif i == 6:
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Missing Heavy Atoms",
                line_start="   Missing heavy atom  ",
                skip_lines=["   Missing heavy atoms detected."],
            )
        elif i == 7:
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Distance Clashes",
            )
        elif i == 8:
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Connectivity",
            )
        else:
            # unknown
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Other",
            )

    return dict(all)


blocks_specs = get_log1_specs(runlog1_headers)


class RunLog1:
    """A class to parse mcce run.log into sections pertaining to step1, and
    process each one of them into a simplified output.
    """

    def __init__(self, pdb: Path) -> None:
        self.pdb = pdb.resolve()
        self.pdbid = self.pdb.stem
        self.s1_dir = self.pdb.parent.joinpath("step1_run")
        # id of block with debug.log mentions, if any:
        self.check_debuglog_idx = [5]
        self.txt_blocks = self.get_blocks()

    def get_debuglog_species(self) -> str:
        fp = self.s1_dir.joinpath("debug.log")
        df = pd.read_csv(fp, sep=r"\s+", header=None, engine="python")
        txt = "Species and properties with assigned default values in debug.log:\n"
        for k in df[1].unique():
            txt = txt + f"{k}: {list(df[df[1] == k][0].unique())}\n"

        return txt

    @staticmethod
    def process_content_block(content: list, loghdr: LogHdr) -> list:
        out = []
        skip = loghdr.skip_lines is not None
        change = loghdr.line_start is not None
        newtpl = None

        # section that lists new.tpl creation:
        if loghdr.idx == 3:
            newtpl = []

        for line in content:
            if not line:
                continue

            if loghdr.idx == 3:
                if line.startswith("   Error! premcce_confname()"):
                    # add conf name & link:
                    conf = line.rsplit(maxsplit=1)[1]
                    if conf.startswith("_"):
                        pchem = get_pubchem_compound_link(conf[1:])
                    else:
                        pchem = get_pubchem_compound_link(conf)
                    newtpl.append(f"{conf}::  {pchem}")

                # TODO
                # elif line.startswith("   Error! The following atoms of residue"):

            if loghdr.idx == 5:
                # flag if 'debug.log' found in line:
                if not loghdr.debuglog:
                    loghdr.has_debuglog(line)

                if line.startswith("   Total deleted cofactors"):
                    if int(line.rsplit(maxsplit=1)[1][:-1]) != 0:
                        line = line.strip()
                    else:
                        continue

            if skip:
                if isinstance(loghdr.skip_lines, tuple):
                    found = False
                    for t in loghdr.skip_lines:
                        found = found or (t in line)
                    if found:
                        continue
                else:
                    if line in loghdr.skip_lines:
                        continue

            if change:
                # remove common start:
                if line.startswith(loghdr.line_start):
                    line = line.removeprefix(loghdr.line_start)

            out.append(line)

        # check if new tpl confs:
        if loghdr.idx == 3 and newtpl:
            out.append("Generic topology file created for:")
            for t in newtpl:
                out.append(t)

        return out

    def get_blocks(self) -> dict:
        """Extract 'processing blocks' from run.log file."""

        with open(self.s1_dir.joinpath("run.log")) as fp:
            text = fp.read()

        block_txt = {}
        for k in blocks_specs:
            lhdr = blocks_specs[k]
            rpt_k = lhdr.rpt_hdr
            content = extract_content_between_tags(text, lhdr.hdr).splitlines()
            if k == 1:
                content = sorted(content)

            if (lhdr.line_start is not None) or (lhdr.skip_lines is not None):
                content = self.process_content_block(content, lhdr)
            # debug
            if k == 6:
                print("get_blocks, k=6 contents:\n{content}")
            block_txt[rpt_k] = [line for line in content if line.strip()]

        # process termini; group res into NTR, CTR
        b2_hdr = blocks_specs[2].rpt_hdr
        if block_txt[b2_hdr]:
            termi = defaultdict(list)
            for line in block_txt[b2_hdr]:
                i = line.index('"', 3) + 1
                termi[line[-3:]].append(line[:i])
            block_txt[b2_hdr] = []
            for k in termi:
                block_txt[b2_hdr].append((k, termi[k]))

        if blocks_specs[5].debuglog:
            # add extra line for each species found:
            rk = blocks_specs[5].rpt_hdr
            for line in self.get_debuglog_species().splitlines():
                block_txt[rk].append(line)

        # collapse dist clashes block 7:
        if block_txt["Distance Clashes"]:
            new7 = []
            new7.append("Clashes found")
            for d in block_txt["Distance Clashes"]:
                new7.append(d.strip())
            new7.append("end_clash")  # tag for formatting section

            block_txt["Distance Clashes"] = new7

        return block_txt


def filter_heavy_atm_section(pdb: Path, s1_info_d: dict) -> dict:
    """Process the 'Missing Heavy Atoms' section to remove
    lines for missing backbone atoms of terminal residues.
    """

    # term values are [2-tuples]
    term = s1_info_d[pdb.stem]["MCCE.Step1"]["Termini"]
    heavy = s1_info_d[pdb.stem]["MCCE.Step1"]["Missing Heavy Atoms"]
    if len(heavy) > 1:
        _ = heavy.pop(-1)
        # == Ignore warning messages if they are in the terminal residues

    hvy_lst = []
    for line in heavy:
        conf, res = line.split(" in ")
        is_bkb = conf.rsplit(maxsplit=1)[1].endswith("BK")
        if is_bkb and (res in T[1] for T in term):
            continue
        hvy_lst.append(line)
    # update dict
    s1_info_d[pdb.stem]["MCCE.Step1"]["Missing Heavy Atoms"] = hvy_lst

    return s1_info_d


def info_s1_log(pdb: Path) -> dict:
    dout = {}
    s1log = RunLog1(pdb)
    # set the section data with dict:
    dout[pdb.stem] = {"MCCE.Step1": s1log.txt_blocks}
    dout = filter_heavy_atm_section(pdb, dout)

    return dout
