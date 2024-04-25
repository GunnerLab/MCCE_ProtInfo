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
from typing import Union


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

    def has_debuglog(self, line:str, calling_key:int=5) -> bool:
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
    "   Strip free cofactors with SAS >   5%...",   # 5
    "   Check missing heavy atoms and complete altLoc conformers...",
    "   Find distance clash (<2.000)...",
    "   Make connectivity network ...",
    ]


def get_log1_specs(loghdrs:list) -> dict:
    """Return a dict of LogHdr classes for processing step1 sections
    in run.log.
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
                            rpt_hdr="Termini:",
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
blocks_specs = get_log1_specs(runlog1_headers)


class RunLog1:
    """A class to parse mcce run.log into sections pertaining to step1, and
    process each one of them into a simplified output.
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

        # section that lists new.tpl creation:
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

                if line.startswith("   Total deleted cofactors"):
                    n_cof = int(line.rsplit(maxsplit=1)[1][:-1])
                    if n_cof == 0:
                        continue

            if skip:
                if isinstance(lhdr.skip_lines, tuple):
                    found = False
                    for t in lhdr.skip_lines:
                        found = found or (t in line)
                    if found:
                        continue
                else:
                    if line in lhdr.skip_lines:
                        continue

            if change:
                # remove common start:
                if line.startswith(lhdr.line_start):
                    line = line.removeprefix(lhdr.line_start)

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
            rpt_k = lhdr.rpt_hdr
            content = extract_content_between_tags(text, lhdr.hdr).splitlines()
            if k == 1:
                content = sorted(content)

            if (lhdr.line_start is not None) or (lhdr.skip_lines is not None):
                content = self.process_content_block(content, lhdr)
            block_txt[rpt_k] = [line for line in content if line.strip()]

        # process termini: '"SER A 1" as NTR'
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
            # add extra line:
            rk = blocks_specs[5].rpt_hdr
            for line in self.get_debuglog_species().splitlines():
                block_txt[rk].append(line)

        return block_txt


def info_s1_log(pdb:Path) -> dict:
    dout = {}
    silog = RunLog1(pdb)
    #match struc in prot info dict:
    dout[pdb.stem] = {"MCCE.Info": silog.txt_blocks}

    return dout
