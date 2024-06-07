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
from protinfo.io_utils import get_path_keys, ENV
import logging
from time import sleep
from typing import Union


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


MSG_KEEP_WATS = "NOTE: Add the '--wet' option at the command line to keep waters and cofactors."


def get_pubchem_compound_link(compound_id: str) -> str:
    """Return the link of the PubChem subtance tab for compound_id.
    The link is prepended with ":  " as it will follow compound_id
    in the report line; this saves checking for empty str.
    """
    if compound_id:
        url_fstr = ":  https://pubchem.ncbi.nlm.nih.gov/#query={}&tab=substance"
        if compound_id.startswith("_"):
            return url_fstr.format(compound_id[1:])
        else:
            return url_fstr.format(compound_id)
    else:
        return ""


def extract_content_between_tags(text: str, tag1: str, tag2: str = "   Done") -> Union[str, None]:
    """Extracts the content between two string tags in a text.
    Args:
      text: A text.
      tag1: The first tag.
      tag2: The second tag, default: "   Done".

    Returns:
      A string containing the content between the two tags if both found;
      A string containing the content from tag1 if tag2 not found;
      None if tag1 is not found.
    """

    pos1 = text.find(tag1)
    if pos1 == -1:
        return None

    start_pos = pos1 + len(tag1)
    end_pos = text.find(tag2, start_pos)
    if end_pos != -1:
        return text[start_pos:end_pos]
    else:
        return text[start_pos:]


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
    get_key: str = None

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


def get_log1_specs(pdb: Path) -> dict:
    """Return a dict of LogHdr classes for processing step1 sections
    in run.log.
    """

    # list of run.log headers returned by mcce step1:
    runlog1_headers = [
        "   Rename residue and atom names...",
        "   Identify NTR and CTR...",
        "   Label backbone, sidechain and altLoc conformers...",
        "   Load pdb lines into data structure...",
        # 5: dynamic header
        "   Strip free cofactors with SAS >  {: .0%}...",
        "   Check missing heavy atoms and complete altLoc conformers...",
        # 7: dynamic header
        "   Find distance clash (<{:.3f})...",
        "   Make connectivity network ...",
    ]

    all = defaultdict(dict)
    for i, hdr in enumerate(runlog1_headers, start=1):
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
            all[i] = LogHdr(
                i,
                hdr,
                rpt_hdr="Labeling",
                line_start="      Labeling ",
                skip_lines=(
                    "Creating temporary parameter file for unrecognized",
                    "Trying labeling again",
                    "Try delete this entry and run MCCE again",
                    "Error! premcce_confname()",
                    "STOP",
                    "is already loaded somewhere else.",
                ),
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
                get_key="H2O_SASCUTOFF",
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
            all[i] = LogHdr(i, hdr, rpt_hdr="Distance Clashes", get_key="CLASH_DISTANCE")
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


class RunLog1:
    """A class to parse mcce run.log into sections pertaining to step1, and
    process each one of them into a simplified output.
    """

    def __init__(self, pdb: Path) -> None:
        self.pdb = pdb.resolve()
        self.pdbid = self.pdb.stem
        self.s1_dir = self.pdb.parent.joinpath("step1_run")
        self.runprm = self.get_runprm()
        self.dry_opt = float(self.runprm["H2O_SASCUTOFF"]) == -0.01
        # id of block with debug.log mentions, if any:
        self.check_debuglog_idx = [5]
        self.blocks_specs = get_log1_specs(pdb)
        self.txt_blocks = self.get_blocks()

    def get_runprm(self) -> dict:
        env = ENV(self.s1_dir)

        return env.runprm

    def get_debuglog_species(self) -> str:
        fp = self.s1_dir.joinpath("debug.log")
        df = pd.read_csv(fp, sep=r"\s+", header=None, engine="python")
        txt = "Species and properties with assigned default values in debug.log:\n"
        for k in df[1].unique():
            txt = txt + f"{k}: {list(df[df[1] == k][0].unique())}\n"

        return txt

    def process_content_block(self, content: list, loghdr: LogHdr) -> list:
        out = []
        skip = loghdr.skip_lines is not None
        change = loghdr.line_start is not None
        newtpl = None
        tpl_mismatch = None
        tpl_err = "   Error! The following atoms of residue "

        # section that lists new.tpl creation:
        if loghdr.idx == 3:
            newtpl = ""

        for line in content:
            if not line:
                continue

            if loghdr.idx == 3:
                if line.startswith("   Error! premcce_confname()"):
                    # add conf name & link:
                    conf = line.rsplit(maxsplit=1)[1]
                    newtpl += f"{conf}{get_pubchem_compound_link(conf)}; "

                elif line.startswith(tpl_err):
                    if tpl_mismatch is None:
                        tpl_mismatch = defaultdict(list)

                    res_info, tpl_conf = line.removeprefix(tpl_err).split(" can not be loaded to conformer type ")
                    res, resloc = res_info.split(maxsplit=1)
                    tpl_mismatch[(res, tpl_conf.strip())].append(resloc)
                    continue
                elif line.startswith("           ") and tpl_mismatch is not None:
                    continue

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

        if loghdr.idx == 5:
            if self.dry_opt:
                out.insert(0, MSG_KEEP_WATS)

        # check if new tpl confs:
        if loghdr.idx == 3:
            if newtpl:
                out.append("Generic topology file created for")
                out.append(newtpl)

            if tpl_mismatch:
                out.append("Unloadable topology")
                msg = ""
                locs = ""
                fmt = "Atoms of residue {} ({}), do not match the topology conformer {}.\n"
                for k in tpl_mismatch:
                    locs = ", ".join(lx for lx in tpl_mismatch[k])
                    msg = msg + fmt.format(k[0], locs, k[1])

                d = get_path_keys(self.pdb)
                msg = msg + " Likely cause: the renaming file (path: " + d["renaming file"]
                msg = msg + ") is missing entries for these species, resulting in unloadable"
                msg = msg + " topology files (path: " + d["topologies"] + "/)."
                out.append(msg)

        return out

    def get_blocks(self) -> dict:
        """Extract 'processing blocks' from run.log file."""

        log_fp = self.s1_dir.joinpath("run.log")
        text = log_fp.read_text()

        block_txt = {}
        for k in self.blocks_specs:
            lhdr = self.blocks_specs[k]
            rpt_k = lhdr.rpt_hdr
            if k in [5, 7]:
                # dynamic headers
                h = lhdr.hdr
                lhdr.hdr = h.format(float(self.runprm[lhdr.get_key]))
            content = extract_content_between_tags(text, lhdr.hdr)
            if content is None:
                continue
            else:
                content = content.splitlines()

            if k == 1:
                content = sorted(content)

            if (lhdr.line_start is not None) or (lhdr.skip_lines is not None):
                content = self.process_content_block(content, lhdr)

            block_txt[rpt_k] = [line for line in content if line.strip()]

        # process termini; group res into NTR, CTR
        b2_hdr = self.blocks_specs[2].rpt_hdr
        if block_txt.get(b2_hdr) is not None:
            if block_txt[b2_hdr]:
                termi = defaultdict(list)
                for line in block_txt[b2_hdr]:
                    i = line.index('"', 3) + 1
                    termi[line[-3:]].append(line[:i])
                block_txt[b2_hdr] = []
                for k in termi:
                    block_txt[b2_hdr].append((k, termi[k]))

        if self.blocks_specs[5].debuglog:
            # add extra line for each species found:
            rk = self.blocks_specs[5].rpt_hdr
            for line in self.get_debuglog_species().splitlines():
                block_txt[rk].append(line)

        # collapse dist clashes block 7:
        if block_txt.get("Distance Clashes") is not None:
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

    # termi values are [2-tuples]
    termi = s1_info_d[pdb.stem]["MCCE.Step1"].get("Termini")
    heavy = s1_info_d[pdb.stem]["MCCE.Step1"].get("Missing Heavy Atoms")
    if heavy is None:
        return s1_info_d

    if termi is None:
        return s1_info_d

    if len(heavy) > 1:
        _ = heavy.pop(-1)
        # == Ignore warning messages if they are in the terminal residues

    hvy_lst = []
    for line in heavy:
        conf, res = line.split(" in ")
        is_bkb = conf.rsplit(maxsplit=1)[1].endswith("BK")
        if is_bkb and (res in T[1] for T in termi):
            continue
        hvy_lst.append(line)

    # update dict
    s1_info_d[pdb.stem]["MCCE.Step1"]["Missing Heavy Atoms"] = hvy_lst

    return s1_info_d


def info_s1_log(pdb: Path) -> dict:
    dout = {}
    sleep(5)

    s1log = RunLog1(pdb)
    # set the section data with dict:
    dout[pdb.stem] = {"MCCE.Step1": s1log.txt_blocks}
    dout = filter_heavy_atm_section(pdb, dout)

    return dout
