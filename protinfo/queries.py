#!/usr/bin/env python

import logging
from pathlib import Path
import protinfo.io_utils as iou
import requests
from typing import Tuple, Union


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def rcsb_download(pdb_fname: str) -> requests.Response:
    url_rscb = "https://files.rcsb.org/download/"
    return requests.get(url_rscb + pdb_fname, allow_redirects=True)


def get_rcsb_pdb(pdbid: str) -> Union[Path, Tuple[None, str]]:
    """Given a pdb id, download the pdb file containing
    the biological assembly from rcsb.org.
    The file is downloaded with a pdb extension.
    """

    pdbid = pdbid.lower()
    bionames = pdbid + ".pdb1", f"{pdbid}-assembly1.cif.gz"
    pdb_file = pdbid + ".pdb"

    content = None
    # list of bool to identify which bio assembly was saved:
    which_ba = [False, False]  # 0: pdb, 1: cif

    # try bio assemblies first:
    r0 = rcsb_download(bionames[0])
    if r0.status_code < 400:
        which_ba[0] = True
        pdb_file = bionames[0][:-1]
        content = r0.content
    else:
        logger.error(f"Error: Could not download the pdb bio assembly:{r0.reason}")

        r1 = rcsb_download(bionames[1])
        if r1.status_code < 400:
            which_ba[1] = True
            pdb_file = bionames[1]
            content = r1.content
        else:
            logger.error(f"Error: Could not download the cif bio assembly:{r1.reason}")

    if which_ba[0] == which_ba[1]:  # both False; last try: legacy pdb
        r2 = rcsb_download(pdb_file)
        if r2.status_code < 400:
            content = r2.content
        else:
            logger.error(f"Error: Could not download the pdb file:{r2.reason}")
            return None, "Error: Could neither download the bio assembly or pdb file."

    # save file:
    with open(pdb_file, "wb") as fo:
        fo.write(content)
    logger.info("Download completed.")

    if which_ba[1]:
        decomp = iou.decompress_gz(Path(pdb_file))
        logger.info(f"{pdb_file} saved & unzipped as {decomp.name}")
        pdb_file = iou.cif2pdb(decomp)

    return Path(pdb_file).resolve()


def get_pubchem_compound_link(compound_id: str) -> str:
    """Return the unvalidated link of the PubChem page for compund_id.
    subtance tab
    """

    if compound_id:
        url_fstr = "https://pubchem.ncbi.nlm.nih.gov/#query={}&tab=substance"
        return url_fstr.format(compound_id.upper())
    else:
        return ""
