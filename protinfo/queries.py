#!/usr/bin/env python

import logging
from pathlib import Path
import requests
import sys
from typing import Union


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_rcsb_pdb(pdbid: str) -> Union[None, Path]:
    """Given a pdb id, download the pdb file containing
    the biological assembly from rcsb.org.
    The file is downloaded with a pdb extension.
    """

    url_rscb = "https://files.rcsb.org/download/"

    pdbid = pdbid.lower()
    bio_file = pdbid + ".pdb1"
    pdb_file = bio_file[:-1]
    try:
        r = requests.get(url_rscb + bio_file, allow_redirects=True)
        r.raise_for_status()
    except Exception:
        logger.exception(f"Error: Could not download {pdb_file}")
        sys.exit(1)

    with open(pdb_file, 'wb') as fo:
        fo.write(r.content)
    logger.info('Download completed.')

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
