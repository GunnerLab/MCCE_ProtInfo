#!/usr/bin/env python

import json
from pathlib import Path
import requests
from typing import Union


def get_rcsb_pdb(pdbid:str) -> Union[None, Path]:
    """Given a pdb id, download the pdb file containing
    the biological assembly from rcsb.org.
    The file is downloaded with a pdb extension.
    """

    url_rscb = "https://files.rcsb.org/download/"

    pdbid = pdbid.lower()
    bio_file = pdbid + ".pdb1"
    pdb_file = bio_file[:-1]
    r = requests.get(url_rscb + bio_file, allow_redirects=True)
    if r.status_code != "200":
        print(f"Could not download {pdb_file}:\n{r.message}")
        return None
    
    with open(pdb_file, 'wb') as fo:
        fo.write(r.content)
    print('Download completed.')

    return Path(pdb_file).resolve()


def get_pubchem_compound_link(compound_id:str) -> str:
    """Return the unvalidated link of the PubChem page for compund_id. """

    if compound_id:
        url_fstr = "https://pubchem.ncbi.nlm.nih.gov/#query={}&tab=compound"
        return url_fstr.format(compound_id.upper())
    else:
        return ""
