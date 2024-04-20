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


def get_rcsb_compound_info(compound_id:str) -> str:
    """Query the RCSB api for compound information."""

    url_chemcomp = "https://data.rcsb.org/rest/v1/core/chemcomp/"
    url_smiles2pic = "https://cactus.nci.nih.gov/chemical/structure/"

    query = url_chemcomp + compound_id
    r = requests.get(query)
    if r.status_code != "200":
        print(f"Query unsuccessful for {compound_id!r}:\n{r.message}")
        return

    data = json.load(r.content)
    info_str = f"\tID: {data['chem_comp']['comp_id']}; Name: {data['chem_comp']['name']}\n"
    info_str = info_str + f"\tFormula: {data['chem_comp']['formula']}\n"
    info_str = info_str + f"Molecular Weight: {data['chem_comp']['formula_weight']}"

    types = ["InChI", "InChIKey", "SMILES", "SMILES_CANONICAL"]
    if data["pdbx_chem_comp_descriptor"]["type"] in types:
        encoding = data["pdbx_chem_comp_descriptor"]["descriptor"]
        info_str = info_str + f"\tLink to image: {url_smiles2pic}{encoding}/image\n"

    return info_str
