# MCCE_ProtInfo
This tool provides information about a protein pdb file on its structure, along with applied transformations and issues flagged in the run.log produced by MCCE's step1.py.  
This information is saved into a file named `ProtInfo.md`.

  * See sample markdown [reports](#Samples)

## USAGE:
```
Command line interface for the MCCE ProtInfo tool, which gathers:
 * Info about the input protein from Bio.PDB parser.
 * Info from MCCE step1 run.log & debug.log when step1 can be run.

This is the 'main' module for the cli, which calls the function
that outputs a single protein report: `get_single_pdb_report(args)`.

Options:
 1. pdb (required): a pdb file name or pdbid (assumed valid).
 2. --fetch (False if not used): If 'pdb' is a pdbid and flag is used,
    the biological assembly is downloaded from rcsb.org.
 Step1 options:
  --wet (False): Keep water molecules.
  --noter (False): Do not label terminal residues (for making ftpl).
  -d (4.0): Protein dielectric constant for delphi.
  -u (''): User selected, comma-separated KEY=var pairs from run.prm; e.g.:
           -u HOME_MCCE=/path/to/mcce_home,EXTRA=./extra.tpl.
  -e (mcce): mcce executable location.
  -h, --help  Show this help message and exit.
  --fetch     Download the biological assembly of given pdb (if not a file).

Usage:
 >ProtInfo 1fat --fetch
 >ProtInfo 1fat.pdb --noter
```

## Basic Info:
Info recorded in the final report:
  * Number of models
  * Number of chains
  * Breaks in chains
  * Number of residues
  * Unknown residues
  * Residues with multiple alternate locations:
    - MCCE can only handle a single location and only considers the 'A' location (which may cause problems)
  * Number of water molecules
  * Number and identity of cofactors

## MCCE Info:
After running MCCE's step1.py, the run.log provides information about errors and repairs applied.  

Info recorded in the final report:
  * Terminal residues
  * Waters: count and buried list
  * Buried Cofactors:
    - Water: buried list & count
    - Other: buried list & count
    - Missing topology
  * __TODO__: Which repairs were applied  ::  Need more info.
  * Fatal errors

---
---

# Samples:

---
# 1fat :: Phytohemagglutinin-L
## ParsedStructure
### Chains:   <strong>4</strong>: A, B, C, D

### Residues:   956

### Waters:   16

### Buried (20% thresh.):
  - <strong>Waters</strong>: 15: A 305, A 306, A 307, A 308, B 301, B 302, B 303, C 313, C 314, C 315, C 316, D 309, D 310, D 311, D 312
  - <strong>Heteros</strong>: 8: A MN 254, A CA 255, B MN 254, B CA 255, C MN 254, C CA 255, D MN 254, D CA 255

### Warnings:
  <strong><font color='red'>Chain discontinuity</font></strong>: A (7975, 8039); B (7991, 8043); C (8007, 8047); D (8023, 8051); 

## MCCE.Step1
### Renamed:
  "CA    CA A 255" to "CA   _CA A 255"
  "CA    CA B 255" to "CA   _CA B 255"
  "CA    CA C 255" to "CA   _CA C 255"
  "CA    CA D 255" to "CA   _CA D 255"
  "MN    MN A 254" to "MN   _MN A 254"
  "MN    MN B 254" to "MN   _MN B 254"
  "MN    MN C 254" to "MN   _MN C 254"
  "MN    MN D 254" to "MN   _MN D 254"

### Termini:
  <strong>NTR</strong>: "SER A   1", "ASN A  38", "SER B   1", "SER C   1", "ASN C  38", "SER D   1", "ASN D  38"
  <strong>CTR</strong>: "ASN A  36", "SER A 233", "SER B 233", "ASN C  36", "SER C 233", "LEU D  35", "SER D 233"

### Labeling:
  - <strong><font color='red'>Generic topology file created for</font></strong>:
  NAG:  https://pubchem.ncbi.nlm.nih.gov/#query=NAG&tab=substance; _MN:  https://pubchem.ncbi.nlm.nih.gov/#query=MN&tab=substance; _CA:  https://pubchem.ncbi.nlm.nih.gov/#query=CA&tab=substance; 

### Missing Heavy Atoms:
  OXT of conf CTR01 in "CTR A  36".
  OXT of conf CTR01 in "CTR A 233".
  OXT of conf CTR01 in "CTR B 233".
  OXT of conf CTR01 in "CTR C  36".
  OXT of conf CTR01 in "CTR C 233".
  OXT of conf CTR01 in "CTR D  35".
  OXT of conf CTR01 in "CTR D 233".

### Distance Clashes:
<details><summary>Clashes found</summary>
  d= 1.53: " CA  NTR A   1" to " CB  SER A   1"
  d= 1.45: " ND2 ASN A  12" to " C1  NAG A 253"
  d= 1.53: " CA  NTR A  38" to " CB  ASN A  38"
  d= 1.52: " CA  NTR B   1" to " CB  SER B   1"
  d= 1.48: " ND2 ASN B  12" to " C1  NAG B 253"
  d= 1.53: " CA  NTR C   1" to " CB  SER C   1"
  d= 1.45: " ND2 ASN C  12" to " C1  NAG C 253"
  d= 1.52: " CA  NTR C  38" to " CB  ASN C  38"
  d= 1.87: " OD1 ASN C 128" to "CA   _CA C 255"
  d= 1.82: " NE2 HIS C 137" to "MN   _MN C 254"
  d= 1.54: " CA  NTR D   1" to " CB  SER D   1"
  d= 1.43: " ND2 ASN D  12" to " C1  NAG D 253"
  d= 1.55: " CA  NTR D  38" to " CB  ASN D  38"
</details>

---

---
# 1aig :: Photosynthetic Reaction Center From Rhodobacter Sphaeroides In The D+Qb-Charge Separated State
## ParsedStructure
### Chains:   <strong>3</strong>: L, M, H

### Residues:   890

### Waters:   53

### Buried (20% thresh.):
  - <strong>Waters</strong>: 16: L 289, L 290, L 292, L 295, L 297, L 298, L 300, L 301, L 304, L 305, M 313, M 314, M 315, M 322, M 327, H 272
  - <strong>Heteros</strong>: 1: M FE2 308

### Warnings:
  <strong><font color='red'>Chain discontinuity</font></strong>: L (6949, 7470); M (7274, 7489); H (7510); 

## MCCE.Step1
### Termini:
  <strong>NTR</strong>: "ALA L   1", "TYR M   3", "ASP H  11"
  <strong>CTR</strong>: "GLY L 281", "HIS M 301", "GLU H 258"

### Labeling:
  - <strong><font color='red'>Generic topology file created for</font></strong>:
  BPH:  https://pubchem.ncbi.nlm.nih.gov/#query=BPH&tab=substance; U10:  https://pubchem.ncbi.nlm.nih.gov/#query=U10&tab=substance; 
  - <strong><font color='red'>Unloadable topology</font></strong>:
  Atoms of residue BCL (L 282, L 283, M 309, M 310), do not match the topology conformer BCL-1.
 Likely cause: the renaming file (path: /home/cat/miniconda3/envs/rpte/name.txt) is missing entries for these species, resulting in unloadable topology files (path: /home/cat/miniconda3/envs/rpte/param/).

---
