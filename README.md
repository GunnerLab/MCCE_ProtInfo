# MCCE_ProtInfo
This tool provides information about a protein pdb file on its structure, along with applied transformations and issues flagged in the run.log produced by MCCE's step1.py.  
This information is saved into a file named `ProtInfo`.[txt, md, toml :: TBD]

  * See sample markdown report for [1FAT](#Sample)

## Basic Info:
Info recorded in the final report:
  * Number of models
  * Number of chains
  * Breaks in chains
  * Number of residues
  * Unknown residues
  * Residues with alternate locations: MCCE can only handle single location
  * Number of water molecules
  * Number and identity of cofactors

## MCCE Info:
After running MCCE's step1.py, the run.log provides information about errors and repairs applied.  

Info recorded in the final report:
  * Terminal residues
  * Waters and their solvent-accessible surface area (SASA)
  * Cofactors:
    - HETATM id, name and, SASA
    - Missing topology
  * Which repairs were applied
  * Fatal errors

## TODO MCCE Recommendations:
This section will list whatever MCCE can recommend.

---
---

# Sample:

---
# 1fat
## ParsedStructure
### Chains
  * <strong>4 </strong> : A, B, C, D

### Residues
  - 956

### Waters
  - 4

### Warnings
  - <strong><font color='red'>Discontinuity</font> </strong>
    - Chain A: ['Line 7975', 'Line 8039']
    - Chain B: ['Line 7991', 'Line 8043']
    - Chain C: ['Line 8007', 'Line 8047']
    - Chain D: ['Line 8023', 'Line 8051']

## MCCE.Step1
### Renamed:
  - "CA    CA A 255" to "CA   _CA A 255"
  - "CA    CA B 255" to "CA   _CA B 255"
  - "CA    CA C 255" to "CA   _CA C 255"
  - "CA    CA D 255" to "CA   _CA D 255"
  - "MN    MN A 254" to "MN   _MN A 254"
  - "MN    MN B 254" to "MN   _MN B 254"
  - "MN    MN C 254" to "MN   _MN C 254"
  - "MN    MN D 254" to "MN   _MN D 254"

### Termini:
  * <strong>NTR </strong> : "SER A   1", "ASN A  38", "SER B   1", "SER C   1", "ASN C  38", "SER D   1", "ASN D  38"
  * <strong>CTR </strong> : "ASN A  36", "SER A 233", "SER B 233", "ASN C  36", "SER C 233", "LEU D  35", "SER D 233"

### Labeling:
  - <strong><font color='red'>Generic topology file created for:</font> </strong>
  - NAG::  https://pubchem.ncbi.nlm.nih.gov/#query=NAG&tab=substance
  - _MN::  https://pubchem.ncbi.nlm.nih.gov/#query=MN&tab=substance
  - _CA::  https://pubchem.ncbi.nlm.nih.gov/#query=CA&tab=substance

### Free Cofactors:
  - Species and properties with assigned default values in debug.log:
  - NAGBK: ['VDW_RAD', 'VDW_EPS']
  - _MNBK: ['VDW_RAD', 'VDW_EPS']
  - _CABK: ['VDW_RAD', 'VDW_EPS']

### Missing Heavy Atoms:
  - OXT of conf CTR01 in "CTR A  36".
  - OXT of conf CTR01 in "CTR A 233".
  - OXT of conf CTR01 in "CTR B 233".
  - OXT of conf CTR01 in "CTR C  36".
  - OXT of conf CTR01 in "CTR C 233".
  - OXT of conf CTR01 in "CTR D  35".
  - OXT of conf CTR01 in "CTR D 233".

### Distance Clashes:
<details><summary>Clashes found</summary>

  -    d= 1.53: " CA  NTR A   1" to " CB  SER A   1"
  -    d= 1.45: " ND2 ASN A  12" to " C1  NAG A 253"
  -    d= 1.53: " CA  NTR A  38" to " CB  ASN A  38"
  -    d= 1.52: " CA  NTR B   1" to " CB  SER B   1"
  -    d= 1.48: " ND2 ASN B  12" to " C1  NAG B 253"
  -    d= 1.53: " CA  NTR C   1" to " CB  SER C   1"
  -    d= 1.45: " ND2 ASN C  12" to " C1  NAG C 253"
  -    d= 1.52: " CA  NTR C  38" to " CB  ASN C  38"
  -    d= 1.87: " OD1 ASN C 128" to "CA   _CA C 255"
  -    d= 1.82: " NE2 HIS C 137" to "MN   _MN C 254"
  -    d= 1.54: " CA  NTR D   1" to " CB  SER D   1"
  -    d= 1.43: " ND2 ASN D  12" to " C1  NAG D 253"
  -    d= 1.55: " CA  NTR D  38" to " CB  ASN D  38"
  -    d= 1.70: "MN   _MN A 254" to " O   HOH A 307"
  -    d= 1.46: "MN   _MN A 254" to " O   HOH A 308"
  -    d= 1.99: "CA   _CA A 255" to " O   HOH A 306"
  -    d= 1.45: "MN   _MN B 254" to " O   HOH B 304"
  -    d= 1.69: "CA   _CA B 255" to " O   HOH B 301"
  -    d= 1.55: "MN   _MN C 254" to " O   HOH C 316"
  -    d= 1.50: "MN   _MN D 254" to " O   HOH D 311"
  -    d= 1.52: "MN   _MN D 254" to " O   HOH D 312"
</details>

---
