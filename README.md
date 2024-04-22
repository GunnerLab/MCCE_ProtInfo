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

## MCCE Recommendations: TODO
This section will list whatever MCCE can recommend.
---
---

# Sample Report:

---
# 1fat
## Input.ParsedStructure
### SingleModel.Chains
  - (4, ['A', 'B', 'C', 'D'])

### SingleModel.Residues
  - 956

### ParsedStructure.Warnings
  - Chain A is discontinuous at line 7975.
  - Chain B is discontinuous at line 7991.
  - Chain C is discontinuous at line 8007.
  - Chain D is discontinuous at line 8023.
  - Chain A is discontinuous at line 8039.
  - Chain B is discontinuous at line 8043.
  - Chain C is discontinuous at line 8047.
  - Chain D is discontinuous at line 8051.

## MCCE.Info
### Renamed:
  - "MN    MN A 254" to "MN   _MN A 254"
  - "CA    CA A 255" to "CA   _CA A 255"
  - "MN    MN B 254" to "MN   _MN B 254"
  - "CA    CA B 255" to "CA   _CA B 255"
  - "MN    MN C 254" to "MN   _MN C 254"
  - "CA    CA C 255" to "CA   _CA C 255"
  - "MN    MN D 254" to "MN   _MN D 254"
  - "CA    CA D 255" to "CA   _CA D 255"

### Terminii:
  - "SER A   1" as NTR
  - "ASN A  38" as NTR
  - "SER B   1" as NTR
  - "SER C   1" as NTR
  - "ASN C  38" as NTR
  - "SER D   1" as NTR
  - "ASN D  38" as NTR
  - "ASN A  36" as CTR
  - "SER A 233" as CTR
  - "SER B 233" as CTR
  - "ASN C  36" as CTR
  - "SER C 233" as CTR
  - "LEU D  35" as CTR
  - "SER D 233" as CTR

### Labeling:
  - <strong><font color='red'>Generic topology file created for: ['NAG', '_MN', '_CA']</font> </strong>

### Missing Heavy atoms:
  - N   of conf SERBK in "SER A   1".
  - CA  of conf SERBK in "SER A   1".
  - C   of conf ASNBK in "ASN A  36".
  - O   of conf ASNBK in "ASN A  36".
  - OXT of conf CTR01 in "CTR A  36".
  - N   of conf ASNBK in "ASN A  38".
  - CA  of conf ASNBK in "ASN A  38".
  - C   of conf SERBK in "SER A 233".
  - O   of conf SERBK in "SER A 233".
  - OXT of conf CTR01 in "CTR A 233".
  - N   of conf SERBK in "SER B   1".
  - CA  of conf SERBK in "SER B   1".
  - C   of conf SERBK in "SER B 233".
  - O   of conf SERBK in "SER B 233".
  - OXT of conf CTR01 in "CTR B 233".
  - N   of conf SERBK in "SER C   1".
  - CA  of conf SERBK in "SER C   1".
  - C   of conf ASNBK in "ASN C  36".
  - O   of conf ASNBK in "ASN C  36".
  - OXT of conf CTR01 in "CTR C  36".
  - N   of conf ASNBK in "ASN C  38".
  - CA  of conf ASNBK in "ASN C  38".
  - C   of conf SERBK in "SER C 233".
  - O   of conf SERBK in "SER C 233".
  - OXT of conf CTR01 in "CTR C 233".
  - N   of conf SERBK in "SER D   1".
  - CA  of conf SERBK in "SER D   1".
  - C   of conf LEUBK in "LEU D  35".
  - O   of conf LEUBK in "LEU D  35".
  - OXT of conf CTR01 in "CTR D  35".
  - N   of conf ASNBK in "ASN D  38".
  - CA  of conf ASNBK in "ASN D  38".
  - C   of conf SERBK in "SER D 233".
  - O   of conf SERBK in "SER D 233".
  - OXT of conf CTR01 in "CTR D 233".
  -    Ignore warning messages if they are in the terminal residues

### Distance Clashes:
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

### Connectivity:
  - Species and properties with assigned default values in debug.log:
  -   NAGBK: ['VDW_RAD', 'VDW_EPS']
  -   _MNBK: ['VDW_RAD', 'VDW_EPS']
  -   _CABK: ['VDW_RAD', 'VDW_EPS']

