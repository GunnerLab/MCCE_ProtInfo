# MCCE_ProtInfo
Provide information about a protein pdb file on its structure, along with applied transformations and issues flagged in the run.log produced by MCCE's step1.py.  
This information is saved into a file named `ProtInfo`.[txt, yaml, toml :: TBD]

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

## MCCE Recommendations:
This section will list whatever MCCE can recommend.
