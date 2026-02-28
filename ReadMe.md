### Similarity Search

----------
#### Overview
This module runs similarity search using three strategies: 
- Maximum Common Substructure (MCS)
- Fingerprint (ecfp4)
- Substructure matching

Given a **query list** and a **library**, the script identifies compounds in the library that are similar to each query compounds.

----------
#### Prerequisites
Two input files are required:
- `input_file_library`: the compound library to search against
- `input_file_query`: the list of query compounds

The script searches `input_file_library` for compounds similar to those in `input_file_query`.

----------

#### Known Issues

1) MCS Calculation in Rdkit (fused/bridged/spiro ring systems)

   The following RDKit MCS call may return an unintuitive or chemically inappropriate maximum common substructure for molecules containing fused rings, bridged bicyclic rings, or spiro bicyclic rings:
   ```python
   mcs_SMARTS = rdFMCS.FindMCS([mol_1, mol_2], ringMatchesRingOnly=True, completeRingsOnly=True)
   ```

   For the two molecules:
- O=C(NC1(O)C(=O)C2=CC=CC=C2C1=O)C1=CC=C(C(F)(F)F)C=C1 
- C12(NC(=O)c3ccc(cc3)C)CC3CC(C1)CC(C3)CC2

   RDKit may return the following MCS SMARTS:
- Cc1ccc(C(=O)NC2CCCCCCCC2)cc1

   This result is likely not the intended shared scaffold and may reflect limitations of the MCS settings for complex ring topologies.