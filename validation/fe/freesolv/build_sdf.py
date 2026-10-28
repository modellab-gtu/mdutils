#!/usr/bin/env python
"""usage: build_sdf.py SMILES out.sdf  -- RDKit AddHs + ETKDG (seed 42) + MMFF; formal charge must be 0."""
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

smi, out = sys.argv[1:3]
m = Chem.AddHs(Chem.MolFromSmiles(smi))
if sum(a.GetFormalCharge() for a in m.GetAtoms()) != 0:
    sys.exit(f"{smi}: formal charge is not 0")
if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
    sys.exit(f"{smi}: ETKDG embedding failed")
if m.GetNumConformers() == 0:
    sys.exit(f"{smi}: no conformer")
if AllChem.MMFFHasAllMoleculeParams(m):
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
else:
    sys.exit(f"{smi}: no MMFF parameters")
w = Chem.SDWriter(out)
w.write(m)
w.close()
