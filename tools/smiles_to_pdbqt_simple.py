#!/usr/bin/env python3
"""Generate PDBQT from SMILES using RDKit"""

import sys
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_pdbqt(smiles: str, output_path: str, name: str = "MOL"):
    """Convert SMILES to PDBQT via PDB"""
    
    # Generate 3D structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"ERROR: Invalid SMILES: {smiles}")
        sys.exit(1)
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    
    # Save as PDB
    pdb_path = output_path.replace('.pdbqt', '.pdb')
    writer = Chem.PDBWriter(pdb_path)
    writer.write(mol)
    writer.close()
    
    # Convert PDB to PDBQT using obabel
    cmd = [
        "obabel",
        pdb_path,
        "-O", output_path,
        "-p", "7.4",  # pH
        "--partialcharge", "gasteiger"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: obabel conversion failed")
        print(result.stderr)
        sys.exit(1)
    
    # Cleanup PDB
    Path(pdb_path).unlink()
    
    print(f"✅ Generated: {output_path}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--name", default="MOL")
    
    args = parser.parse_args()
    
    smiles_to_pdbqt(args.smiles, args.output, args.name)
