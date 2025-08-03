from prody import *
import numpy as np
import sys

def compute_backbone_rmsd(pdb1_path, pdb2_path):
    ag1_full = parsePDB(pdb1_path)
    ag2_full = parsePDB(pdb2_path)

    # Select protein atoms from residue 140 onward
    ag1 = ag1_full.select('protein and chain A and resnum 140 to 999')
    ag2 = ag2_full.select('protein and chain A and resnum 140 to 999')

    if ag1 is None or ag2 is None:
        print("❌ Protein atoms not found in the specified range (resnum ≥140).")
        return

    # Select CA atoms for alignment
    ca1 = ag1.select('name CA')
    ca2 = ag2.select('name CA')

    if ca1 is None or ca2 is None:
        print("❌ No CA atoms found in selection for alignment.")
        return

    # Align structure 2 to structure 1
    transform = calcTransformation(ca2, ca1)
    transform.apply(ag2)

    # Now select backbone atoms (N, CA, C, O) only
    backbone1 = ag1.select('name N CA C O')
    backbone2 = ag2.select('name N CA C O')

    if backbone1 is None or backbone2 is None:
        print("❌ No backbone atoms (N, CA, C, O) found after residue 140.")
        return

    # Match atoms by residue number and atom name
    matched_atoms1 = []
    matched_atoms2 = []

    common_resnums = set(backbone1.getResnums()).intersection(backbone2.getResnums())

    for resnum in sorted(common_resnums):
        res1 = backbone1.select(f'resnum {resnum}')
        res2 = backbone2.select(f'resnum {resnum}')
        if res1 is None or res2 is None:
            continue

        common_atoms = set(res1.getNames()).intersection(res2.getNames())
        for atom_name in common_atoms:
            a1 = res1.select(f'name {atom_name}')
            a2 = res2.select(f'name {atom_name}')
            if a1 is not None and a2 is not None:
                matched_atoms1.append(a1)
                matched_atoms2.append(a2)

    if not matched_atoms1 or not matched_atoms2:
        print("❌ No matching backbone atoms found.")
        return

    ag1_final = matched_atoms1[0]
    ag2_final = matched_atoms2[0]
    for a1, a2 in zip(matched_atoms1[1:], matched_atoms2[1:]):
        ag1_final += a1
        ag2_final += a2

    # Final RMSD calculation
    rmsd = calcRMSD(ag1_final, ag2_final)
    print(f"✅ Matched backbone RMSD (residues 140+): {rmsd:.3f} Å")
    print(f"Atoms compared: {ag1_final.numAtoms()}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python backbone_rmsd_140plus.py <pdb1> <pdb2>")
        sys.exit(1)

    pdb1 = sys.argv[1]
    pdb2 = sys.argv[2]
    compute_backbone_rmsd(pdb1, pdb2)
