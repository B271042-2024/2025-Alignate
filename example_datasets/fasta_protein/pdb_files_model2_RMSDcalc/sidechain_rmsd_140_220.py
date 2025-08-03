from prody import *
import numpy as np
import sys

def compute_matched_sidechain_rmsd(pdb1_path, pdb2_path):
    ag1_full = parsePDB(pdb1_path)
    ag2_full = parsePDB(pdb2_path)

    # Filter only structured protein residues from 140+
    ag1 = ag1_full.select('protein and chain A and resnum 140 to 220')
    ag2 = ag2_full.select('protein and chain A and resnum 140 to 220')

    if ag1 is None or ag2 is None:
        print("❌ One of the selections returned nothing (check chain ID or residue range).")
        return

    # Align on CA atoms
    ca1 = ag1.select('name CA')
    ca2 = ag2.select('name CA')

    if ca1 is None or ca2 is None:
        print("❌ No CA atoms found in selected region for alignment.")
        return

    transform = calcTransformation(ca2, ca1)
    transform.apply(ag2)

    matched_atoms1 = []
    matched_atoms2 = []

    common_resnums = set(ag1.getResnums()).intersection(ag2.getResnums())
    common_resnums = [r for r in common_resnums if r >= 140]

    for resnum in sorted(common_resnums):
        res1 = ag1.select(f'resnum {resnum} and not name N CA C O H')
        res2 = ag2.select(f'resnum {resnum} and not name N CA C O H')
        if res1 is None or res2 is None:
            continue

        common = set(res1.getNames()).intersection(res2.getNames())
        for name in common:
            a1 = res1.select(f'name {name}')
            a2 = res2.select(f'name {name}')
            if a1 is not None and a2 is not None:
                matched_atoms1.append(a1)
                matched_atoms2.append(a2)

    if not matched_atoms1 or not matched_atoms2:
        print("❌ No matched side-chain atoms found.")
        return

    ag1_final = matched_atoms1[0]
    ag2_final = matched_atoms2[0]
    for a1, a2 in zip(matched_atoms1[1:], matched_atoms2[1:]):
        ag1_final += a1
        ag2_final += a2

    rmsd = calcRMSD(ag1_final, ag2_final)
    print(f"✅ Matched side-chain RMSD (residues 140+): {rmsd:.3f} Å")
    print(f"Atoms compared: {ag1_final.numAtoms()}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sidechain_rmsd_matched.py <pdb1> <pdb2>")
        sys.exit(1)

    pdb1 = sys.argv[1]
    pdb2 = sys.argv[2]
    compute_matched_sidechain_rmsd(pdb1, pdb2)
