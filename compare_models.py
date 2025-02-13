#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis import align
import subprocess
import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# Define directories
hp_model_dir = "hp_model_structures"
af_model_dir = "af_model_structures"

# Proteins (same as before)
proteins = {
    "Ubiquitin": "P0CG47",
    "Lysozyme": "P00698"
}

def compute_rmsd(pdb1, pdb2):
    """
    Computes RMSD (root mean square deviation) between two structures using CA atoms.
    Aligns pdb2 (second structure) onto pdb1 (first structure) before calculation.
    """
    try:
        u1 = mda.Universe(pdb1)
        u2 = mda.Universe(pdb2)
        align.AlignTraj(u2, u1, select="name CA", in_memory=True).run()
        rmsd_calc = RMSD(u1, u2, select="name CA")
        rmsd_calc.run()
        rmsd_value = rmsd_calc.results.rmsd[-1, 2]
        print(f"RMSD between {pdb1} and {pdb2}: {rmsd_value:.3f} Å")
        return rmsd_value
    except Exception as e:
        print(f"Error computing RMSD for {pdb1} and {pdb2}: {e}")
        return None

def compare_secondary_structure(pdb1, pdb2):
    """
    Compares the secondary structures of two proteins via DSSP using an external call.
    Note: Ensure that DSSP (mkdssp) is installed and on your PATH.
    """
    try:
        dssp_executable = "mkdssp"
        dssp1_out = subprocess.run([dssp_executable, pdb1], capture_output=True, text=True)
        dssp2_out = subprocess.run([dssp_executable, pdb2], capture_output=True, text=True)
        dssp1_lines = dssp1_out.stdout.splitlines()
        dssp2_lines = dssp2_out.stdout.splitlines()
        # Extract secondary structure information (third column)
        sec_struct1 = [line.split()[2] for line in dssp1_lines if len(line.split()) > 2]
        sec_struct2 = [line.split()[2] for line in dssp2_lines if len(line.split()) > 2]
        if len(sec_struct1) == 0:
            print("No secondary structure data from DSSP for first structure.")
            return None
        # Compute similarity as percentage of matching residues
        similarity = sum(1 for a, b in zip(sec_struct1, sec_struct2) if a == b) / len(sec_struct1) * 100
        print(f"Secondary Structure Similarity: {similarity:.2f}%")
        return similarity
    except Exception as e:
        print(f"Error comparing secondary structure: {e}")
        return None

if __name__ == "__main__":
    # Compare each protein
    for name, uniprot_id in proteins.items():
        hp_pdb = os.path.join(hp_model_dir, f"{name}_HP.pdb")
        # Uncomment the next line to compare the GA version instead:
        # hp_pdb = os.path.join(hp_model_dir, f"{name}_HP_GA.pdb")
        af_pdb = os.path.join(af_model_dir, f"{name}_AF.pdb")
        if os.path.exists(hp_pdb) and os.path.exists(af_pdb):
            print(f"\nComparing structures for {name} ({uniprot_id})")
            compute_rmsd(hp_pdb, af_pdb)
            compare_secondary_structure(hp_pdb, af_pdb)
        else:
            print(f"Missing PDB files for {name}, skipping comparison.")