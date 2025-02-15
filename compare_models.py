#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis import align
from Bio.PDB import PDBParser
import subprocess
import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# Define directories
hp_model_dir = "hp_model_structures"
af_model_dir = "af_model_structures"

# Proteins
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
        print(f"RMSD: {rmsd_value:.3f} Å")
        return rmsd_value
    except Exception as e:
        print(f"Error computing RMSD: {e}")
        return None

def compare_secondary_structure(pdb1, pdb2):
    """Compares the secondary structure of two proteins via DSSP using an external call."""
    try:
        dssp_executable = "mkdssp"
        dssp1_out = subprocess.run([dssp_executable, pdb1], capture_output=True, text=True)
        dssp2_out = subprocess.run([dssp_executable, pdb2], capture_output=True, text=True)

        sec_struct1 = [line.split()[2] for line in dssp1_out.stdout.splitlines() if len(line.split()) > 2]
        sec_struct2 = [line.split()[2] for line in dssp2_out.stdout.splitlines() if len(line.split()) > 2]

        if not sec_struct1 or not sec_struct2:
            print("No valid secondary structure data from DSSP.")
            return None
        
        # Similarity as percentage of matching residues
        similarity = sum(1 for a, b in zip(sec_struct1, sec_struct2) if a == b) / len(sec_struct1) * 100
        print(f"Secondary Structure Similarity: {similarity:.2f}%")
        return similarity
    except Exception as e:
        print(f"Error comparing secondary structure: {e}")
        return None
    
def compute_ca_distances(pdb_file):
    """Extracts consecutive CA-CA distances from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_coords.append(residue["CA"].get_coord())
    ca_coords = np.array(ca_coords)
    # Calculate distances between consecutive CA atoms
    distances = np.linalg.norm(np.diff(ca_coords, axis=0), axis=1)
    return distances

def compare_ca_distances(pdb1, pdb2):
    """Compares the average CA-CA bond distances between two PDB files."""
    dists1 = compute_ca_distances(pdb1)
    dists2 = compute_ca_distances(pdb2)
    # If lengths differ, compare up to the shorter length
    min_len = min(len(dists1), len(dists2))
    diff = np.abs(dists1[:min_len] - dists2[:min_len])
    print(f"Mean difference in consecutive CA distances: {np.mean(diff):.3f} Å")

def compute_ca_angles(pdb_file):
    """Extracts consecutive CA–CA–CA bond angles (in degrees) from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_coords.append(residue["CA"].get_coord())
    ca_coords = np.array(ca_coords)
    angles = []
    for i in range(1, len(ca_coords)-1):
        v1 = ca_coords[i] - ca_coords[i-1]
        v2 = ca_coords[i+1] - ca_coords[i]
        # Angle between v1 and v2
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Clip to avoid numerical issues
        angle = np.degrees(np.arccos(cos_angle))
        angles.append(angle)
    return np.array(angles)

def compare_ca_angles(pdb1, pdb2):
    """Compares the average difference in CA–CA–CA bond angles between two PDB files."""
    angles1 = compute_ca_angles(pdb1)
    angles2 = compute_ca_angles(pdb2)
    min_len = min(len(angles1), len(angles2))
    diff = np.abs(angles1[:min_len] - angles2[:min_len])
    print(f"Mean difference in CA bond angles: {np.mean(diff):.3f} degrees")

def compute_dihedral(p0, p1, p2, p3):
    """Computes the dihedral angle (in degrees) given four points."""
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def compute_radius_of_gyration(pdb_file):
    """Computes the radius of gyration for CA atoms in a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_coords.append(residue["CA"].get_coord())
    ca_coords = np.array(ca_coords)
    center = np.mean(ca_coords, axis=0)
    rg = np.sqrt(np.mean(np.sum((ca_coords - center)**2, axis=1)))
    return rg

def compare_radius_of_gyration(pdb1, pdb2):
    """Compares the radius of gyration between two PDB files."""
    rg1 = compute_radius_of_gyration(pdb1)
    rg2 = compute_radius_of_gyration(pdb2)
    print(f"Radius of Gyration Difference: {abs(rg1 - rg2):.3f} Å")

def compute_contact_map(pdb_file, threshold=8.0):
    """
    Computes a binary contact map for CA atoms in a PDB file.
    Two CA atoms are considered in contact if their distance is less than the threshold.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_coords.append(residue["CA"].get_coord())
    ca_coords = np.array(ca_coords)
    n = len(ca_coords)
    contact_map = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if np.linalg.norm(ca_coords[i] - ca_coords[j]) < threshold:
                contact_map[i, j] = 1
                contact_map[j, i] = 1
    return contact_map

def compare_contact_maps(pdb1, pdb2, threshold=8.0):
    """Compares the similarity of contact maps between two PDB files."""
    cm1 = compute_contact_map(pdb1, threshold)
    cm2 = compute_contact_map(pdb2, threshold)
    n1, n2 = cm1.shape[0], cm2.shape[0]
    n = min(n1, n2)
    diff = np.sum(np.abs(cm1[:n, :n] - cm2[:n, :n])) / (n*n)
    similarity = (1 - diff) * 100
    print(f"Contact map similarity (threshold {threshold} Å): {similarity:.2f}%")

def main():
    """Main function to compare protein models."""
    for protein_name, uniprot_id in proteins.items():
        # Select the HP model file (uncomment for GA version)
        hp_pdb = os.path.join(hp_model_dir, f"{protein_name}_HP.pdb")
        # hp_pdb = os.path.join(HP_MODEL_DIR, f"{protein_name}_HP_GA.pdb")
        af_pdb = os.path.join(af_model_dir, f"{protein_name}_AF.pdb")
        if os.path.exists(hp_pdb) and os.path.exists(af_pdb):
            print(f"\nComparing structures for {protein_name} ({uniprot_id})")
            compute_rmsd(hp_pdb, af_pdb)
            compare_secondary_structure(hp_pdb, af_pdb)
            compare_ca_distances(hp_pdb, af_pdb)
            compare_ca_angles(hp_pdb, af_pdb)
            compare_radius_of_gyration(hp_pdb, af_pdb)
            compare_contact_maps(hp_pdb, af_pdb, threshold=8.0)
        else:
            print(f"Missing PDB files for {protein_name}, skipping comparison.")

if __name__ == "__main__":
    main()