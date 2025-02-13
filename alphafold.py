#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB
from Bio.PDB.alphafold_db import get_structural_models_for
import os

# Define proteins (Uniprot IDs)
proteins = {
    "Ubiquitin": "P0CG47",
    "Lysozyme": "P00698"
}

# Output directories
pdb_output_dir = "af_model_structures"
plot_output_dir = "af_model_plots"
os.makedirs(pdb_output_dir, exist_ok=True)
os.makedirs(plot_output_dir, exist_ok=True)

def download_alphafold_pdb(protein_name, uniprot_id):
    """Retrieves the AlphaFold structure and saves the first model as a PDB file."""
    structures = list(get_structural_models_for(uniprot_id, directory=pdb_output_dir))
    if not structures:
        print(f"No structures found for {protein_name} ({uniprot_id}).")
        return None
    structure = structures[0]
    io = PDB.PDBIO()
    io.set_structure(structure)
    pdb_filename = os.path.join(pdb_output_dir, f"{protein_name}_AF.pdb")
    io.save(pdb_filename)
    print(f"AlphaFold structure saved: {pdb_filename}")
    return pdb_filename

def extract_backbone_coordinates(pdb_file):
    """Extracts CA (alpha carbon) atom coordinates from a PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_coords.append(residue["CA"].get_coord())
    return np.array(ca_coords)

def plot_protein_structure(coords, protein_name):
    """Plots the CA backbone of the protein in 3D and saves the plot as an image."""
    output_path = os.path.join(plot_output_dir, f"{protein_name}_AF.png")
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(coords[:,0], coords[:,1], coords[:,2], marker="o", linestyle="-", color="blue")
    ax.set_title(f"AlphaFold Structure: {protein_name}")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved: {output_path}")

if __name__ == "__main__":
    # Process proteins: download structure, extract CA coordinates and plot
    for name, uniprot_id in proteins.items():
        pdb_path = download_alphafold_pdb(name, uniprot_id)
        if pdb_path and os.path.exists(pdb_path):
            backbone_atoms = extract_backbone_coordinates(pdb_path)
            print(f"Extracted {len(backbone_atoms)} residues for {name}.")
            if backbone_atoms.size > 0:
                plot_protein_structure(backbone_atoms, name)
            else:
                print(f"No CA atoms found for {name}.")
        else:
            print(f"Skipping visualization for {name} due to missing file.")