#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom
from uniprot_to_seq import get_amino_acid_sequence
import random
import os

# Random seed for reproducibility
random.seed(42)

# Hydrophobic residues: Glycine (G), Alanine (A), Valine (V), Leucine (L), Isoleucine (I), Proline (P), Phenylalanine (F), Methionine (M), Tryptophan (W)
hydrophobic_residues = {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W'}

# Define proteins with UniProt IDs
proteins = {
    "Ubiquitin": "P0CG47",
    "Lysozyme": "P00698"
}

# Output directories
pdb_output_dir = "hp_model_structures"
plot_output_dir = "hp_model_plots"
os.makedirs(pdb_output_dir, exist_ok=True)
os.makedirs(plot_output_dir, exist_ok=True)

def convert_to_hp(sequence):
    """Converts an amino acid sequence into hydrophobic-polar (HP) notation."""
    sequence = sequence.upper()
    return ''.join('H' if aa in hydrophobic_residues else 'P' for aa in sequence)

class HPModel3D:
    """
    A 3D HP folding model using a greedy algorithm.
    
    This class generates a coarse-grained lattice-based fold from an HP sequence. 
    It maximizes nonbonded hydrophobic contacts (each contact reduces the energy). 
    The resulting fold is saved as a PDB file and visualized in 3D.
    
    Attributes:
        protein_name (str): Name of the protein.
        hp_sequence (str): Protein sequence in HP notation.
        lattice (dict): Maps lattice coordinates (tuple) to (residue type, sequence index).
        energy (int): Total energy (lower is better).
        path (list): List of lattice coordinates representing the fold.
        pdb_filename (str): Path where the PDB file is saved.
    """
    def __init__(self, protein_name, hp_sequence):
        self.protein_name = protein_name
        self.hp_sequence = hp_sequence
        self.lattice = {}
        self.energy = 0 
        self.directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
        self.pdb_filename = os.path.join(pdb_output_dir, f"{protein_name}_HP.pdb")
        
        self.fold_protein()
        self.visualize_3d_structure()

    def fold_protein(self):
        """
        Folds the HP sequence onto a 3D lattice using a greedy algorithm.
        
        Starting at the origin, at each step the algorithm chooses the available lattice position that
        maximizes the number of non-covalent hydrophobic contacts (ignoring immediate neighbors).
        """
        # Initialize the fold at the origin
        start = (0, 0, 0)
        self.lattice[start] = (self.hp_sequence[0], 0)
        path = [start]
        current_position = start

        for i in range(1, len(self.hp_sequence)):
            best_pos = None
            best_contacts = -1
            
            # Shuffle available directions to randomize tie-breakers
            random.shuffle(self.directions)
            for dx, dy, dz in self.directions:
                new_pos = (current_position[0] + dx,
                           current_position[1] + dy,
                           current_position[2] + dz)
                if new_pos in self.lattice:
                    continue  # Skip if already occupied

                # Count non-adjacent hydrophobic contacts
                contacts = 0
                for ddx, ddy, ddz in self.directions:
                    neighbor = (new_pos[0] + ddx, new_pos[1] + ddy, new_pos[2] + ddz)
                    if neighbor in self.lattice:
                        res_type, idx = self.lattice[neighbor]
                        # Exclude immediate neighbors (covalent bonds)
                        if res_type == 'H' and abs(idx - i) > 1:
                            contacts += 1

                if contacts > best_contacts:
                    best_contacts = contacts
                    best_pos = new_pos

            if best_pos is None:
                print(f"Chain growth blocked at residue {i} for {self.protein_name}. Terminating fold.")
                break
            
            # Place residue at the best position
            self.lattice[best_pos] = (self.hp_sequence[i], i)
            current_position = best_pos
            path.append(best_pos)
            self.energy -= best_contacts  # Lower energy (more contacts) is favorable

        self.path = path
        print(f"{self.protein_name}: Folded with energy {self.energy} and {len(self.path)} residues.")
        self.save_to_pdb()

    def save_to_pdb(self):
        """
        Saves the HP fold as a PDB file.

        Backbone atoms (N, CA, C, O) are added using approximate bond distances:
            N-CA ~1.45 Å, CA-C ~1.52 Å, C-O ~1.23 Å.
        """
        structure = Structure.Structure("HP")
        model = Model.Model(0)
        chain = Chain.Chain("A")

        for i, (x, y, z) in enumerate(self.path, start=1):
            residue = Residue.Residue((" ", i, " "), "GLY", " ")
            # Approximate positions along z-axis; these are placeholders.
            atom_positions = {
                "N":  np.array([x, y - 0.3, z - 1.45], dtype=float),
                "CA": np.array([x, y, z], dtype=float),
                "C":  np.array([x, y + 0.3, z + 1.52], dtype=float),
                "O":  np.array([x + 0.5, y + 0.3, z + 1.52 + 1.23], dtype=float)
            }
            for atom_name, pos in atom_positions.items():
                element = atom_name[0]
                atom = Atom.Atom(atom_name, pos, 1.0, 1.0, " ", atom_name, i, element)
                residue.add(atom)
            chain.add(residue)

        model.add(chain)
        structure.add(model)
        io = PDBIO()
        io.set_structure(structure)
        io.save(self.pdb_filename)
        print(f"HP model structure saved as {self.pdb_filename}")

    def visualize_3d_structure(self):
        """
        Generates and saves a 3D plot of the HP folding path.
        
        Residues are colored based on hydrophobicity: red for hydrophobic ('H') and blue for polar ('P').
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        x_coords, y_coords, z_coords = zip(*self.path)
        colors = ["red" if self.hp_sequence[i] == 'H' else "blue" for i in range(len(self.path))]
        ax.scatter(x_coords, y_coords, z_coords, c=colors, s=100)
        ax.plot(x_coords, y_coords, z_coords, color="black", linestyle="-", linewidth=1)
        ax.set_title(f"3D HP Model Folding: {self.protein_name}\nEnergy: {self.energy}")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.set_zlabel("Z-axis")
        output_path = os.path.join(plot_output_dir, f"{self.protein_name}_HP.png")
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"3D visualization saved as {output_path}")

if __name__ == "__main__":
    # Process proteins
    for name, uniprot_id in proteins.items():
        try:
            sequence = get_amino_acid_sequence(uniprot_id)
            hp_sequence = convert_to_hp(sequence)
            hp_model = HPModel3D(name, hp_sequence)
        except Exception as e:
            print(f"Error processing {name}: {e}")