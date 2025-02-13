#!/usr/bin/env python3
import numpy as np
import random
import os
import matplotlib.pyplot as plt
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom
from uniprot_to_seq import get_amino_acid_sequence

# Set random seed for reproducibility
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

class HPGeneticFolding:
    """
    A 3D HP folding model using a genetic algorithm.
    
    The algorithm evolves a population of folds (represented as self-avoiding walks on a lattice)
    to minimize energy calculated from non-covalent hydrophobic contacts. 
    It applies elitism, crossover and mutation to generate new folds.
    
    Attributes:
        protein_name (str): Name of the protein.
        hp_sequence (str): Protein sequence in HP notation.
        best_fold (list): The best fold (list of 3D coordinates) found.
        best_energy (float): Energy of the best fold.
    """
    def __init__(self, protein_name, hp_sequence, population_size=100, generations=200, mutation_rate=0.1):
        self.protein_name = protein_name
        self.hp_sequence = hp_sequence
        self.population_size = population_size
        self.generations = generations
        self.mutation_rate = mutation_rate
        self.directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]  # 3D moves
        self.best_fold = None
        self.best_energy = float('inf')

        self.run_genetic_algorithm()
        self.save_to_pdb()
        self.visualize_3d_structure()

    def generate_random_fold(self):
        """Generates a random self-avoiding walk on the lattice representing a fold."""
        position = (0, 0, 0)
        path = [position]
        for i in range(1, len(self.hp_sequence)):
            # Find available moves that don't lead to an already occupied position
            available = [(position[0] + dx, position[1] + dy, position[2] + dz)
                         for dx, dy, dz in self.directions 
                         if (position[0] + dx, position[1] + dy, position[2] + dz) not in path]
            if not available:
                break  # Terminate if no available moves
            position = random.choice(available)
            path.append(position)
        return path

    def fitness(self, path):
        """
        Calculates the fitness of a fold based on non-covalent hydrophobic contacts.
        
        Each non-covalent contact between hydrophobic residues (excluding adjacent residues)
        contributes a penalty of -6. The fitness value is divided by 2 to avoid double counting.
        """
        score = 0
        for i, pos in enumerate(path):
            if self.hp_sequence[i] != 'H':
                continue
            for dx, dy, dz in self.directions:
                neighbor = (pos[0] + dx, pos[1] + dy, pos[2] + dz)
                if neighbor in path:
                    j = path.index(neighbor)
                    # Exclude covalent (adjacent) contacts
                    if abs(i - j) > 1:
                        score -= 6
        return score / 2

    def crossover(self, parent1, parent2):
        """Combines two parent folds to create offspring."""
        cut = random.randint(1, len(self.hp_sequence) - 1)
        child = parent1[:cut] + parent2[cut:]
        return child

    def mutate(self, path):
        """Mutates a fold by shifting one random position (if available)."""
        if random.random() < self.mutation_rate and len(path) > 1:
            idx = random.randint(1, len(path)-1)
            dx, dy, dz = random.choice(self.directions)
            new_pos = (path[idx][0] + dx, path[idx][1] + dy, path[idx][2] + dz)
            if new_pos not in path:
                path[idx] = new_pos

    def run_genetic_algorithm(self):
        """Runs the genetic algorithm to evolve a population of folds."""
        population = [self.generate_random_fold() for _ in range(self.population_size)]
        for gen in range(self.generations):
            population = sorted(population, key=self.fitness)
            new_population = population[:10]  # Elitism: keep top 10 folds
            while len(new_population) < self.population_size:
                try:
                    parent1, parent2 = random.sample(population[:50], 2)
                except ValueError:
                    parent1, parent2 = population[0], population[1]
                child = self.crossover(parent1, parent2)
                self.mutate(child)
                new_population.append(child)
            population = new_population
        self.best_fold = population[0]
        self.best_energy = self.fitness(self.best_fold)
        print(f"Best Fold Found for {self.protein_name} with Energy {self.best_energy}")

    def save_to_pdb(self):
        """
        Saves the best fold as a PDB file.

        Backbone atoms (N, CA, C, O) are added using approximate bond distances:
            N-CA ~1.45 Å, CA-C ~1.52 Å, C-O ~1.23 Å.
        """
        structure = Structure.Structure("HP")
        model = Model.Model(0)
        chain = Chain.Chain("A")
        for i, (x, y, z) in enumerate(self.best_fold, start=1):
            residue = Residue.Residue((" ", i, " "), "GLY", " ")
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
        filename = os.path.join(pdb_output_dir, f"{self.protein_name}_HP_GA.pdb")
        io.save(filename)
        print(f"HP GA model structure saved as {filename}")

    def visualize_3d_structure(self):
        """
        Generates and saves a 3D plot of the HP folding path.
        
        Residues are colored based on hydrophobicity: red for hydrophobic ('H') and blue for polar ('P').
        """
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        x_coords, y_coords, z_coords = zip(*self.best_fold)
        colors = ["red" if self.hp_sequence[i] == 'H' else "blue" for i in range(len(self.best_fold))]
        ax.scatter(x_coords, y_coords, z_coords, c=colors, s=100)
        ax.plot(x_coords, y_coords, z_coords, color="black", linestyle="-", linewidth=1)
        ax.set_title(f"3D HP Model Folding (GA): {self.protein_name}")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.set_zlabel("Z-axis")
        output_path = os.path.join(plot_output_dir, f"{self.protein_name}_HP_GA.png")
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"3D visualization saved as {output_path}")

if __name__ == "__main__":
    # Process proteins using the genetic algorithm
    for name, uniprot_id in proteins.items():
        try:
            protein_sequence = get_amino_acid_sequence(uniprot_id)
            hp_sequence = convert_to_hp(protein_sequence)
            hp_model = HPGeneticFolding(name, hp_sequence)
        except Exception as e:
            print(f"Error processing {name}: {e}")