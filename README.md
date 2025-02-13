# Protein Folding as a Computational Mechanism: Comparing AlphaFold and the HP Model

This repository is a Python framework for exploring and comparing protein folding models. The project implements coarse-grained HP (hydrophobic–polar) models using both greedy and genetic algorithms, and compares these simplified folds with high-resolution AlphaFold predictions.

## Features

- **HP Model (Greedy Algorithm):**  
  Generates a 3D lattice-based folding model from an amino acid sequence by greedily optimizing hydrophobic contacts.

- **HP Model (Genetic Algorithm):**  
  Evolves folding pathways using a genetic algorithm that optimizes non-covalent hydrophobic contacts.

- **AlphaFold Integration:**  
  Downloads protein structures from the AlphaFold Protein Structure Database using Biopython.

- **Structural Comparison:**  
  Calculates RMSD and secondary structure similarity between the coarse-grained HP models and AlphaFold predictions.

- **Visualization:**  
  Creates 3D plots of the folding paths for easy visual comparison.

## Installation

### Prerequisites

- Python 3.7 or later
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [Biopython](https://biopython.org/)
- [MDAnalysis](https://www.mdanalysis.org/)
- [Requests](https://docs.python-requests.org/)

> **Note:**  
> For secondary structure analysis, ensure that DSSP (e.g. `mkdssp`) is installed and available in your PATH.  
> On macOS, you can install DSSP via Homebrew:
> ```bash
> brew install dssp
> ```
> On Ubuntu/Debian:
> ```bash
> sudo apt-get update && sudo apt-get install dssp
> ```
> Alternatively, install via Conda:
> ```bash
> conda install -c salilab dssp
> ```

### Setup

1. **Clone the repository:**

   ```bash
   git clone https://github.com/yourusername/FoldCompare.git
   cd FoldCompare
   ```

2. **Install dependencies:**

   ```bash
   pip install numpy matplotlib biopython MDAnalysis requests
   ```

## Usage
This repository contains several scripts. Here’s a brief overview:

- **HP Model (Greedy Algorithm):**
Generate an HP model using the greedy algorithm and produce a PDB file and 3D visualization.
 
	```bash
  	python hp_model.py
  	```

- **HP Model (Genetic Algorithm):**
Generate an HP model using a genetic algorithm.

	```bash
  	python hp_model_genetic.py
  	```

- **AlphaFold Integration:**
Download and visualize AlphaFold predicted structures.

	```bash
	python alphafold.py
	```

- **Structural Comparison:**
Compare the HP model with AlphaFold structures (computing RMSD and secondary structure similarity).

	```bash
	python compare_models.py
 	```

*Note on Scale*:
Since the HP models are highly simplified and lattice-based, their spatial scales do not match the detailed AlphaFold predictions. When comparing plots, it is best to present them side by side and discuss differences in topology and contact patterns rather than absolute distances.

## Project Structure

```
FoldCompare/
├── af_model_structures/      (AlphaFold PDB files generated)
├── af_model_plots/           (AlphaFold structure visualizations)
├── comparison_plots/         (Optional comparison plots between models)
├── hp_model_structures/      (HP model PDB files generated)
├── hp_model_plots/           (HP model folding visualizations)
├── hp_model.py               (Greedy algorithm for HP model folding)
├── hp_model_genetic.py       (Genetic algorithm for HP model folding)
├── alphafold.py              (AlphaFold model retrieval and visualization)
├── compare_models.py         (Structural comparison tools)
├── uniprot_to_seq.py         (Utility to fetch protein sequences from UniProt)
└── README.md
```

## Contributing

Contributions are welcome! Please follow these steps:
	1.	Fork the repository.
	2.	Create a new branch (e.g., “git checkout -b feature/YourFeature”).
	3.	Commit your changes.
	4.	Push to your branch.
	5.	Open a pull request.

## Acknowledgements
- AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/
- Biopython: https://biopython.org/
- MDAnalysis: https://www.mdanalysis.org/
- UniProt: https://www.uniprot.org/
