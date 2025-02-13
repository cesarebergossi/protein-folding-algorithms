#!/usr/bin/env python3
import requests

def get_amino_acid_sequence(uniprot_id):
    """
    Retrieves the amino acid sequence (in FASTA format) for a given UniProt ID.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    # Check if the request was successful
    if response.status_code == 200:
        lines = response.text.splitlines()
        # Exclude the header and join the sequence lines
        sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
        return sequence
    else:
        raise ValueError(f"Error fetching sequence for {uniprot_id}. Status: {response.status_code}")

# Example usage  
if __name__ == "__main__":
    print(get_amino_acid_sequence("P0CG47"))