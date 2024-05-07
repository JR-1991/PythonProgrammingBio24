"""
This script produces the dataset used in Exercise 004.
"""

import pandas as pd

from typing import List, Tuple
from utils import CODON_TABLE, to_triplets

def read_fasta_file(file: str) -> List[str]:
    """Read a FASTA file and return a list of sequences.

    Example:

        >> read_fasta_file("example.fasta")
        [('header1', 'ATCGAA'), ('header2', 'ATCGAA')]

    Args:
        file (str): Path to the FASTA file.

    Returns:
        List[str]: List of sequences.
    """

    with open(file) as f:
        lines = f.readlines()

    sequences = []
    for line in lines:
        if not line.startswith(">"):
            sequences.append(line.strip())

    return sequences

def calculate_gc(sequence: str) -> float:
    """Calculate the GC content of a nucleotide sequence.

    Example:

        >> calculate_gc("ATCGAA")
        0.5

    Args:
        sequence (str): Nucleotide sequence.

    Returns:
        float: GC content.
    """

    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence)

def get_codon_usage(sequence: str) -> dict[str, int]:
    """Calculate the codon usage of a nucleotide sequence."""

    triplets = to_triplets(sequence)
    n_triplets = len(triplets)
    codon_usage = {}

    for codon in CODON_TABLE:
        codon_usage[codon] = triplets.count(codon) / n_triplets

    return codon_usage

def analyse_sequences(sequences: List[str], label: str) -> pd.DataFrame:
    """Analyse a list of sequences and return a DataFrame with the GC content."""

    data = []
    for sequence in sequences:
        gc = calculate_gc(sequence)
        codon_usage = get_codon_usage(sequence)
        data.append({
            "organism": label,
            "lens": len(sequence),
            "gc": gc,
            **codon_usage,
        })

    df = pd.DataFrame(data)
    return df

if __name__ == "__main__":
    archaea = read_fasta_file("../data/archaea_sequences.fasta")
    coli = read_fasta_file("../data/ecoli_sequences.fasta")
    plants = read_fasta_file("../data/plant_sequences.fasta")

    archaea_df = analyse_sequences(archaea, "archaea")
    coli_df = analyse_sequences(coli, "coli")
    plants_df = analyse_sequences(plants, "plants")

    df = pd.concat([archaea_df, coli_df, plants_df])
    df.to_csv("../data/gc_len_data.csv", index=False)
