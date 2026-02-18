#!/usr/bin/env python

import pandas as pd
from Bio import SearchIO
import sys
import os
from Bio import SeqIO

# File paths from Snakemake
input_hmm_file = snakemake.input.hmm  # Input HMMscan output file
input_fasta_file = snakemake.input.fasta  # Input FASTA file with all sequence IDs
output_csv_file = snakemake.output.csv  # Output parsed CSV file
log_file = snakemake.log[0]  # Log file

def parse_hmmscan(input_hmm_file, input_fasta_file, output_csv_file):
    """
    Parses an HMMscan output file in hmmer3-tab format, ensuring all sequences are represented.
    
    Parameters:
        input_hmm_file (str): Path to the input HMMscan output file.
        input_fasta_file (str): Path to the input FASTA file containing all sequence IDs.
        output_csv_file (str): Path to the output CSV file.
    """
    # Load all sequence IDs from the input FASTA file
    all_sequences = {record.id: {"HMMscan_hit_name": "No HMMscan hits", "HMMscan_hit_eval": "No HMMscan hits"} 
                     for record in SeqIO.parse(input_fasta_file, "fasta")}

    # Parse the HMMscan output
    with open(input_hmm_file, "r") as HMM_output:
        HMM_results = SearchIO.parse(HMM_output, "hmmer3-tab")
        for record in HMM_results:
            seqID = str(record.id)
            if len(record.hits) > 0:
                # Use the best hit (first hit) for each sequence
                top_hit = record.hits[0]
                all_sequences[seqID] = {
                    "HMMscan_hit_name": top_hit.id,
                    "HMMscan_hit_eval": top_hit.hsps[0].evalue,
                }

    # Convert to DataFrame and save as CSV
    df = pd.DataFrame.from_dict(all_sequences, orient="index")
    df.index.name = "Sequence ID"
    df.to_csv(output_csv_file)
    print(f"Parsed results saved to {output_csv_file}")

# Redirect stdout and stderr to the log file
with open(log_file, "w") as log_file:
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting parse_hmmscan", flush=True)
    parse_hmmscan(input_hmm_file, input_fasta_file, output_csv_file)
    print("Hmmscan parsing done.", flush=True)
