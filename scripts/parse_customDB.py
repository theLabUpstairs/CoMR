#!/usr/bin/env python

import pandas as pd
from Bio import SearchIO
from Bio import BiopythonDeprecationWarning
import warnings
import sys
import os

# File paths from Snakemake
custom_output = snakemake.input[0]  # Input CustomDB BLAST output
output_file = snakemake.output[0]  # Output parsed CSV file
log_file = snakemake.log[0]  # Log file

def parse_custom(custom_output):
    """
    Parses the CustomDB BLAST output and extracts the top hit information for each sequence.

    Parameters:
        custom_output (str): Path to the CustomDB BLAST output file.

    Returns:
        pd.DataFrame: A DataFrame containing the sequence ID, top hit, bitscore, and e-value.
    """
    results = SearchIO.parse(custom_output, "blast-tab")
    rows = []

    for record in results:
        seq_id = str(record.id)
        top_hit_acc = "None"
        top_hit_bitscore = "None"
        top_hit_evalue = "None"

        # Extract top hit if available
        if len(record.hits) > 0:
            top_hit_acc = record.hits[0].id
            top_hit_bitscore = record.hsps[0].bitscore
            top_hit_evalue = record.hsps[0].evalue

        rows.append({
            "Sequence ID": seq_id,
            "CustomDB_header": top_hit_acc,
            "CustomDB_hit_bitscore": top_hit_bitscore,
            "CustomDB_hit_evalue": top_hit_evalue
        })

    return pd.DataFrame(rows)

def main():
    # Parse the CustomDB output
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
        df = parse_custom(custom_output)

    # Replace NaN or missing values with "None" for consistency
    df.fillna("None", inplace=True)

    # Save the parsed DataFrame to CSV
    df.to_csv(output_file, index=False)

# Redirect stdout and stderr to the Snakemake log file
with open(log_file, "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting parse_customDB", flush=True)
    main()
    print("CustomDB parsing done.", flush=True)
