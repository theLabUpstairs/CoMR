
#!/usr/bin/env python
import pandas as pd
import argparse
import os, subprocess, re
from Bio import SearchIO, SeqIO
import random, string
import sys
import time

input_file_1 = snakemake.input.input1 #whole fasta
input_file_2 = snakemake.input.input2 #M start only fasta
output_file_1 = snakemake.output.output1 #parsed csv
output_file_2 = snakemake.output.output2 #targetp original output
software_path = snakemake.params.software_path


def mitofates(input_file_1, input_file_2, output_file_1, output_file_2):
    input_fasta_bn = os.path.basename(input_file_2).split(".fasta")[0]
    your_command = f"{software_path} {input_file_2} fungi > {input_fasta_bn}.MFout"
    subprocess.run(your_command, shell=True)
    summary_file = f"{input_fasta_bn}.MFout"
    mitofates_open = open(summary_file).readlines()
    # this runs mitofates on metstart seq only

    # Parse the indexed fasta file to get all sequence identifiers
    all_indices = set()
    with open(input_file_1, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):  # FASTA headers start with ">"
                header = line[1:].strip()  # Remove the ">" and any trailing newline characters
                index = header.split()[0]  # Assuming the first part of the header is the identifier
                all_indices.add(index)

    # Initialize rows list (converted in panda dataframe later)
    rows = []

    # Process indices from mitofates_open
    processed_indices = set()  # To keep track of indices processed
    for line in mitofates_open:
        if not line.startswith("Sequence"):
            parts = line.strip().split("\t")
            index = parts[0]
            processed_indices.add(index)  # Mark this index as processed and add it to processed_indices
            mf_score, mf_seq, mf_cleavage_pos = parts[1].replace(" ","_"), parts[2].replace(" ","_"), parts[3].replace(" ","_")
            rows.append({
                'Index': index, 'mf_score': mf_score, 'mf_seq': mf_seq, 'mf_cleavage_pos': mf_cleavage_pos
            })
    
    # Process indices from non Met starting sequences
    for index in all_indices - processed_indices:
        rows.append({
            'Index': index, 'mf_score': "incomplete_seq", 'mf_seq': "incomplete_seq", 
            'mf_cleavage_pos': "incomplete_seq"
        })

    # Create DataFrame and write it to CSV 
    df = pd.DataFrame(rows)
    df.to_csv(output_file_1, index=False) # save parse targetp output
    os.rename(summary_file, output_file_2) # save original targetp output


# Redirect low-level file descriptors (Perl)

with open(snakemake.log[0], "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    # Debug message to confirm redirection
    print("Starting mitofates", flush=True)
    # Call the mitofates function
    mitofates(input_file_1, input_file_2, output_file_1, output_file_2)
    print("Mitofates done.", flush=True)

