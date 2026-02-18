#!/usr/bin/env python
import pandas as pd
import os, subprocess
from Bio import SeqIO
import sys

input_file_1 = snakemake.input.input1 #whole fasta
input_file_2 = snakemake.input.input2 #M start only fasta
output_file_1 = snakemake.output.output1 #parsed csv
output_file_2 = snakemake.output.output2 #targetp original output
software_path = snakemake.params.software_path
batch = snakemake.params.batch_path

def targetp(input_file_1, input_file_2, output_file_1, output_file_2):
    input_fasta_bn = os.path.basename(input_file_2).split(".fasta")[0]
    your_command = f"{software_path} -org non-pl -batch {batch} -fasta {input_file_2} > /dev/null"
    subprocess.run(your_command, shell=True)
    summary_file = f"{input_fasta_bn}_summary.targetp2"
    targetp_open = open(summary_file).readlines()
    # this runs targetp on metstart seq only

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

    # Process indices from targetp_open
    processed_indices = set()  # To keep track of indices processed
    for line in targetp_open:
        if not line.startswith("#"):
            parts = line.strip().split("\t")
            index = parts[0]
            processed_indices.add(index)  # Mark this index as processed and add it to processed_indices
            try:
                prediction, mts_prob, cleavage_pos = parts[1], parts[4], parts[5].split(".")[0]
            except IndexError:
                prediction, mts_prob, cleavage_pos = parts[1], parts[4], "No_cleavage_site_predicted"
            rows.append({
                'Index': index, 'tp_prediction': prediction, 'tp_score': mts_prob, 'tp_cleavage': cleavage_pos
            })

    # Process indices from non Met starting sequences
    for index in all_indices - processed_indices:
        rows.append({
            'Index': index, 'tp_prediction': "incomplete_seq", 'tp_score': "incomplete_seq", 
            'tp_cleavage': "incomplete_seq"
        })

    # Create DataFrame and write it to CSV 
    df = pd.DataFrame(rows)
    df.to_csv(output_file_1, index=False) # save parse targetp output
    os.rename(summary_file, output_file_2) # save original targetp output


# Redirect stdout and stderr to the Snakemake log file
with open(snakemake.log[0], "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting targetp", flush=True)
    targetp(input_file_1, input_file_2, output_file_1, output_file_2)
    print("Targetp done.", flush=True)
