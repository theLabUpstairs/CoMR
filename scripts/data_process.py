#!/usr/bin/env python
import pandas as pd
import argparse
import os, subprocess, re
from Bio import SearchIO, SeqIO
import random, string
import sys
import time

# Redirect stdout and stderr to the Snakemake log file
with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = log_file
    print("Starting fasta parsing", flush=True)
    
    input_file = snakemake.input[0]  # If there's only one input file
    output_file = snakemake.output[0]  # If there's only one output file

    fasta_index = SeqIO.index(input_file, "fasta")
    #input_fasta_bn = os.path.basename(input_file).split(".fasta")[0]

    fasta_keys=list(fasta_index.keys())
    df=pd.DataFrame(fasta_keys,columns=['seqID'])
    df=df.set_index("seqID")
    metlist=[]
    for record in SeqIO.parse(input_file, "fasta"):
        seqID = str(record.id)
        seq = str(record.seq)
        if seq.startswith("M"):
            df.at[seqID,"MetStart"] = "M"
            metlist.append(record)
        else:
            df.at[seqID,"MetStart"] = "P"
    #output_file="%s_metstart.fasta" % (input_fasta_bn)

    # Instead of SeqIO.write, manually write the records to ensure single-line sequences
    with open(output_file, 'w') as output_handle:
        for record in metlist:
            # Manually format the FASTA record for single-line sequence output
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")
    
    print("Fasta parsing done.", flush=True)


