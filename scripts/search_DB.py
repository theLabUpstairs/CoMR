#!/usr/bin/env python
import subprocess
import sys
import os

input_file = snakemake.input[0]
mito_output = snakemake.output.mito_output
sub_output = snakemake.output.sub_output
subdb_path = snakemake.params.subdb_path
mitodb_path = snakemake.params.mitodb_path
threads = snakemake.threads
software_path = snakemake.params.software_path

def sub_mitoDB(input_file, subdb_path, mitodb_path, threads, sub_output, mito_output, software_path):
    # Command for subdb
    subdb_cmd = [
        software_path, "blastp", "-p", str(threads), "-d", subdb_path, "-q", input_file, 
        "-o", sub_output, "--outfmt", "6", "--more-sensitive", "-k", "500"
    ]
    # Running subdb_cmd
    subprocess.run(subdb_cmd, check=True)

    # Command for mitodb
    mitodb_cmd = [
        software_path, "blastp", "-p", str(threads), "-d", mitodb_path, "-q", input_file,
        "-o", mito_output, "--outfmt", "6", "--more-sensitive", "-k", "500"
    ]
    # Running mitodb_cmd
    subprocess.run(mitodb_cmd, check=True)

# Redirect stdout and stderr to the Snakemake log file
with open(snakemake.log[0], "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting searchDB", flush=True)
    sub_mitoDB(input_file, subdb_path, mitodb_path, threads, sub_output, mito_output, software_path)
    print("Search Sub/MitoDB done.", flush=True)
