#!/usr/bin/env python
import subprocess
import sys
import os

input_file = snakemake.input[0]
output_file = snakemake.output[0]
customdb_path = snakemake.params.customdb_path
threads = snakemake.threads
software_path = snakemake.params.software_path

def customDB(input_file, customdb_path, threads, output_file, software_path):
    # Command for subdb
    customdb_cmd = [
        software_path, "blastp", "-p", str(threads), "-d", customdb_path, "-q", input_file, 
        "-o", output_file, "--outfmt", "6", "--more-sensitive", "-k", "500"
    ]
    # Running subdb_cmd
    subprocess.run(customdb_cmd, check=True)

# Call the function
# Redirect stdout and stderr to the Snakemake log file
with open(snakemake.log[0], "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting searchCustomDB", flush=True)
    customDB(input_file, customdb_path, threads, output_file, software_path)
    print("Search_CustomDB done.", flush=True)
