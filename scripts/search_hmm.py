#!/usr/bin/env python
import subprocess
import os 

input_file = snakemake.input[0]
output_file = snakemake.output[0]
software_path = snakemake.params.software_path
db_path = snakemake.params.db_path
threads = snakemake.params.threads

def HMMscan(input_file, db_path, threads, software_path, output_file):
    hmmscan_cmd = [
        software_path,
        "--cpu", str(threads),
        "--tblout", output_file,
        db_path,
        input_file
    ]
    subprocess.run(hmmscan_cmd, check=True)

# Redirect stdout and stderr to the Snakemake log file
with open(snakemake.log[0], "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting hmmscan", flush=True)
    HMMscan(input_file, db_path, threads, software_path, output_file)
    print("Hmmscan done.", flush=True)