#!/usr/bin/env python
import os
import subprocess
import pandas as pd
from Bio import SeqIO
import sys

# Snakemake inputs and outputs
input_file = snakemake.input.input2  # Multi-FASTA input
output_csv = snakemake.output.output1  # Parsed CSV
output_summary = snakemake.output.output2  # MitoProt summary
log_file = snakemake.log[0]
software_path = snakemake.params.software_path  # Path to wrapper script

def split_fasta(input_file, output_dir):
    """Split a multi-FASTA file into individual FASTA files."""
    os.makedirs(output_dir, exist_ok=True)
    for record in SeqIO.parse(input_file, "fasta"):
        output_file = os.path.join(output_dir, f"{record.id}.fasta")
        with open(output_file, "w") as out_f:
            SeqIO.write(record, out_f, "fasta")
    return [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith(".fasta")]

def run_mitoprot(input_file, software_path):
    """Run MitoProt on a single sequence file."""
    output_dir = os.path.dirname(input_file)
    cmd = f"{software_path} {input_file} {output_dir}"  # Removed '-f a'
    subprocess.run(cmd, shell=True, check=True)

def parse_mitoprot_results(fasta_file, output_dir, parsed_results):
    """Parse MitoProt results for a single sequence."""
    # Construct the expected MitoProt output file name
    mitoprot_output = os.path.join(output_dir, f"{os.path.basename(fasta_file).split('.')[0]}.mitoprot")

    # Check if the output file exists
    if not os.path.isfile(mitoprot_output):
        raise FileNotFoundError(f"MitoProt output file not found: {mitoprot_output}")

    # Append the raw output to the summary file
    with open(mitoprot_output, "r") as file:
        lines = file.readlines()

    # Parse relevant lines
    cleavagesite, mitoprot_score, mp_mts_seq = 0, "n.a.", ""
    for line in lines:
        if line.startswith("CleavSite   :     "):
            cleavagesite = line.replace("CleavSite   :     ", "").strip()
        if line.startswith("DFM         :          "):
            line = line.replace("DFM         :          ", "").strip()
            mitoprot_score = int(float(line.split()[1]) * 100)
            break

    # Extract MTS sequence
    if int(cleavagesite) > 0:
        for record in SeqIO.parse(fasta_file, "fasta"):
            mp_mts_seq = str(record.seq)[:int(cleavagesite)].strip()
            break

    # Append results
    parsed_results.append({
        "SequenceID": os.path.basename(fasta_file).split(".")[0],
        "mp_cleavagesite": cleavagesite,
        "mp_MTS_sequence": mp_mts_seq,
        "mp_score": mitoprot_score
    })

    return lines  # Return the raw output for concatenation


def main():
    # Redirect output to log file
    with open(log_file, "w") as log:
        sys.stdout = log
        sys.stderr = log
        os.dup2(log.fileno(), 1)
        os.dup2(log.fileno(), 2)

        # Step 1: Split the input multi-FASTA file
        temp_dir = os.path.join(os.path.dirname(output_csv), "temp_sequences")
        fasta_files = split_fasta(input_file, temp_dir)
        print(f"Running MitoProt on {len(fasta_files)} sequence(s).", flush=True)

        # Step 2: Run MitoProt on each sequence and concatenate outputs
        print("Starting Mitoprot")
        parsed_results = []
        with open(output_summary, "w") as summary_handle:
            for fasta_file in fasta_files:
                run_mitoprot(fasta_file, software_path)
                # Append raw MitoProt output to the summary
                raw_output = parse_mitoprot_results(fasta_file, os.path.dirname(fasta_file), parsed_results)
                summary_handle.writelines(raw_output)
                summary_handle.write("\n")

        # Step 3: Aggregate results into a single CSV
        df = pd.DataFrame(parsed_results)
        df.to_csv(output_csv, index=False)

        # Cleanup temporary directory (optional)
        for file in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir, file))
        os.rmdir(temp_dir)
        print("Mitoprot done.")

if __name__ == "__main__":
    main()
