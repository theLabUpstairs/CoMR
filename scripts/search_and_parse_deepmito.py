#!/usr/bin/env python
import math
import os
import shlex
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

# Snakemake variables
input_file_1 = snakemake.input.input1  # whole fasta
input_file_2 = snakemake.input.input2  # M start only fasta
output_file_1 = snakemake.output.output1  # parsed CSV
output_file_2 = snakemake.output.output2  # DeepMito summary
software_path = snakemake.params.software_path  # Docker command
db_path = snakemake.params.db_path  # Database path
log_file = snakemake.log[0]
SUPPRESS_PATTERNS = {
    "tensorflow_warning": "WARNING:tensorflow",
    "syntax_warning": "SyntaxWarning:",
}

def _normalize_command(cmd):
    """Return a list command regardless of how it was provided."""
    if isinstance(cmd, (list, tuple)):
        return list(cmd)
    if isinstance(cmd, str):
        return shlex.split(cmd)
    raise ValueError(f"Unsupported command specification: {cmd!r}")


def run_deepmito(input_file, output_file, command, db_path):
    """Run DeepMito CLI with the packaged environment."""
    abs_input = os.path.abspath(input_file)
    abs_output = os.path.abspath(output_file)
    abs_db = os.path.abspath(db_path)

    cmd = command + [
        "multi-fasta",
        "-f",
        abs_input,
        "-d",
        abs_db,
        "-o",
        abs_output,
    ]
    print(f"Running: {' '.join(cmd)}", flush=True)
    env = os.environ.copy()
    env.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
    suppressed_counts = {key: 0 for key in SUPPRESS_PATTERNS}
    seen_patterns = {key: False for key in SUPPRESS_PATTERNS}
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    ) as process:
        assert process.stdout is not None
        for raw_line in process.stdout:
            line = raw_line.rstrip()
            matched_key = next(
                (
                    key
                    for key, pattern in SUPPRESS_PATTERNS.items()
                    if pattern in line
                ),
                None,
            )
            if matched_key:
                if seen_patterns[matched_key]:
                    suppressed_counts[matched_key] += 1
                    continue
                seen_patterns[matched_key] = True
            print(line, flush=True)
        retcode = process.wait()
    if retcode != 0:
        raise subprocess.CalledProcessError(retcode, cmd)

    for key, count in suppressed_counts.items():
        if count:
            print(
                f"Suppressed {count} additional log lines containing '{SUPPRESS_PATTERNS[key]}'",
                flush=True,
            )


def count_fasta_sequences(fasta_path):
    count = 0
    with open(fasta_path, "r") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def split_fasta(fasta_path, chunk_size, temp_dir):
    """Split FASTA file into chunks with up to chunk_size sequences."""
    chunk_paths = []
    chunk_index = 0
    current_handle = None
    seqs_in_chunk = 0

    def open_new_chunk():
        nonlocal chunk_index, current_handle, seqs_in_chunk
        chunk_path = os.path.join(temp_dir, f"chunk_{chunk_index:04d}.fasta")
        chunk_index += 1
        chunk_paths.append(chunk_path)
        current_handle = open(chunk_path, "w")
        seqs_in_chunk = 0

    current_header = None
    current_seq = []

    with open(fasta_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    if current_handle is None:
                        open_new_chunk()
                    current_handle.write(current_header + "\n")
                    current_handle.write("\n".join(current_seq) + "\n")
                    seqs_in_chunk += 1
                    if seqs_in_chunk >= chunk_size:
                        current_handle.close()
                        current_handle = None
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        if current_handle is None:
            open_new_chunk()
        current_handle.write(current_header + "\n")
        current_handle.write("\n".join(current_seq) + "\n")
        seqs_in_chunk += 1

    if current_handle is not None:
        current_handle.close()

    return chunk_paths


def merge_summaries(chunk_summaries, final_output):
    """Concatenate DeepMito summaries keeping header only once."""
    if not chunk_summaries:
        open(final_output, "w").close()
        return

    first_file = True
    with open(final_output, "w") as out_handle:
        for summary in chunk_summaries:
            with open(summary, "r") as in_handle:
                for line in in_handle:
                    if first_file or not line.startswith("##"):
                        out_handle.write(line)
            first_file = False


def run_deepmito_in_parallel(input_file, output_file, command, db_path, workers):
    """Split a FASTA and run DeepMito chunks concurrently."""
    total_sequences = count_fasta_sequences(input_file)
    if total_sequences == 0:
        print("Input FASTA has no sequences; creating empty outputs.", flush=True)
        open(output_file, "w").close()
        return

    if workers <= 1 or total_sequences == 1:
        run_deepmito(input_file, output_file, command, db_path)
        return

    chunk_size = max(1, math.ceil(total_sequences / workers))
    print(
        f"Splitting {total_sequences} sequences into chunks of {chunk_size} "
        f"with up to {workers} workers.",
        flush=True,
    )

    with tempfile.TemporaryDirectory(prefix="deepmito_chunks_") as temp_dir:
        chunk_paths = split_fasta(input_file, chunk_size, temp_dir)
        chunk_outputs = [f"{chunk}.summary" for chunk in chunk_paths]
        pool_workers = min(workers, len(chunk_paths))
        print(
            f"Running {len(chunk_paths)} DeepMito chunks using {pool_workers} workers.",
            flush=True,
        )

        with ThreadPoolExecutor(max_workers=pool_workers) as executor:
            futures = [
                executor.submit(
                    run_deepmito,
                    chunk_path,
                    chunk_output,
                    command,
                    db_path,
                )
                for chunk_path, chunk_output in zip(chunk_paths, chunk_outputs)
            ]
            for future in as_completed(futures):
                future.result()

        merge_summaries(chunk_outputs, output_file)



def parse_deepmito(summary_file, output_csv):
    """Parse the DeepMito summary output and save as CSV."""
    rows = []
    with open(summary_file, 'r') as file:
        for line in file:
            if not line.startswith("##"):
                results = line.strip().split("\t")
                rows.append({
                    'Sequence': results[0],
                    'dm_location': results[2],
                    'dm_score': results[5],
                    'dm_GO_term': results[8]
                })
    
    # Save results to CSV
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)

# Redirect output to log file
with open(log_file, "w") as log:
    # Redirect Python-level stdout and stderr
    sys.stdout = log
    sys.stderr = log

    # Redirect system-level stdout and stderr
    os.dup2(log.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log.fileno(), 2)  # Redirect fd=2 (stderr)

    # Debug message to confirm redirection
    print("Logging initialized. Starting DeepMito process", flush=True)

    threads = max(1, getattr(snakemake, "threads", 1))
    command = _normalize_command(software_path)

    print(f"Starting DeepMito with {threads} worker(s)", flush=True)
    run_deepmito_in_parallel(input_file_2, output_file_2, command, db_path, threads)
    parse_deepmito(output_file_2, output_file_1)
    print("DeepMito done.", flush=True)
