import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# Snakemake inputs and outputs
fasta_file = snakemake.input.fasta
csv_file = snakemake.input.csv
align_path = snakemake.params.align_path
alignments_dir = snakemake.output.alignments_dir
log_file = snakemake.log[0]
global_threads = max(
    1, getattr(snakemake, "threads", snakemake.params.threads)
)
mafft_threads = max(
    1, min(global_threads, snakemake.params.align_threads)
)
min_trimmed_residues = max(
    0, getattr(snakemake.params, "min_trimmed_residues", 1)
)

GAP_CHARACTERS = {"-", "."}

# Parse the CSV file
csv_data = pd.read_csv(csv_file)
sequence_hits = csv_data[csv_data["HMMscan_hit_name"] != "No HMMscan hits"]

# Parse the multi-FASTA file
sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Ensure output directory exists
os.makedirs(alignments_dir, exist_ok=True)


def ungapped_length(sequence) -> int:
    """Return the count of non-gap characters in a sequence."""
    return sum(1 for char in str(sequence) if char not in GAP_CHARACTERS)


def clean_trimal_output(trimal_file: str, target_id: str, min_residues: int):
    """
    Remove gap-only sequences from the trimmed alignment. If the target sequence
    is removed, delete the entire alignment to prevent tree building.
    """
    records = list(SeqIO.parse(trimal_file, "fasta"))
    trimal_path = Path(trimal_file)

    if not records:
        if trimal_path.exists():
            trimal_path.unlink()
        return [], False

    kept_records = []
    removed_ids = []
    target_removed = False

    for record in records:
        non_gap_len = ungapped_length(record.seq)
        if non_gap_len >= min_residues:
            kept_records.append(record)
        else:
            removed_ids.append(record.id)
            if record.id == target_id:
                target_removed = True

    if not kept_records or target_removed:
        if trimal_path.exists():
            trimal_path.unlink()
        return removed_ids, False

    if removed_ids:
        with open(trimal_file, "w") as handle:
            SeqIO.write(kept_records, handle, "fasta")

    return removed_ids, True


def process_sequence(sequence_id, hit_name):
    if sequence_id not in sequences:
        return f"Sequence {sequence_id} not found in FASTA. Skipping."

    # Accept hit names with or without the ".aligned" suffix.
    if hit_name.endswith(".aligned.fas"):
        alignment_name = hit_name
    elif hit_name.endswith(".aligned"):
        alignment_name = f"{hit_name}.fas"
    else:
        alignment_name = f"{hit_name}.aligned.fas"
    alignment_file = os.path.join(align_path, alignment_name)
    if not os.path.exists(alignment_file):
        raise FileNotFoundError(
            f"Alignment template {alignment_name} missing for {sequence_id}."
        )

    sequence_fasta = os.path.join(alignments_dir, f"{sequence_id}.fasta")
    with open(sequence_fasta, "w") as temp_fasta:
        SeqIO.write(sequences[sequence_id], temp_fasta, "fasta")

    run_file = os.path.join(alignments_dir, f"{sequence_id}.RUN")
    with open(run_file, "w") as run_fh:
        with open(alignment_file, "r") as align_fh:
            run_fh.write(align_fh.read())
        with open(sequence_fasta, "r") as temp_fasta:
            run_fh.write(temp_fasta.read())

    mafft_output = os.path.join(alignments_dir, f"{sequence_id}.new.aln")
    trimal_output = os.path.join(alignments_dir, f"{sequence_id}.trimal")

    mafft_cmd = [
        "mafft",
        "--auto",
        "--thread",
        str(mafft_threads),
        run_file,
    ]

    print(
        f"Running MAFFT for {sequence_id} with template {alignment_name}",
        flush=True,
    )
    with open(mafft_output, "w") as mafft_out:
        subprocess.run(mafft_cmd, stdout=mafft_out, check=True)

    trimal_cmd = [
        "trimal",
        "-in",
        mafft_output,
        "-out",
        trimal_output,
        "-gappyout",
    ]
    print(f"Running trimal for {sequence_id}", flush=True)
    subprocess.run(trimal_cmd, check=True)

    removed_ids, valid_alignment = clean_trimal_output(
        trimal_output, sequence_id, min_trimmed_residues
    )
    if removed_ids:
        print(
            f"Removed {len(removed_ids)} gap-only sequence(s) from {sequence_id}: "
            f"{', '.join(removed_ids)}",
            flush=True,
        )
    if not valid_alignment:
        return (
            f"Trimmed alignment for {sequence_id} skipped because it did not "
            f"contain at least {min_trimmed_residues} ungapped residue(s)."
        )

    return f"Alignment processing for {sequence_id} is complete."


def main():
    if sequence_hits.empty:
        print("No HMMscan hits to align. Exiting.", flush=True)
        return

    parallelism = max(1, global_threads // mafft_threads)
    print(
        f"Running alignments with {parallelism} worker(s); "
        f"{mafft_threads} thread(s) per MAFFT job.",
        flush=True,
    )

    with ThreadPoolExecutor(max_workers=parallelism) as executor:
        futures = {}
        for _, row in sequence_hits.iterrows():
            sequence_id = row["Sequence ID"]
            hit_name = row["HMMscan_hit_name"]
            future = executor.submit(process_sequence, sequence_id, hit_name)
            futures[future] = sequence_id

        for future in as_completed(futures):
            sequence_id = futures[future]
            try:
                msg = future.result()
            except Exception as exc:  # pragma: no cover
                raise RuntimeError(
                    f"Alignment failed for {sequence_id}"
                ) from exc
            print(msg, flush=True)

    for pattern in ["*.RUN", "*.fasta", "*.new.aln"]:
        for file in Path(alignments_dir).glob(pattern):
            file.unlink()
            print(f"Deleted intermediate file: {file}", flush=True)
    print("Alignments done.", flush=True)


log_handle = open(log_file, "w")
sys.stdout = log_handle
sys.stderr = log_handle
os.dup2(log_handle.fileno(), 1)
os.dup2(log_handle.fileno(), 2)
try:
    main()
finally:
    log_handle.close()
