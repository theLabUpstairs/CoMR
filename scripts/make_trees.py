import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple

# Snakemake inputs and outputs
alignments_dir = snakemake.input.alignments_dir  # Input directory with alignments
trees_dir = snakemake.output.trees_dir  # Output directory for trees
log_file = snakemake.log[0]  # Log file
global_threads = max(
    1, getattr(snakemake, "threads", snakemake.params.threads)
)
tree_threads = max(
    1, min(global_threads, snakemake.params.tree_threads)
)
min_non_gap_fraction = 0.15
min_informative_sites = 5

# Ensure the output directory exists
os.makedirs(trees_dir, exist_ok=True)


def assess_alignment(alignment_file: Path) -> Tuple[bool, str]:
    """
    Return whether an alignment has enough signal for tree inference.
    """
    sequences = []
    current = []
    with alignment_file.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                continue
            current.append(line)
    if current:
        sequences.append("".join(current))

    if not sequences:
        return False, "no sequences found"

    lengths = {len(seq) for seq in sequences}
    if len(lengths) != 1:
        return False, "inconsistent sequence lengths"
    aln_len = lengths.pop()
    if aln_len == 0:
        return False, "alignment length is zero"

    gap_chars = {"-", ".", "?"}
    total_chars = len(sequences) * aln_len
    non_gap_chars = 0
    informative_sites = 0

    for col in range(aln_len):
        col_chars = []
        for seq in sequences:
            c = seq[col]
            if c not in gap_chars:
                non_gap_chars += 1
                col_chars.append(c.upper())
        if len(col_chars) >= 2 and len(set(col_chars)) >= 2:
            informative_sites += 1

    if non_gap_chars == 0:
        return False, "all characters are gaps/missing"

    non_gap_fraction = non_gap_chars / total_chars
    if non_gap_fraction < min_non_gap_fraction:
        return False, (
            "too gap-heavy "
            f"(non-gap fraction {non_gap_fraction:.3f} < {min_non_gap_fraction:.3f})"
        )

    if informative_sites < min_informative_sites:
        return False, (
            "too few informative sites "
            f"({informative_sites} < {min_informative_sites})"
        )

    return True, "ok"


def build_tree(alignment_file: Path):
    alignment_path = str(alignment_file)
    output_prefix = os.path.join(trees_dir, alignment_file.stem)
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    is_usable, reason = assess_alignment(alignment_file)
    if not is_usable:
        return f"Skipped {alignment_file.name}: {reason}"

    iqtree_cmd = [
        "iqtree2",
        "-s",
        alignment_path,
        "-T",
        str(tree_threads),
        "-m",
        "LG+G",
        "-fast",
        "--prefix",
        output_prefix,
    ]

    print(
        f"Running IQ-TREE for {alignment_file.name} "
        f"using {tree_threads} thread(s)",
        flush=True,
    )
    try:
        subprocess.run(iqtree_cmd, check=True)
    except subprocess.CalledProcessError as exc:
        # Some alignments still fail in IQ-TREE despite pre-checks.
        # Keep processing the batch and report the failure as skipped.
        return (
            f"Skipped {alignment_file.name}: IQ-TREE failed "
            f"(exit code {exc.returncode})"
        )
    return f"IQ-TREE finished for {alignment_file.name}"


def main():
    alignment_files = sorted(Path(alignments_dir).glob("*.trimal"))
    if not alignment_files:
        print(f"No alignment files found in {alignments_dir}. Exiting.", flush=True)
        sys.exit(1)

    worker_count = max(1, global_threads // tree_threads)
    worker_count = min(worker_count, len(alignment_files))
    print(
        f"Building trees for {len(alignment_files)} alignments with "
        f"{worker_count} worker(s) x {tree_threads} thread(s).",
        flush=True,
    )

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(build_tree, alignment_file): alignment_file.name
            for alignment_file in alignment_files
        }
        for future in as_completed(futures):
            msg = future.result()
            print(msg, flush=True)

    print(f"Tree construction done. Results saved in {trees_dir}", flush=True)


log_handle = open(log_file, "w")
sys.stdout = log_handle
sys.stderr = log_handle
os.dup2(log_handle.fileno(), 1)
os.dup2(log_handle.fileno(), 2)
try:
    main()
finally:
    log_handle.close()
