import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

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

# Ensure the output directory exists
os.makedirs(trees_dir, exist_ok=True)


def build_tree(alignment_file: Path):
    alignment_path = str(alignment_file)
    output_prefix = os.path.join(trees_dir, alignment_file.stem)
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

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
    subprocess.run(iqtree_cmd, check=True)
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
            try:
                msg = future.result()
            except Exception as exc:  # pragma: no cover
                raise RuntimeError(
                    f"IQ-TREE failed for {futures[future]}"
                ) from exc
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
