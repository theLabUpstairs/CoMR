import os
from pathlib import Path

import pandas as pd

output_file = snakemake.output.output_file
log_file = snakemake.log[0]
input_files = [Path(str(path)) for path in snakemake.input]

common_column_names = [
    "Sequence",
    "Sequence ID",
    "Query",
    "Index",
    "SequenceID",
    "Target",
]
standard_name = "SequenceID"
columns_to_keep = [
    "SequenceID",
    "RealSeqName",
    "tp_prediction",
    "tp_score",
    "tp_cleavage",
    "mf_score",
    "mf_seq",
    "mf_cleavage_pos",
    "mp_cleavagesite",
    "mp_MTS_sequence",
    "mp_score",
    "dm_location",
    "dm_score",
    "dm_GO_term",
    "Subtractive_DB_score",
    "Subtractive_DB_best_hit",
    "CustomDB_header",
    "CustomDB_hit_bitscore",
    "CustomDB_hit_evalue",
    "top_blast_hit_acc",
    "top_blast_hit_desc",
    "top_blast_hit_evalue",
    "top_blast_hit_bitscore",
    "top_blast_hit_pident",
    "top_blast_hit_taxid",
    "top_blast_hit_scientific_name",
    "top_blast_hit_kingdom",
    "top_nonhypo_hit_acc",
    "top_nonhypo_hit_desc",
    "top_nonhypo_hit_evalue",
    "top_nonhypo_hit_bitscore",
    "top_nonhypo_hit_pident",
    "top_nonhypo_hit_taxid",
    "top_nonhypo_hit_scientific_name",
    "top_nonhypo_hit_kingdom",
    "top_mito_hit_acc",
    "top_mito_hit_desc",
    "top_mito_hit_evalue",
    "top_mito_hit_bitscore",
    "top_mito_hit_pident",
    "top_mito_hit_taxid",
    "top_mito_hit_scientific_name",
    "top_mito_hit_kingdom",
    "Tree_result",
]


def log(message: str):
    with open(log_file, "a") as handle:
        handle.write(f"{message}\n")


def load_and_normalize(file_path: Path):
    log(f"Processing file: {file_path}")
    try:
        df = pd.read_csv(file_path)
    except Exception as exc:  # pragma: no cover
        log(f"Error reading {file_path}: {exc}")
        return None

    for col in common_column_names:
        if col in df.columns:
            df = df.rename(columns={col: standard_name})
            break
    else:
        log(f"Warning: {file_path} lacks a column that maps to {standard_name}. Skipping.")
        return None

    return df


def merge_tables(dataframes):
    if not dataframes:
        raise ValueError("No valid input files were parsed; cannot build table.")

    merged = dataframes[0]
    for df in dataframes[1:]:
        merged = pd.merge(merged, df, on=standard_name, how="outer")
    return merged


with open(log_file, "w") as log_handle:
    log_handle.write("Starting table generation\n")
    log_handle.write("Files provided by Snakemake:\n")
    for path in input_files:
        log_handle.write(f" - {path}\n")

dataframes = []
for csv_file in input_files:
    normalized = load_and_normalize(csv_file)
    if normalized is not None:
        dataframes.append(normalized)

result = merge_tables(dataframes)
result.fillna("No data", inplace=True)
result = result.reindex(columns=columns_to_keep, fill_value="No data")

try:
    result.to_csv(output_file, index=False)
    log(f"Table generation completed successfully. Output file: {output_file}")
except Exception as exc:
    log(f"Error writing output file: {exc}")
    raise
