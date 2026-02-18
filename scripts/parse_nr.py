#!/usr/bin/env python

import pandas as pd
import sys
import re
import os

# File paths from Snakemake
diamond_output = snakemake.input[0]  # Input DIAMOND BLAST output
taxonomy_output = snakemake.input[1]  # Input taxonomy output (if needed)
output_file = snakemake.output[0]  # Output parsed CSV file
log_file = snakemake.log[0]  # Log file
taxonomy_enabled = snakemake.params.taxonomy_enabled
excluded_taxids = set()
excluded_taxids_file = snakemake.params.get("excluded_taxids_file", "")
if excluded_taxids_file and os.path.exists(excluded_taxids_file):
    with open(excluded_taxids_file) as f:
        excluded_taxids = set(line.strip() for line in f if line.strip())
else:
    excluded_taxids = set(str(tid) for tid in snakemake.params.get("excluded_taxids", []))



def parse_diamond(diamond_output, taxonomy_output, taxonomy_enabled, excluded_taxids):
    """
    Parses the DIAMOND BLAST output and extracts top hit information, including
    taxonomic information, top non-hypothetical hits, and mitochondrial hits.

    Parameters:
        diamond_output (str): Path to the DIAMOND BLAST output file.
        taxonomy_output (str): Path to the taxonomy output file.
        taxonomy_enabled (bool): Whether taxonomy is included in the DIAMOND output.

    Returns:
        pd.DataFrame: A DataFrame containing the parsed results with taxonomy.
    """
    # Define expected column names based on DIAMOND --outfmt
    col_names = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'staxids', 'sscinames', 'sskingdoms', 'stitle'
    ]

    # Load DIAMOND output or taxonomy mapping
    try:
        if not taxonomy_enabled:
            # Load the taxonomy mapping file
            diamond_df = pd.read_csv(taxonomy_output, sep="\t", header=None, names=col_names)
            print(f"Loaded taxonomy mapping with {len(diamond_df)} rows.", flush=True)
        else:
            # Load the DIAMOND BLAST output directly
            diamond_df = pd.read_csv(diamond_output, sep="\t", header=None, names=col_names)
            print(f"Loaded DIAMOND output with {len(diamond_df)} rows.", flush=True)
    except Exception as e:
        print(f"Error loading file: {e}", flush=True)
        sys.exit(1)
    print("Sample rows from loaded DIAMOND output:", flush=True)
    print(diamond_df.head(), flush=True)

    if excluded_taxids:
        def has_excluded_taxid(staxids):
            if pd.isna(staxids):
                return False
            for tid in str(staxids).split(";"):
                if tid in excluded_taxids:
                    return True
            return False

        before_count = len(diamond_df)
        diamond_df = diamond_df[~diamond_df['staxids'].apply(has_excluded_taxid)]
        print(
            f"Filtered out {before_count - len(diamond_df)} hits by excluded taxids.",
            flush=True,
        )

    # Ensure consistent ordering so the first hit per query is the best hit.
    diamond_df["evalue"] = pd.to_numeric(diamond_df["evalue"], errors="coerce")
    diamond_df["bitscore"] = pd.to_numeric(diamond_df["bitscore"], errors="coerce")
    diamond_df.sort_values(
        by=["qseqid", "bitscore", "evalue"],
        ascending=[True, False, True],
        kind="mergesort",
        inplace=True,
    )


    # Initialize empty rows for results
    rows = []

    # Group by query sequence (qseqid) to handle hits for each query
    for seq_id, group in diamond_df.groupby('qseqid'):
        top_hit = group.iloc[0]  # First row is the top hit
        top_hypo, mito_found = 0, 0

        # Initialize default values
        top_non_hypo_acc, top_non_hypo_desc, top_non_hypo_evalue, top_non_hypo_bitscore, top_non_hypo_pident = "None", "None", "None", "None", "None"
        top_non_hypo_taxid, top_non_hypo_scientific_name, top_non_hypo_kingdom = "None", "None", "None"
        top_mito_acc, top_mito_desc, top_mito_evalue, top_mito_bitscore, top_mito_pident = "None", "None", "None", "None", "None"
        top_mito_taxid, top_mito_scientific_name, top_mito_kingdom = "None", "None", "None"

        # Iterate through hits to find non-hypothetical and mitochondrial hits
        for _, hit in group.iterrows():
            # Safely extract hit description
            hit_desc = str(hit.get('stitle', "None")).lower() if pd.notna(hit.get('stitle')) else "none"

            # Debugging hit_desc
            # print(f"Debug: Processing hit description: {hit_desc}", flush=True)

            # Preprocess hit_desc to remove unwanted characters
            hit_desc_cleaned = re.sub(r"[^\w\s]", " ", hit_desc) # Remove punctuation for better matching
            hit_desc_cleaned = re.sub(r"\s+", " ", hit_desc_cleaned).strip()

            #print(f"Debug: Cleaned hit description: {hit_desc_cleaned}", flush=True)

            # Top non-hypothetical hit
            if "hypothetical" not in hit_desc_cleaned and top_hypo == 0:
                top_non_hypo_acc = hit['sseqid']
                top_non_hypo_desc = hit_desc
                top_non_hypo_evalue = hit['evalue']
                top_non_hypo_bitscore = hit['bitscore']
                top_non_hypo_pident = hit['pident']

                # if taxonomy_enabled:
                top_non_hypo_taxid = hit.get('staxids', "None")
                top_non_hypo_scientific_name = hit.get('sscinames', "None")
                top_non_hypo_kingdom = hit.get('sskingdoms', "None")
                # else:
                #     taxonomy_row = taxonomy_df[taxonomy_df["sseqid"] == hit['sseqid']]
                #     if not taxonomy_row.empty:
                #         top_non_hypo_taxid = taxonomy_row["staxids"].values[0]
                #         top_non_hypo_scientific_name = taxonomy_row["sscinames"].values[0]
                #         top_non_hypo_kingdom = taxonomy_row["sskingdoms"].values[0]
                top_hypo = 1

            # Top mitochondrial or hydrogenosomal hit
            if re.search(r"\b(mitochondrion|mitochondria|mitochondrial|hydrogenosome|hydrogenosomal|mitosome|mitosomal)\b", hit_desc_cleaned) and mito_found == 0:
                #print(f"Debug: Found mitochondrial hit: {hit_desc_cleaned}", flush=True)
                top_mito_acc = hit['sseqid']
                top_mito_desc = hit_desc
                top_mito_evalue = hit['evalue']
                top_mito_bitscore = hit['bitscore']
                top_mito_pident = hit['pident']

                # if taxonomy_enabled:
                top_mito_taxid = hit.get('staxids', "None")
                top_mito_scientific_name = hit.get('sscinames', "None")
                top_mito_kingdom = hit.get('sskingdoms', "None")
                # else:
                #     taxonomy_row = taxonomy_df[taxonomy_df["sseqid"] == hit['sseqid']]
                #     if not taxonomy_row.empty:
                #         top_mito_taxid = taxonomy_row["staxids"].values[0]
                #         top_mito_scientific_name = taxonomy_row["sscinames"].values[0]
                #         top_mito_kingdom = taxonomy_row["sskingdoms"].values[0]
                mito_found = 1


        # Append results for the current sequence
        rows.append({
            "Sequence ID": seq_id,
            "top_blast_hit_acc": top_hit['sseqid'],
            "top_blast_hit_desc": top_hit.get('stitle', "None"),
            "top_blast_hit_evalue": top_hit['evalue'],
            "top_blast_hit_bitscore": top_hit['bitscore'],
            "top_blast_hit_pident": top_hit['pident'],
            "top_blast_hit_taxid": top_hit.get('staxids', "None"),
            "top_blast_hit_scientific_name": top_hit.get('sscinames', "None"),
            "top_blast_hit_kingdom": top_hit.get('sskingdoms', "None"),
            "top_nonhypo_hit_acc": top_non_hypo_acc,
            "top_nonhypo_hit_desc": top_non_hypo_desc,
            "top_nonhypo_hit_evalue": top_non_hypo_evalue,
            "top_nonhypo_hit_bitscore": top_non_hypo_bitscore,
            "top_nonhypo_hit_pident": top_non_hypo_pident,
            "top_nonhypo_hit_taxid": top_non_hypo_taxid,
            "top_nonhypo_hit_scientific_name": top_non_hypo_scientific_name,
            "top_nonhypo_hit_kingdom": top_non_hypo_kingdom,
            "top_mito_hit_acc": top_mito_acc,
            "top_mito_hit_desc": top_mito_desc,
            "top_mito_hit_evalue": top_mito_evalue,
            "top_mito_hit_bitscore": top_mito_bitscore,
            "top_mito_hit_pident": top_mito_pident,
            "top_mito_hit_taxid": top_mito_taxid,
            "top_mito_hit_scientific_name": top_mito_scientific_name,
            "top_mito_hit_kingdom": top_mito_kingdom,

        })

    return pd.DataFrame(rows)

def main():
    # Parse the DIAMOND output with taxonomy
    df = parse_diamond(diamond_output, taxonomy_output, taxonomy_enabled, excluded_taxids)

    # Replace NaN or missing values with "None" for consistency
    df.fillna("None", inplace=True)

    # Save the parsed DataFrame to CSV
    df.to_csv(output_file, index=False)

# Redirect logging
with open(log_file, "w") as log_file:
    # Redirect Python-level stdout and stderr
    sys.stdout = log_file
    sys.stderr = log_file

    # Redirect system-level stdout and stderr
    os.dup2(log_file.fileno(), 1)  # Redirect fd=1 (stdout)
    os.dup2(log_file.fileno(), 2)  # Redirect fd=2 (stderr)

    print("Starting parse_nr", flush=True)
    main()
    print("Nr parsing done.", flush=True)
