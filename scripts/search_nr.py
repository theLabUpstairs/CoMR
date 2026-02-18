#!/usr/bin/env python

import pandas as pd
import os
import sys
import subprocess
import tempfile
from tqdm import tqdm  # Progress bar
import re

# Parameters from Snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]
taxonomy_output = snakemake.output[1]
nr_path = snakemake.params.nr_path
threads = snakemake.threads
block_size = snakemake.params.block_size
max_hits = snakemake.params.max_hits
sensitivity = snakemake.params.sensitivity
software_path = snakemake.params.software_path
taxonomy_enabled = snakemake.params.taxonomy_enabled
accession2taxid_path = snakemake.params.accession2taxid
nodes_path = snakemake.params.nodes
names_path = snakemake.params.names

def validate_file(file_path, description):
    """
    Validates that a file exists.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{description} not found: {file_path}")
    print(f"{description} found: {file_path}", flush=True)

def searchNR(input_file, nr_path, threads, output_file, software_path, max_hits, block_size, sensitivity, taxonomy_enabled):
    """
    Runs DIAMOND BLAST against the specified database. Includes taxonomy fields if enabled.
    """
    # Determine output format
    if taxonomy_enabled:
        fmt6_custom = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms stitle'
    else:
        fmt6_custom = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'

    nr_cmd = f"{software_path} blastp -p {threads} -d {nr_path} -q {input_file} " \
             f"-o {output_file} --outfmt {fmt6_custom} {sensitivity} -k {max_hits} --block-size {block_size}"

    print(f"Running DIAMOND command: {nr_cmd}", flush=True)
    result = os.system(nr_cmd)
    if result != 0:
        raise RuntimeError(f"DIAMOND command failed with exit code {result}")


def mapTaxonomy(output_file, accession2taxid, nodes, names, taxonomy_output):
    """
    Optimized version: Maps taxonomy to DIAMOND results, including kingdom-level information.
    """
    try:
        # Load DIAMOND results
        col_names = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle' 
        ]

        diamond_df = pd.read_csv(output_file, sep="\t", header=None, names=col_names)

        print(f"Loaded DIAMOND output with {len(diamond_df)} entries.", flush=True)
        print("Sample sseqid values:")
        print(diamond_df['sseqid'].head(10), flush=True)

        diamond_df['accession'] = diamond_df['sseqid'].astype(str).apply(
            lambda x: x.split('.')[0] if '.' in x else x
        )
        print("Extracted accessions:", diamond_df['accession'].head(10), flush=True)
        unique_accessions = set(diamond_df['accession'].dropna().unique())
        print(f"Unique accessions extracted: {len(unique_accessions)}", flush=True)

        # Exit early if no unique accessions
        if not unique_accessions:
            print("No unique accessions found. Exiting script.", flush=True)
            return

        # Create a temporary file for filtered accession-to-taxid
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".tsv") as temp_file:
            output_temp_file = temp_file.name
        print(f"Using temporary file: {output_temp_file}", flush=True)

        # Define data types for accession and taxid
        dtype_spec = {'accession': str, 'taxid': str}

        # Process chunks and write filtered results to a temporary file
        header_written = False
        chunk_count = 0

        tqdm_disable = not sys.stdout.isatty()
        for chunk in tqdm(pd.read_csv(
            accession2taxid,
            sep="\t",
            header=None,
            names=['accession', 'accession_version', 'taxid', 'gi'],
            usecols=['accession', 'taxid'],
            dtype=dtype_spec,
            chunksize=1000000  # Process 1,000,000 rows per chunk
        ), desc="Processing accession2taxid chunks", disable=tqdm_disable):
            chunk_count += 1
            chunk_filtered = chunk.query('accession in @unique_accessions')
            
            if not chunk_filtered.empty:
                chunk_filtered.to_csv(output_temp_file, sep="\t", mode='a', header=not header_written, index=False)
                header_written = True
                print(f"Filtered {len(chunk_filtered)} rows in chunk {chunk_count}.", flush=True)

        # Reload the final filtered file for further processing
        accession_df = pd.read_csv(output_temp_file, sep="\t", dtype=dtype_spec)
        print(f"Loaded filtered accession-to-taxid file with {len(accession_df)} entries.", flush=True)

        # Merge DIAMOND output with accession-to-taxid
        merged_df = diamond_df.merge(accession_df, on='accession', how='left')

        # Load nodes.dmp for taxonomic hierarchy
        nodes_df = pd.read_csv(nodes, sep="|", header=None, usecols=[0, 1, 2], 
                               names=['taxid', 'parent_taxid', 'rank'], engine="python")
        nodes_df['taxid'] = nodes_df['taxid'].astype(str).str.strip()
        nodes_df['parent_taxid'] = nodes_df['parent_taxid'].astype(str).str.strip()
        nodes_df['rank'] = nodes_df['rank'].astype(str).str.strip()
        print(f"Loaded nodes.dmp with {len(nodes_df)} entries.", flush=True)

        # Build taxonomic tree dictionaries
        tax_tree = nodes_df.set_index('taxid')['parent_taxid'].to_dict()
        rank_dict = nodes_df.set_index('taxid')['rank'].to_dict()

        # Precompute all paths to superkingdoms
        def precompute_kingdoms(tax_tree, rank_dict):
            """
            Precomputes the kingdom for each taxid in the taxonomic tree.
            """
            kingdom_map = {}

            def trace_to_kingdom(taxid):
                current_taxid = taxid
                visited = set()
                while current_taxid in tax_tree:
                    if current_taxid in visited:  # Break cycles (shouldn't happen in a valid tree)
                        return "Unknown"
                    visited.add(current_taxid)
                    rank = rank_dict.get(current_taxid, "no_rank")
                    if rank == 'superkingdom':  # Kingdom-level taxonomy
                        return current_taxid
                    current_taxid = tax_tree[current_taxid]
                return "Unknown"

            # Process all taxids in the tree
            for taxid in tax_tree:
                kingdom_map[taxid] = trace_to_kingdom(taxid)

            return kingdom_map


        # Build the precomputed kingdom map
        print("Precomputing kingdom mappings...", flush=True)
        kingdom_map = precompute_kingdoms(tax_tree, rank_dict)
        print(f"Precomputed kingdoms for {len(kingdom_map)} taxids.", flush=True)

        # Map taxids in merged_df to their kingdoms
        merged_df['kingdom_taxid'] = merged_df['taxid'].apply(lambda x: kingdom_map.get(str(x), "Unknown"))

        # Load names.dmp for taxonomic names
        names_df = pd.read_csv(
            names, sep="|", header=None, usecols=[0, 1, 3],
            names=['taxid', 'scientific_name', 'name_type'], engine="python"
        )

        # Clean up the columns
        names_df['taxid'] = names_df['taxid'].astype(str).str.strip()
        names_df['scientific_name'] = names_df['scientific_name'].astype(str).str.strip()
        names_df['name_type'] = names_df['name_type'].astype(str).str.strip()

        # Filter to keep only scientific names
        names_df = names_df[names_df['name_type'] == 'scientific name']

        # Build a dictionary for taxid to scientific_name mapping
        kingdom_names = names_df.set_index('taxid')['scientific_name'].to_dict()

        print(f"Loaded and filtered names.dmp with {len(names_df)} scientific names.", flush=True)


        # Replace kingdom_taxid with names
        merged_df['sskingdoms'] = merged_df['kingdom_taxid'].apply(lambda x: kingdom_names.get(x, "Unknown"))
        merged_df.drop(columns=['kingdom_taxid'], inplace=True)

        # Add the sscinames column
        merged_df['sscinames'] = merged_df['taxid'].apply(lambda x: kingdom_names.get(str(x), "Unknown"))

        # Add the staxids column explicitly
        merged_df['staxids'] = merged_df['taxid'].fillna("Unknown")

        print("Added kingdom-level information and scientific names to the merged data.", flush=True)

        # Fill missing values
        merged_df['staxids'] = merged_df['staxids'].fillna("Unknown")
        merged_df['sscinames'] = merged_df['sscinames'].fillna("Unknown")
        merged_df['sskingdoms'] = merged_df['sskingdoms'].fillna("Unknown")

        # Reorder the DataFrame to match the desired column order
        col_names = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
            'staxids', 'sscinames', 'sskingdoms', 'stitle'
        ]
        merged_df = merged_df[col_names]

        # Save the DataFrame to CSV
        merged_df.to_csv(taxonomy_output, sep="\t", index=False, header=False)
        print(f"Annotated results with kingdom saved to {taxonomy_output}.", flush=True)


        # Cleanup temporary file
        os.remove(output_temp_file)
        print(f"Temporary file {output_temp_file} removed.", flush=True)

    except Exception as e:
        print(f"Error in mapTaxonomy: {e}", flush=True)
        raise

# Redirect logging
with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = log_file

    os.dup2(log_file.fileno(), 1)  # Redirect stdout
    os.dup2(log_file.fileno(), 2)  # Redirect stderr

    try:
        # Validate input files
        validate_file(input_file, "Input FASTA file")
        validate_file(nr_path, "NR DIAMOND database")
        if not taxonomy_enabled:
            validate_file(accession2taxid_path, "Accession to TaxID mapping file")
            validate_file(nodes_path, "NCBI nodes.dmp file")
            validate_file(names_path, "NCBI names.dmp file")

        # Run DIAMOND
        print("Starting searchNR", flush=True)
        searchNR(input_file, nr_path, threads, output_file, software_path, max_hits, block_size, sensitivity, taxonomy_enabled)
        print("Search NR done.", flush=True)

        # Handle taxonomy mapping or placeholder output
        if not taxonomy_enabled:
            print("Starting taxonomy mapping", flush=True)
            mapTaxonomy(output_file, accession2taxid_path, nodes_path, names_path, taxonomy_output)
            print("Finished taxonomy mapping", flush=True)
        else:
            print("Skipping taxonomy mapping as taxonomy is included in the DIAMOND database.", flush=True)
            try:
                # Write placeholder taxonomy file
                with open(taxonomy_output, "w") as f:
                    f.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\t"
                            "send\tevalue\tbitscore\tstaxids\tsscinames\tsskingdoms\tstitle\n")
                    f.write("# Taxonomy included directly in the DIAMOND output.\n")
                print(f"Placeholder taxonomy file created: {taxonomy_output}", flush=True)
            except Exception as e:
                print(f"Error creating placeholder taxonomy file: {e}", flush=True)
                sys.exit(1)
        print("Taxonomy handling done.", flush=True)
    except Exception as e:
        print(f"Critical error: {e}", flush=True)
        sys.exit(1)
