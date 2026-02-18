#!/usr/bin/env python

import pandas as pd
from Bio import SearchIO
from Bio import BiopythonDeprecationWarning
import warnings
import sys
import os
 
mito_db_output = snakemake.input[0]
subtractive_db_output = snakemake.input[1]
output_file = snakemake.output[0]
log_file = snakemake.log[0]

def parse_blast_results(blast_output, categories, prefix):
    results = SearchIO.parse(blast_output, "blast-tab")
    rows = []

    for record in results:
        query = record.id
        row = {"Query": query}

        for cat in categories:
            row[f"{prefix}_{cat}_name"] = "None"
            row[f"{prefix}_{cat}_evalue"] = None
            row[f"{prefix}_{cat}_bitscore"] = None

        for hsp in record.hsps:
            category = hsp.hit_id.split("_")[0].upper()
            if category in categories and row[f"{prefix}_{category}_name"] == "None":
                row[f"{prefix}_{category}_name"] = hsp.hit_id
                row[f"{prefix}_{category}_evalue"] = hsp.evalue
                row[f"{prefix}_{category}_bitscore"] = hsp.bitscore

        rows.append(row)

    df = pd.DataFrame(rows).set_index("Query")

    for cat in categories:
        for suffix in ["name", "evalue", "bitscore"]:
            col = f"{prefix}_{cat}_{suffix}"
            if col not in df.columns:
                df[col] = None

    df = df.where(pd.notnull(df), None)  # Replace nan with None
    return df


def compare_hits(mito_df, sub_df, categories):
    """
    Compares bitscores for each pair of model organisms between MITO and SUB.
    """
    query_list = list(set(mito_df.index).union(set(sub_df.index)))
    result_rows = []

    # Handle missing queries
    for missing_query in set(query_list) - set(mito_df.index):
        missing_row = pd.DataFrame([{col: None for col in mito_df.columns}], index=[missing_query])
        mito_df = pd.concat([mito_df, missing_row])

    for missing_query in set(query_list) - set(sub_df.index):
        missing_row = pd.DataFrame([{col: None for col in sub_df.columns}], index=[missing_query])
        sub_df = pd.concat([sub_df, missing_row])

    for query in query_list:
        try:
            score = 0
            row = {"Query": query, "Subtractive_DB_score": 0, "Subtractive_DB_best_hit": "None"}

            # Compare bitscores
            for cat in categories:
                mito_bitscore = mito_df.loc[query, f"MITO_{cat}_bitscore"] if f"MITO_{cat}_bitscore" in mito_df.columns else None
                sub_bitscore = sub_df.loc[query, f"SUB_{cat}_bitscore"] if f"SUB_{cat}_bitscore" in sub_df.columns else None

                # if pd.notna(mito_bitscore) and pd.notna(sub_bitscore):
                #     if mito_bitscore > sub_bitscore:
                #         score += 1
                # elif pd.notna(mito_bitscore):
                #     score += 1

                if pd.notna(mito_bitscore) and pd.notna(sub_bitscore):
                    if mito_bitscore >= sub_bitscore:
                        score += 1
                elif pd.notna(mito_bitscore):
                    score += 1
                elif pd.notna(sub_bitscore):
                    score += 0


            row["Subtractive_DB_score"] = score

            if query in mito_df.index:
                mito_bitscores = mito_df.loc[query].filter(like="_bitscore")
                if mito_bitscores.notna().any():  # Ensure there are valid scores
                    mito_best_bitscore = mito_bitscores.max()
                    bitscore_col = mito_bitscores.idxmax()
                    # Translate bitscore column to name column
                    name_col = bitscore_col.replace("_bitscore", "_name")
                    # print(f"DEBUG: bitscore_col = {bitscore_col} for query {query}")
                    # print(f"DEBUG: name_col = {name_col} for query {query}")
                    # print(f"DEBUG: Available name columns = {mito_df.loc[query].filter(like='_name').index.tolist()}")
                    if name_col in mito_df.loc[query].filter(like="_name").index:
                        mito_best_name = mito_df.loc[query, name_col]
                    else:
                        # print(f"DEBUG: Name column {name_col} not found for query {query}.")
                        mito_best_name = "None"
                else:
                    mito_best_bitscore = None
                    mito_best_name = "None"
            else:
                mito_best_bitscore = None
                mito_best_name = "None"

            if query in sub_df.index:
                sub_bitscores = sub_df.loc[query].filter(like="_bitscore")
                if sub_bitscores.notna().any():  # Ensure there are valid scores
                    sub_best_bitscore = sub_bitscores.max()
                    bitscore_col = sub_bitscores.idxmax()
                    # Translate bitscore column to name column
                    name_col = bitscore_col.replace("_bitscore", "_name")
                    # print(f"DEBUG: bitscore_col = {bitscore_col} for query {query} (SUB)")
                    # print(f"DEBUG: name_col = {name_col} for query {query} (SUB)")
                    # print(f"DEBUG: Available name columns = {sub_df.loc[query].filter(like='_name').index.tolist()}")
                    if name_col in sub_df.loc[query].filter(like="_name").index:
                        sub_best_name = sub_df.loc[query, name_col]
                    else:
                        # print(f"DEBUG: Name column {name_col} not found for query {query}.")
                        sub_best_name = "None"
                else:
                    sub_best_bitscore = None
                    sub_best_name = "None"
            else:
                sub_best_bitscore = None
                sub_best_name = "None"


            # # Debug final values
            # print(f"DEBUG: Final score for Query {query} = {score}")
            # print(f"DEBUG: Best hit for Query {query} = {row['Subtractive_DB_best_hit']}")

            # Final comparison for best hit
            if mito_best_bitscore is not None and (sub_best_bitscore is None or mito_best_bitscore >= sub_best_bitscore):
                # Select the best hit from MitoDB if its bitscore is higher or equal
                row["Subtractive_DB_best_hit"] = mito_best_name
                # print(f"DEBUG: Assigned mito_best_name = {mito_best_name} for query {query}")
            elif sub_best_bitscore is not None:
                # Otherwise, select the best hit from SubDB
                row["Subtractive_DB_best_hit"] = sub_best_name
                # print(f"DEBUG: Assigned sub_best_name = {sub_best_name} for query {query}")
            else:
                # If neither has a valid hit, assign "None"
                row["Subtractive_DB_best_hit"] = "None"
                # print(f"DEBUG: No best hit found for query {query}")

            result_rows.append(row)
        except Exception as e:
            print(f"Error processing Query {query}: {e}")
            raise

    result_df = pd.DataFrame(result_rows).set_index("Query")

    return result_df


def main():
    categories = ["YEAST", "TETRA", "ARABI", "HUMAN", "ANDAGO", "ACANTH"]

    # Parse the MitoDB and SubtractiveDB outputs
    print("Starting parse_Sub/MitoDB")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
        mito_df = parse_blast_results(mito_db_output, categories, "MITO")
        sub_df = parse_blast_results(subtractive_db_output, categories, "SUB")

    # # Debug: Print parsed DataFrames
    # print("DEBUG: MITO DataFrame Columns:", mito_df.columns)
    # print("DEBUG: MITO DataFrame Head:")
    # print(mito_df.head())
    # print("DEBUG: SUB DataFrame Columns:", sub_df.columns)
    # print("DEBUG: SUB DataFrame Head:")
    # print(sub_df.head())

    # Compare hits and calculate scores
    result_df = compare_hits(mito_df, sub_df, categories)

    # Merge all results
    final_df = pd.concat([mito_df, sub_df, result_df], axis=1)

    # Replace NaN values with "None"
    final_df = final_df.where(final_df.notna(), "None")

    # Save final DataFrame to CSV
    final_df.to_csv(output_file)
    print("Sub/MitoDB parsing done.")

# Redirect stdout and stderr to log file
with open(log_file, "w") as log_file:
    sys.stdout = log_file
    sys.stderr = log_file

    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        raise
