#!/usr/bin/env python

import pandas as pd

# Snakemake parameters (defined at the beginning for clarity)
input_file = snakemake.input[0]  # Input CSV file
scorecard_file = snakemake.params.scorecard  # Path to the scorecard file
output_file = snakemake.output[0]  # Output CSV file
log_file = snakemake.log[0]  # Log file

# Define the columns to keep
columns_to_keep = [
    "SequenceID", "RealSeqName", "tp_prediction", "tp_score", "tp_cleavage", "mf_score", "mf_seq", 
    "mf_cleavage_pos", "mp_cleavagesite", "mp_MTS_sequence", "mp_score", "dm_location", "dm_score", 
    "dm_GO_term", "Subtractive_DB_score", "Subtractive_DB_best_hit", "CustomDB_header", 
    "CustomDB_hit_bitscore", "CustomDB_hit_evalue", "top_blast_hit_acc", "top_blast_hit_desc", 
    "top_blast_hit_evalue", "top_blast_hit_bitscore", "top_blast_hit_pident", "top_blast_hit_taxid", 
    "top_blast_hit_scientific_name", "top_blast_hit_kingdom", "top_nonhypo_hit_acc", "top_nonhypo_hit_desc", 
    "top_nonhypo_hit_evalue", "top_nonhypo_hit_bitscore", "top_nonhypo_hit_pident", 
    "top_nonhypo_hit_taxid", "top_nonhypo_hit_scientific_name", "top_nonhypo_hit_kingdom", "top_mito_hit_acc", 
    "top_mito_hit_desc", "top_mito_hit_evalue", "top_mito_hit_bitscore", "top_mito_hit_pident", 
    "top_mito_hit_taxid", "top_mito_hit_scientific_name", "top_mito_hit_kingdom", "Tree_result"
]


def scorecard(df, scorecard_file):
    # Read scorecard file
    scorecard = pd.read_csv(scorecard_file, sep="\t")
    scorecard.set_index("category", inplace=True)

    # Initialize scoring column
    df["CoMR_Score"] = "incomplete"

    component_counts = {
        "subtractiveDB": 0,
        "targetp": 0,
        "customDB": 0,
        "blast_mito": 0,
        "mitoprot": 0,
        "mitofates": 0,
        "tree": 0,
    }
    log_entries = []

    for index, row in df.iterrows():
        current_score = 0
        contributions = []

        # Subtractive DB score
        subtractive_score = safe_float(row["Subtractive_DB_score"])
        subDB_value = safe_float(scorecard.at["subtractiveDB", "value"])
        subDB_threshold = scorecard.at["subtractiveDB", "threshold"]

        if subDB_threshold == "DEFAULT":
            if subtractive_score >= 3:
                current_score += subDB_value
                contributions.append(f"subtractiveDB(+{subDB_value})")
                component_counts["subtractiveDB"] += 1
        else:
            if subtractive_score >= float(subDB_threshold):
                current_score += subDB_value
                contributions.append(f"subtractiveDB(+{subDB_value})")
                component_counts["subtractiveDB"] += 1

        # TargetP score
        targetp_value = safe_float(scorecard.at["targetp", "value"])
        targetp_score = row.get("tp_prediction", "")
        if targetp_score == "mTP":
            current_score += targetp_value
            contributions.append(f"targetp(+{targetp_value})")
            component_counts["targetp"] += 1

        # CustomDB score
        customDB_value = safe_float(scorecard.at["customDB", "value"])
        customDB_threshold = scorecard.at["customDB", "threshold"]

        if customDB_threshold == "DEFAULT":
            customDB_score = row.get("CustomDB_header", "No data")
            if customDB_score != "No data":
                current_score += customDB_value
                contributions.append(f"customDB(+{customDB_value})")
                component_counts["customDB"] += 1
        else:
            customDB_score = safe_float(row.get("CustomDB_hit_bitscore", 0)) 
            if customDB_score > float(customDB_threshold):
                current_score += customDB_value
                contributions.append(f"customDB(+{customDB_value})")
                component_counts["customDB"] += 1

        # Blast Mito score
        blastmito_value = safe_float(scorecard.at["blast_mito", "value"])
        blastmito_threshold = scorecard.at["blast_mito", "threshold"]

        if blastmito_threshold == "DEFAULT":
            blastmito_score = safe_float(row.get("top_mito_hit_bitscore", 0))
            if blastmito_score > 0:  # Only score if blastmito_score exists
                current_score += blastmito_value
                contributions.append(f"blast_mito(+{blastmito_value})")
                component_counts["blast_mito"] += 1

        # MitoProt score
        mitoprot_score = safe_float(row["mp_score"])
        mitoprot_value = safe_float(scorecard.at["mitoprot", "value"])
        mitoprot_threshold = scorecard.at["mitoprot", "threshold"]

        if mitoprot_threshold == "DEFAULT":
            if mitoprot_score >= 70:
                current_score += mitoprot_value
                contributions.append(f"mitoprot(+{mitoprot_value})")
                component_counts["mitoprot"] += 1
        else:
            if mitoprot_score >= float(mitoprot_threshold):
                current_score += mitoprot_value
                contributions.append(f"mitoprot(+{mitoprot_value})")
                component_counts["mitoprot"] += 1

        # MitoFates score
        mitofates_value = safe_float(scorecard.at["mitofates", "value"])
        mitofates_threshold = scorecard.at["mitofates", "threshold"]

        if mitofates_threshold == "DEFAULT":
            mitofates_score = row.get("mf_seq", "")
            if mitofates_score == "Possessing_mitochondrial_presequence":
                current_score += mitofates_value
                contributions.append(f"mitofates(+{mitofates_value})")
                component_counts["mitofates"] += 1
        else:
            mitofates_score = safe_float(row.get("mf_seq", 0))
            if mitofates_score >= float(mitofates_threshold):
                current_score += mitofates_value
                contributions.append(f"mitofates(+{mitofates_value})")
                component_counts["mitofates"] += 1

        # Tree search score
        tree_value = safe_float(scorecard.at["tree", "value"])
        tree_threshold = scorecard.at["tree", "threshold"]

        if tree_threshold == "DEFAULT":
            tree_score = row.get("Tree_result", "")
            if tree_score == "MITO":
                current_score += tree_value
                contributions.append(f"tree(+{tree_value})")
                component_counts["tree"] += 1

        # Assign final score to the DataFrame
        df.at[index, "CoMR_Score"] = current_score

        log_entries.append({
            "SequenceID": row.get("SequenceID", "NA"),
            "RealSeqName": row.get("RealSeqName", "NA"),
            "CoMR_Score": current_score,
            "Contributions": ", ".join(contributions) if contributions else "None",
        })

    # Reorder columns to place CoMR_Score after SequenceID and RealSeqName
    reordered_columns = (
        ["SequenceID", "RealSeqName", "CoMR_Score"] +
        [col for col in df.columns if col not in {"SequenceID", "RealSeqName", "CoMR_Score"}]
    )
    df = df[reordered_columns]
    return df, log_entries, component_counts


def safe_float(value, default=0.0):
    """
    Safely convert a value to float, returning a default value if conversion fails.
    """
    if value == "No data" or pd.isna(value):
        return default
    try:
        return float(value)
    except ValueError:
        return default

def main():
    # Read input data
    df = pd.read_csv(input_file)

    # Reindex DataFrame to ensure all columns in `columns_to_keep` are present
    df = df.reindex(columns=columns_to_keep, fill_value="No data")

    # Apply scoring
    scored_df, log_entries, component_counts = scorecard(df, scorecard_file)

    # Save the output
    scored_df.to_csv(output_file, index=False)

    # Write per-sequence contributions to the log without inflating file size
    with open(log_file, "w") as log:
        log.write("Scoring complete. Output written to: {}\n".format(output_file))
        log.write("Per-sequence summary (component contributions shown only when they added to the score):\n")
        for entry in log_entries:
            log.write(
                "{seq} ({name}): Score={score}; contributors={contributors}\n".format(
                    seq=entry["SequenceID"],
                    name=entry["RealSeqName"],
                    score=entry["CoMR_Score"],
                    contributors=entry["Contributions"],
                )
            )

        log.write("\nComponent usage summary:\n")
        for component, count in component_counts.items():
            log.write(f"{component}: {count} sequences\n")

if __name__ == "__main__":
    main()
