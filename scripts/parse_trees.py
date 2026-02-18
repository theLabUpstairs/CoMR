#!/usr/bin/env python

import pandas as pd
import os
import sys
from ete3 import Tree

# Snakemake inputs and outputs
trees_dir = snakemake.input.trees_dir  # Input directory with tree files
output_csv = snakemake.output[0]  # Output CSV file
log_file = snakemake.log[0]  # Log file

# Redirect stdout and stderr to the log file
with open(log_file, "a") as log:
    sys.stdout = log
    sys.stderr = log
    os.dup2(log.fileno(), 1)  # Redirect system stdout
    os.dup2(log.fileno(), 2)  # Redirect system stderr

    def parse_trees(treefile):
        """Parse a single tree file and return the classification of the target sequence."""
        # Extract target sequence from the filename
        target_seq = os.path.basename(treefile).split(".treefile")[0]

        t = Tree(treefile)
        all_leaves = [node.name for node in t.traverse(strategy="postorder")]
        if target_seq not in all_leaves:
            return target_seq, "ERROR"

        t.set_outgroup(target_seq)

        for leaf in t:
            if leaf.name != target_seq:
                if leaf.name.split("x")[1].split("_")[0] == "MITO":
                    leaf.add_feature("typing", "1")
                elif leaf.name.split("x")[1].split("_")[0] == "OTHER":
                    leaf.add_feature("typing", "0")
            else:
                leaf.add_feature("typing", ["1", "0"])

        node_num = 1
        for node in t.traverse(strategy="postorder"):
            if not node.is_leaf() and "typing" not in node.features:
                state_list = [child.typing for child in node.get_children() if "typing" in child.features]
                if len(state_list) == 2:
                    if "0" in state_list and "1" in state_list:
                        node.add_feature("typing", ["1", "0"])
                    elif "0" in state_list:
                        node.add_feature("typing", "0")
                    else:
                        node.add_feature("typing", "1")
                elif len(state_list) == 3:
                    if state_list.count("1") >= 2:
                        node.add_feature("typing", "1")
                    elif state_list.count("0") >= 2:
                        node.add_feature("typing", "0")
                    else:
                        node.add_feature("typing", ["1", "0"])
                elif len(state_list) > 3:
                    if state_list.count("1") > state_list.count("0"):
                        node.add_feature("typing", "1")
                    elif state_list.count("0") > state_list.count("1"):
                        node.add_feature("typing", "0")
                    else:
                        node.add_feature("typing", ["1", "0"])

        for node in t.traverse(strategy="preorder"):
            if node.is_root():
                if node.typing == "1":
                    return target_seq, "MITO"
                elif node.typing == "0":
                    return target_seq, "OTHER"
                else:
                    return target_seq, "UNDEFINED"

    def main():
        # Aggregate results
        results = []

        # Process all *.treefile in the input directory
        print(f"Starting parse_trees")
        for treefile in os.listdir(trees_dir):
            if treefile.endswith(".treefile"):
                filepath = os.path.join(trees_dir, treefile)
                target_seq, result = parse_trees(filepath)
                if result != "ERROR":
                    results.append({"TreeFile": treefile, "Target": target_seq, "Tree_result": result})
                else:
                    print(f"Error: Target sequence ({target_seq}) not found in tree {treefile}.")

        # Write consolidated results to CSV
        pd.DataFrame(results).to_csv(output_csv, index=False)
        print(f"Parsing trees done.")

    if __name__ == "__main__":
        main()



# def parse_trees(treefile):

# ### MacClade algorithm ###
# # downpass - assigning ancestral states in a Fitch way
# # however, only the root as info from all of the tree
# # because we only care about the query sequence:
# # root the tree at the query sequence, then do a simple Fitch assignment
# # the root (query sequence) will have all info from the tree
# # after postorder traversal of the tree, see what state the query is given