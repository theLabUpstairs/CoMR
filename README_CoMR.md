# CoMR Workflow Guide

CoMR (Comprehensive Mitochondrial Reconstructor) is a Snakemake workflow that
aggregates evidence layers to identify and score mitochondrial or
mitochondrial-related organelle (MRO) proteins from predicted proteomes. 
Together, these steps let CoMR combine targeting predictions, curated
homology, large-scale BLAST-like evidence, and phylogenetic validation into a
single ranked list of mitochondrial/MRO candidates ready for expert review.


Instructions for installing dependencies or running the workflow live 
in the main README file.

Below is a description step-by-step of how CoMR proceeds.

<!-- ![CoMR pipeline](CoMR.png) -->

## Inputs, Outputs, and Directory Layout

- **Inputs**: One or more FASTA files `.fasta` placed in `data/`. 
- **Outputs**: For every FASTA basename the workflow produces a scored evidence
  table at `05_CoMR_<FASTA>/<FASTA>_CoMR_scored.csv`. The unscored aggregate
  (`05_CoMR_<FASTA>/<FASTA>_CoMR.csv`) plus intermediate logs/parses are kept
  alongside so you can trace how any candidate was supported.
- **Intermediate directories** (per FASTA):
  - `00_data_format_<FASTA>/`: Indexed FASTA plus ID-to-original name map and a
    version filtered to methionine-starting sequences.
  - `01_analysis_original_<FASTA>/`: Raw outputs from every external program
    (TargetP/Mitoprot/MitoFates/DeepMito summary tables, DIAMOND `.out` files, `hmmscan` tables).
  - `02_analysis_parsed_<FASTA>/`: Parsed CSVs that normalize each evidence
    source to a consistent schema expected by the merge/scoring rules.
  - `03_alignments_<FASTA>/` and `04_trees_<FASTA>/`: MAFFT/trimAl alignments
    and IQ-TREE phylogenies generated for every `hmmscan` hit family.
  - `logs_<FASTA>/`: Per-rule logs referenced by Snakemake's `log:` directive.

## Step 1 - FASTA Normalization

1. **`Fasta_format` (`scripts/data_index.py`)** reindexes every header with a
   random unique identifier and saves the mapping to
   `00_data_format_<FASTA>/<FASTA>_index_mapping.csv`. This isolates downstream
   tools from whitespace or special characters and allows reproducible joins.
2. **`Metstart` (`scripts/data_process.py`)** filters the indexed sequences to
   those that begin with methionine, which is required by the signal peptide
   predictors. The filtered FASTA is used whenever predictors need N-terminal
   context while the fully indexed FASTA stays available for searches.

## Step 2 - Targeting and Localization Predictors

For each FASTA the workflow runs four independent predictors and parses their
key scores:

- **`Targetp`**: predicts and score the presence of mitochondrial targeting
  peptides (mTP) and cleavage positions.
- **`Mitofates`**: reports probabilities that a presequence exists,
  predicted cleavage positions, and sequence motifs.
- **`Mitoprot`**: provides export probability, MTS sequence, and
  predicted cleavage sites.
- **`Deepmito`**: deep-learning localization predictions supported by
  the bundled UniProt FASTA.

Each script writes both the original output into `01_analysis_original_*`
and a simplified CSV into `02_analysis_parsed_*` so later steps can treat every
predictor uniformly.

## Step 3 - Homology and Profile Evidence

CoMR integrates multiple homology-based signals:

- **`HMMscan` + `HMM_parse`**: sequences are searched against the
  curated HMM library (`CoMR_mito.hmmdb`). Hits are parsed, linked back to the
  indexed IDs, and used later to decide which sequences enter the alignment/tree
  pipeline.
- **DIAMOND searches**:
  - `SubtractiveDB_blastp` runs a dual-subject search against both the curated
    MitoDB and SubtractedDB; `SubtractiveDB_parse` summarizes the best hits and
    scores into one CSV.
  - `CustomDB_blastp/parse` is optional and activated when
    `database.customdb` exists and `enable_customdb=true`. It lets you search a
    project-specific DIAMOND database (e.g., taxon-focused proteomes).
  - `NR_blastp/NR_parse` performs a large-scale search of the NCBI NR database
    with taxonomy reconstruction.

All DIAMOND rules honor the global `threads_diamond` limit and a shared
`diamond_slot` resource so you can cap how many heavy searches run at once.

## Step 4 - Profile-Guided Alignments and Phylogenies

- **`HMM_align`** takes the `hmmscan` hits, pulls the matching reference
  alignment template from `database.alignments`, and appends each query
  sequence. MAFFT (multi-threaded per job) builds the alignment, trimAl removes
  gappy regions. Temporary FASTA/MUSCLE files are cleaned automatically.
- **`HMM_trees`** runs IQ-TREE 2 on every trimmed alignment (`*.trimal`)
  produced above, using LG+G with the `-fast` heuristic. Parallelization is
  handled by dividing the global core budget among tree jobs.
- **`Trees_parse`** reads all `.treefile` results with ete3. By rooting each
  tree on the query sequence and propagating MITO/OTHER states via a Fitch-style
  traversal, it labels the placement as `MITO`, `OTHER`, or `UNDEFINED`. These
  classifications become another evidence layer during scoring.

## Step 5 - Evidence Integration and Scoring

1. **`Results` (`scripts/get_table.py`)** merges every parsed CSV into a single
   table. A sequence only appears once, and any source that lacked data is
   filled with "No data".
2. **`CoMR_scores` (`scripts/get_scores.py`)** applies the rubric defined in
   `CoMR_scorecards/*.score` (the default is `equal.score`). Each component can
   contribute a weighted value if its evidence exceeds the configured threshold:
   - SubtractiveDB score >= threshold (default: >=3)
   - TargetP predicts `mTP`
   - DIAMOND mito hit exists
   - MitoProt >=70% export probability 
   - MitoFates reports "Possessing_mitochondrial_presequence"
   - Tree parsing designates the sequence as `MITO`

   The script logs every per-sequence contribution and writes the final score in
   `CoMR_Score`. 

## Optional Workflows and Tuning Hooks

- **CustomDB**: Provide a `.dmnd` file under `database.customdb.path` in the
  runtime config and set `enable_customdb=true`. This branch is useful for
  lineage-specific proteomes not covered by the bundled DBs.
- **NR**: Enabled by default when `database.nr` points to a DIAMOND-indexed NR.
  Disable temporarily with `enable_nr=false` if you want a fast signal-peptide
  pass without whole-NR searches.
- **Threading**: `misc.threads`, `misc.threads_diamond`, `misc.align_threads_per_job`,
  and `misc.tree_threads_per_job` control the concurrency per stage. A global
  `diamond_slot` resource (default 1) protects disk-bound DIAMOND runs. 
  Adapt the config.yalm according to your system and needs.
- **Scorecards**: Drop additional `.score` files into `CoMR_scorecards/` to
  encode different weighting schemes (e.g., emphasize biochemical predictors or
  phylogenetic support). Point `scores.scorecard` in the config.yalm to the desired
  file before launching Snakemake.



