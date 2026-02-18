# Snakefile

from snakemake.exceptions import WorkflowError

configfile: "config/config.yaml"


def _coerce_list(value):
    """Normalize comma-delimited strings or iterables into a clean list."""
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    raise WorkflowError(f"Unsupported FASTA specification type: {type(value)}")


def _resolve_fasta_samples(cfg):
    """Return the FASTA sample names coming from CLI flags or runtime config."""
    cli_samples = _coerce_list(cfg.get("fasta"))
    if cli_samples:
        return cli_samples

    samples = _coerce_list(cfg.get("fasta_files"))
    if samples:
        return samples

    raise WorkflowError(
        "No FASTA input specified. Provide `--config fasta=sample1,sample2` "
        "or configure `fasta_files` in the runtime config."
    )


def _coerce_bool(value):
    """Interpret heterogeneous truthy/falsy config values as booleans."""
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "on"}:
            return True
        if lowered in {"0", "false", "no", "off"}:
            return False
    raise WorkflowError(f"Unable to interpret boolean value from {value!r}")


def _resolve_customdb_settings(cfg):
    """Extract CustomDB path/enabled state while honoring CLI overrides."""
    db_cfg = cfg.get("database", {})
    raw_custom = db_cfg.get("customdb")

    custom_path = None
    default_enabled = False

    if isinstance(raw_custom, dict):
        custom_path = raw_custom.get("path")
        default_enabled = raw_custom.get("enabled", False)
    else:
        custom_path = raw_custom
        default_enabled = False

    override = cfg.get("enable_customdb")
    if override is None:
        override = cfg.get("customdb_enabled")

    if override is not None:
        enabled = _coerce_bool(override)
    else:
        enabled = default_enabled

    if enabled and not custom_path:
        raise WorkflowError(
            "CustomDB is enabled but `database.customdb` is not configured with a path."
        )

    return custom_path, enabled


def _resolve_nr_settings(cfg):
    """Extract NR database path/enabled flag with validation."""
    db_cfg = cfg.get("database", {})
    nr_path = db_cfg.get("nr")

    override = cfg.get("enable_nr")
    if override is None:
        override = cfg.get("nr_enabled")

    if override is not None:
        enabled = _coerce_bool(override)
    else:
        enabled = True

    if enabled and not nr_path:
        raise WorkflowError(
            "NR workflow is enabled but `database.nr` is not configured. "
            "Set `enable_nr=false` or provide the NR DIAMOND path."
        )

    return nr_path, enabled


def _fasta_targets():
    """Expose FASTA targets for use inside expand/target rules."""
    return _resolve_fasta_samples(config)


def _customdb_path():
    """Convenience accessor so rules can grab the CustomDB path."""
    return _resolve_customdb_settings(config)[0]


def _customdb_enabled():
    """True when the optional CustomDB search should be triggered."""
    return _resolve_customdb_settings(config)[1]


def _nr_path():
    """Resolve the NR DIAMOND database path."""
    return _resolve_nr_settings(config)[0]


def _nr_enabled():
    """True when NR searches are requested."""
    return _resolve_nr_settings(config)[1]


def _diamond_threads():
    """Number of threads reserved for DIAMOND rules."""
    return config.get("misc", {}).get("threads_diamond", 1)


def _diamond_slot_limit():
    """Maximum number of parallel DIAMOND jobs allowed."""
    return config.get("misc", {}).get("diamond_slots", 1)


if "workflow" in globals():  # pragma: no cover - provided by Snakemake runtime
    workflow.global_resources = getattr(workflow, "global_resources", {})
    workflow.global_resources.setdefault("diamond_slot", _diamond_slot_limit())


def _results_inputs(fasta):
    """Build dynamic input list for Results rule with optional CustomDB output."""
    inputs = [
        f"00_data_format_{fasta}/{fasta}_index_mapping.csv",
        f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_targetp.csv",
        f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_mitofates.csv",
        f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_mitoprot.csv",
        f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_deepmito.csv",
    ]

    if _customdb_enabled():
        inputs.append(
            f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_CustomDB.csv"
        )

    if _nr_enabled():
        inputs.append(
            f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_nr.csv"
        )

    inputs.extend(
        [
            f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_DB.csv",
            f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_hmmscan.csv",
            f"02_analysis_parsed_{fasta}/{fasta}_CoMR_parse_trees.csv",
        ]
    )
    return inputs


# Final target rule: require a scored CoMR table for every FASTA requested.
rule all:
    input:
        lambda wildcards: expand(
            "05_CoMR_{FASTA}/{FASTA}_CoMR_scored.csv",
            FASTA=_fasta_targets(),
        )

# Normalize input FASTA files by indexing headers and building an ID map.
rule Fasta_format:
    input:
        "data/{FASTA}.fasta"
    output:
        output1="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta",
        output2="00_data_format_{FASTA}/{FASTA}_index_mapping.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_data_index.log"
    script:
        "scripts/data_index.py"

# Filter sequences to those that start with methionine so signal peptides are comparable.
rule Metstart:
    input:
        "00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta"
    output:
        "00_data_format_{FASTA}/{FASTA}_CoMR_metstart.fasta"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_data_process.log"
    script:
        "scripts/data_process.py"

# Predict targeting peptides with TargetP2 and summarize raw + parsed outputs.
rule Targetp:
    input:
        input1="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta",
        input2="00_data_format_{FASTA}/{FASTA}_CoMR_metstart.fasta"
    output:
        output1="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_targetp.csv",
        output2="01_analysis_original_{FASTA}/{FASTA}_CoMR_summary.targetp2"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_and_parse_targetp.log"
    params:
        software_path=config["software"]["targetp"],
        batch_path = config["misc"]["batch_targetp"]
    script:
        "scripts/search_and_parse_targetp.py"

# Score sequences with MitoFates and capture both parsed probabilities and the raw summary.
rule Mitofates:
    input:
        input1="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta",
        input2="00_data_format_{FASTA}/{FASTA}_CoMR_metstart.fasta"
    output:
        output1="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_mitofates.csv",
        output2="01_analysis_original_{FASTA}/{FASTA}_CoMR_summary.mitofates"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_and_parse_mitofates.log"
    params:
        software_path=config["software"]["mitofates"]
    script:
        "scripts/search_and_parse_mitofates.py"

# Evaluate mitochondrial localization using MitoProt and capture parsed metrics.
rule Mitoprot:
    input:
        input1="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta",
        input2="00_data_format_{FASTA}/{FASTA}_CoMR_metstart.fasta"
    output:
        output1="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_mitoprot.csv",
        output2="01_analysis_original_{FASTA}/{FASTA}_CoMR_summary.mitoprot"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_and_parse_mitoprot.log"
    params:
        software_path=config["software"]["mitoprot"]
    script:
        "scripts/search_and_parse_mitoprot.py"

# Run DeepMito using the supplied UniProt DB; record parsed predictions and raw output.
rule Deepmito:
    threads: config["misc"].get("deepmito_threads", config["misc"].get("threads", 1))
    input:
        input1="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta",
        input2="00_data_format_{FASTA}/{FASTA}_CoMR_metstart.fasta"
    output:
        output1="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_deepmito.csv",
        output2="01_analysis_original_{FASTA}/{FASTA}_CoMR_summary.deepmito"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_and_parse_deepmito.log"
    params:
        software_path=config["software"]["deepmito"],
        db_path=config["database"]["uniprot"],
    script:
        "scripts/search_and_parse_deepmito.py"

# Search proteins against the curated HMM library.
rule HMMscan:
    threads: config["misc"].get("threads", 1)
    input:
        "00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta"
    output:
        "01_analysis_original_{FASTA}/{FASTA}_CoMR_hmmscan.out"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_hmm.log"
    params:
        software_path=config["software"]["hmmscan"],
        db_path=config["database"]["hmmdb"],
        threads = config["misc"]["threads"],
    script:
        "scripts/search_hmm.py"


# Perform DIAMOND blastp against mito/subtractive databases in one pass.
rule SubtractiveDB_blastp:
    threads: _diamond_threads()
    resources:
        diamond_slot=1
    input:
        "00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta"
    output:
        mito_output="01_analysis_original_{FASTA}/{FASTA}_CoMR_MitoDB.out",
        sub_output="01_analysis_original_{FASTA}/{FASTA}_CoMR_SubstractiveDB.out"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_DB.log"
    params:
        software_path=config["software"]["diamond"],
        subdb_path=config["database"]["subdb"],
        mitodb_path=config["database"]["mitodb"],
    script:
        "scripts/search_DB.py"

# Optionally search a user-provided DIAMOND database.
rule CustomDB_blastp:
    threads: _diamond_threads()
    resources:
        diamond_slot=1
    input:
        "00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta"
    output:
        "01_analysis_original_{FASTA}/{FASTA}_CoMR_CustomDB.out"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_customDB.log"
    params:
        software_path=config["software"]["diamond"],
        customdb_path=lambda wildcards: _customdb_path(),
    script:
        "scripts/search_customDB.py"

# Convert hmmscan table to a CSV that links hits back to the indexed FASTA IDs.
rule HMM_parse:
    input:
        hmm="01_analysis_original_{FASTA}/{FASTA}_CoMR_hmmscan.out",
        fasta="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta"

    output:
        csv="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_hmmscan.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_parse_hmm.log"
    script:
        "scripts/parse_hmm.py"


# Summarize both mito and subtractive DIAMOND searches into one CSV.
rule SubtractiveDB_parse:
    input:
        "01_analysis_original_{FASTA}/{FASTA}_CoMR_MitoDB.out",
        "01_analysis_original_{FASTA}/{FASTA}_CoMR_SubstractiveDB.out"
    output:
        "02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_DB.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_parse_DB.log"
    script:
        "scripts/parse_DB.py"

# Translate optional CustomDB BLAST hits to the standard parsed CSV schema.
rule CustomDB_parse:
    input:
        "01_analysis_original_{FASTA}/{FASTA}_CoMR_CustomDB.out"
    output:
        "02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_CustomDB.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_parse_CustomDB.log"
    script:
        "scripts/parse_customDB.py"

# Launch DIAMOND against the NCBI NR database with optional taxonomy tracing.
rule NR_blastp:
    threads: _diamond_threads()
    resources:
        diamond_slot=1
    input:
        "00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta"
    output:
        output_file="01_analysis_original_{FASTA}/{FASTA}_CoMR_nr.out",
        taxonomy_output="01_analysis_original_{FASTA}/{FASTA}_CoMR_nr_taxonomy.tsv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_search_nr.log"
    params:
        software_path=config["software"]["diamond"],
        nr_path=lambda wildcards: _nr_path(),
        block_size=config['diamond_search']['block_size'],
        max_hits=config['diamond_search']['max_hits'],
        sensitivity=config["diamond_search"]["sensitivity"],
        taxonomy_enabled=config["diamond_search"]["taxonomy_enabled"],
        accession2taxid=config["taxonomy"]["accession2taxid"],
        nodes=config["taxonomy"]["nodes"],
        names=config["taxonomy"]["names"]
    script:
        "scripts/search_nr.py"

# Convert NR DIAMOND (and optional taxonomy) results into the parsed CSV.
rule NR_parse:
    input:
        diamond_output="01_analysis_original_{FASTA}/{FASTA}_CoMR_nr.out",
        taxonomy_output="01_analysis_original_{FASTA}/{FASTA}_CoMR_nr_taxonomy.tsv"
    output:
        output_file="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_nr.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_parse_nr.log"
    params:
        taxonomy_enabled=config["diamond_search"]["taxonomy_enabled"],
        excluded_taxids=config["diamond_search"].get("excluded_taxids", []),
        excluded_taxids_file=config["diamond_search"].get("excluded_taxids_file", "")
    script:
        "scripts/parse_nr.py"



# Align all hmmscan-supported families so downstream tree building can use trimmed MSAs.
rule HMM_align:
    threads: config["misc"].get("threads", 1)
    input:
        fasta="00_data_format_{FASTA}/{FASTA}_CoMR_indexed.fasta",
        csv="02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_hmmscan.csv"
    output:
        alignments_dir=directory("03_alignments_{FASTA}/")
    params:
        align_path=config["database"]["alignments"],
        threads=config["misc"].get("threads", 1),
        align_threads=config["misc"].get("align_threads_per_job", 1),
        min_trimmed_residues=config["misc"].get("min_trimmed_residues", 1)
    log:
        "logs_{FASTA}/{FASTA}_CoMR_align.log"
    script:
        "scripts/make_align.py"


# Build phylogenetic trees for each aligned family to validate orthology.
rule HMM_trees:
    threads: config["misc"].get("threads", 1)
    input:
        alignments_dir="03_alignments_{FASTA}/"
    output:
        trees_dir=directory("04_trees_{FASTA}/")
    params:
        threads=config["misc"].get("threads", 1),
        tree_threads=config["misc"].get("tree_threads_per_job", 1)
    log:
        "logs_{FASTA}/{FASTA}_CoMR_trees.log"
    script:
        "scripts/make_trees.py"


# Extract placement/monophyly summaries from every generated treefile.
rule Trees_parse:
    """
    Iterates through all *.treefile in 04_trees_{FASTA} and concatenates results into a single CSV.
    """
    input:
        trees_dir="04_trees_{FASTA}"
    output:
        "02_analysis_parsed_{FASTA}/{FASTA}_CoMR_parse_trees.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_parse_trees.log"
    script:
        "scripts/parse_trees.py"


# Merge every parsed evidence source into the master CoMR table per FASTA.
rule Results:
    input:
        lambda wildcards: _results_inputs(wildcards.FASTA)
    output:
        output_file = "05_CoMR_{FASTA}/{FASTA}_CoMR.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_table.log"
    script:
        "scripts/get_table.py"



# Apply the scoring rubric onto the merged table to rank mitochondrial candidates.
rule CoMR_scores:
    input:
        "05_CoMR_{FASTA}/{FASTA}_CoMR.csv"
    output:
        "05_CoMR_{FASTA}/{FASTA}_CoMR_scored.csv"
    log:
        "logs_{FASTA}/{FASTA}_CoMR_scored.log"
    params:
        scorecard = config["scores"]["scorecard"]
    script:
        "scripts/get_scores.py"
