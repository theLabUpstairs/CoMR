# CoMR — Quick Setup and Run Guide

CoMR is a mitochondrial / mitochondrial-related organelle (MRO) proteome
prediction and reconstruction workflow for model organisms and eukaryotes with
atypical targeting signals. It combines targeting predictors, HMM searches,
DIAMOND homology searches, and downstream parsing to produce scored candidate
lists.

## Choose the right document

| If you want to... | Read |
| --- | --- |
| Install and run the packaged workflow | `README.md` |
| Understand pipeline logic, outputs, and scoring | `README_CoMR.md` |
| Build the container images yourself | `README_BUILD.md` |

## Before you start

The CoMR container includes Snakemake and the runtime software environment, but
you still need to provide a few host-side assets:

| Required on host | Why |
| --- | --- |
| CoMR repository clone | Workflow code and config |
| TargetP 2.0 binary | Licensed external dependency |
| CoMR database bundle | HMMs, alignments, MitoDB, SubtractedDB, UniProt |
| DIAMOND-formatted NR database | Optional but recommended NR homology search |
| NCBI taxonomy files | Needed if taxonomy is not embedded in `nr.dmnd` |

## Quick start

### 1. Clone CoMR

```bash
git clone https://github.com/theLabUpstairs/CoMR.git
cd CoMR
```

### 2. Fetch bundled third-party tools

```bash
bash scripts/fetch_third_party.sh
```

This retrieves MitoFates and Mitoprot II.

### 3. Install TargetP 2.0

Request a license from [DTU Health Tech](https://services.healthtech.dtu.dk/services/TargetP-2.0/),
download the Linux archive, and unpack it somewhere readable on the host, for
example `/your/path/to/targetp-2.0`.

```bash
tar -xzf targetp-2.0.Linux.tar.gz
chmod -R 755 targetp-2.0
```

CoMR expects this directory to be mounted inside the container as
`/mnt/software/targetp-2.0`.

### 4. Prepare the databases

#### 4.1 CoMR bundle

Download and extract the CoMR database bundle `CoMR_DB_hmm` from
[Figshare](https://doi.org/10.17044/scilifelab.31361839), including:

- `Alignments/`
- `Hmm_profile/`
- `SMD_MitoDB.fasta`
- `SMD_SubtractedDB.fasta`
- `uniprot_sprot.fasta`

Example host location: `/your/path/to/CoMR_DB_hmm`. If you wish to use an optional Custom Database, place it at the same location.

#### 4.2 NR database for DIAMOND

If you want the NR search stage, prepare `nr.dmnd` and record its path.
Example host location: `/your/path/to/blastdb/nr.dmnd`

Build NR with embedded taxonomy:

```bash
mkdir -p /your/path/to/blastdb
cd /your/path/to/blastdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -xzf new_taxdump.tar.gz names.dmp nodes.dmp

diamond makedb \
  --in nr.gz \
  --db nr \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnodes nodes.dmp \
  --taxonnames names.dmp
```

If you already have a DIAMOND NR database without taxonomy support, also keep
these NCBI taxonomy files available:

- `prot.accession2taxid.gz`
- `nodes.dmp`
- `names.dmp`

Example host location: `/your/path/to/taxonomy`

If you do not have an NR database, CoMR can still run with NR searches disabled using
`enable_nr=false`.

### 5. Create the runtime config

Copy the template:

```bash
cp config/config.yaml config/config_runtime.yaml
```

The template is intended to run with minimal edits. In practice, these are the
settings users most often change:

#### 5.1 Optional search settings

In case your Diamond-indexed NR database was built without taxonomy:

```yaml
diamond_search:
  taxonomy_enabled: False
```

Enable an optional CustomDB FASTA:

```yaml
database:
  customdb: "/mnt/databases/your_custom_db.fasta"
```

If `customdb` is set, also run CoMR with `enable_customdb=true` and be sure your Custom DB exists at `/your/path/to/CoMR_DB_hmm`.

Exclude specific taxa from NR hits by taxid:

```yaml
diamond_search:
  excluded_taxids:
    - 9606
    - 10090
```

Or exclude taxa from a one-taxid-per-line file:

```yaml
diamond_search:
  excluded_taxids_file: "config/exclusions/listoftaxids.txt"
```

Use taxon exclusion only when `taxonomy_enabled: True`.

#### 5.2 Resource settings

Adjust the `misc` section to match your hardware:

| Node profile | Total cores / RAM | `threads` | `threads_diamond` | `block_size` | `diamond_slots` | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Workstation | 8 cores / 32 GB | 4-8 | 4-8 | 1 | 1 | Often best to disable NR |
| Mid-size HPC | 32 cores / 128 GB | 24-32 | 8-16 | 1-2 | 1 | Leave CPU headroom for MAFFT and DeepMito |
| Large HPC | 64+ cores / 256+ GB | 32-64 | 32-64 | 4+ | 2 | Only raise `diamond_slots` if storage/RAM can keep up |

Rule of thumb: DIAMOND typically needs about 2-3 GB RAM per thread plus strong
I/O.

#### 5.3 Input and output settings

You can define input FASTA files in the config:

```yaml
fasta_files:
  - /path/to/sample1.fasta
  - /path/to/sample2.fa
```

Or pass them dynamically at runtime:

```bash
--config fasta=/path/to/sample.fasta
--config fasta=/path/to/sample1.fasta,/path/to/sample2.fasta
```

Accepted input extensions:

- `.fasta`
- `.fa`
- `.fas`
- `.aa`
- `.pep`

Each input must resolve to a unique sample basename. For example,
`sample.fasta` and `sample.pep` will collide.

To send declared outputs outside the repository, set:

```yaml
output_dir: "/path/to/comr_results"
```

This redirects workflow outputs and logs such as:

- `00_data_format_<FASTA>/`
- `01_analysis_original_<FASTA>/`
- `02_analysis_parsed_<FASTA>/`
- `03_alignments_<FASTA>/`
- `04_trees_<FASTA>/`
- `05_CoMR_<FASTA>/`
- `logs_<FASTA>/`

Internal runtime/cache directories such as `.snakemake/` and `.inline_cache/`
remain under the CoMR installation directory.

## Run CoMR

Set your local paths once:

```bash
COMR_ROOT=/path/to/CoMR
DB_DIR=/path/to/CoMR_DB_hmm
NR_DMND=/path/to/blastdb/nr.dmnd
TAXONOMY=/path/to/taxonomy
TARGETP=/path/to/targetp-2.0
CORES=32
FASTA_INPUT=/path/to/proteins.pep
OUTPUT_DIR=/path/to/comr_results
```

Keep the container-internal paths consistent with the config:

- databases under `/mnt/databases`
- NR under `/mnt/blastdb`
- taxonomy under `/mnt/taxonomy`
- TargetP under `/mnt/software/targetp-2.0`

### Docker

Pull the image:

```bash
docker pull ghcr.io/thelabupstairs/comr:latest
```

Run:

```bash
IMAGE=comr:latest

docker run --rm \
  --user "$(id -u)":"$(id -g)" \
  -e HOME=/opt/CoMR \
  -v "$DB_DIR:/mnt/databases:ro" \
  -v "$NR_DMND:/mnt/blastdb/nr.dmnd:ro" \
  -v "$TAXONOMY:/mnt/taxonomy:ro" \
  -v "$TARGETP:/mnt/software/targetp-2.0:ro" \
  -v "$COMR_ROOT:/opt/CoMR" \
  -w /opt/CoMR \
  "$IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_INPUT" output_dir="$OUTPUT_DIR"
```

Notes:

- Add `enable_customdb=true` inside the same `--config` block to enable CustomDB
- Add `enable_nr=false` inside the same `--config` block to skip NR
- If you override multiple keys, keep them after a single `--config`
- Use repeated `--group-add <gid>` if your filesystem permissions require
  additional groups

Example with multiple overrides:

```bash
snakemake --cores "$CORES" \
  --configfile config/config_runtime.yaml \
  --config fasta="$FASTA_INPUT" output_dir="$OUTPUT_DIR" enable_customdb=true
```

### Singularity / Apptainer (recommended on HPCs)

Download the image `CoMR.sif` from [Figshare](https://doi.org/10.17044/scilifelab.31361839)

Set:

```bash
COMR_IMAGE=/path/to/CoMR.sif
cd "$COMR_ROOT"
mkdir -p "$COMR_ROOT/.inline_cache"
```

On Slurm systems:

```bash
srun singularity exec \
  --bind "$DB_DIR:/mnt/databases:ro" \
  --bind "$NR_DMND:/mnt/blastdb:ro" \
  --bind "$TAXONOMY:/mnt/taxonomy:ro" \
  --bind "$TARGETP:/mnt/software/targetp-2.0:ro" \
  --bind "$COMR_ROOT:/opt/CoMR" \
  --bind "$COMR_ROOT/.inline_cache:/opt/software/MitoFates/bin/modules/_Inline" \
  "$COMR_IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_INPUT" output_dir="$OUTPUT_DIR"
```

On systems without Slurm:

```bash
singularity exec \
  --bind "$DB_DIR:/mnt/databases:ro" \
  --bind "$NR_DMND:/mnt/blastdb:ro" \
  --bind "$TAXONOMY:/mnt/taxonomy:ro" \
  --bind "$TARGETP:/mnt/software/targetp-2.0:ro" \
  --bind "$COMR_ROOT:/opt/CoMR" \
  --bind "$COMR_ROOT/.inline_cache:/opt/software/MitoFates/bin/modules/_Inline" \
  "$COMR_IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_INPUT" output_dir="$OUTPUT_DIR"
```

If your system uses Apptainer, replace `singularity exec` with `apptainer exec`.
If you are on an HPC system, check your local site documentation or ask your
system administrators for the correct module setup, scheduler integration, and
container invocation pattern.

## Practical reminders

- Record which version and location of every external asset you use
- Keep host paths and container paths consistent
- If NR taxonomy is not embedded in `nr.dmnd`, make sure taxonomy files are
  mounted separately
- If you disable NR, CoMR can still run, but NR-based search/parse stages will
  be skipped
- If you enable CustomDB, make sure the FASTA exists inside the mounted
  database directory and pass `enable_customdb=true`
