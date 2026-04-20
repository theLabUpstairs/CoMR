# CoMR — Comprehensive Mitochondrial Reconstructor

CoMR is a mitochondrial/mitochondrial-related organelles (MRO) proteome
prediction and reconstruction pipeline for both model organisms and eukaryotes
with atypical mitochondrial targeting signals. It combines multiple targeting
predictors, HMM searches, DIAMOND homology searches, and automated downstream
parsing to produce scored candidate lists.

## How to use this repository

| If you need… | See |
| --- | --- |
| Steps to **install** and **run** the workflow via prebuilt container images (Docker or Singularity), including databases/TargetP setup | Continue below |
| Understand how CoMR works and which steps the pipeline takes, and what outputs to expect | `README_CoMR.md`|
| Instructions to **build** the Docker/Singularity images yourself | `README_BUILD.md` |

---

## How to install and run CoMR

Before launching the containerized pipeline, you must fetch the licensed TargetP binary, 
download the CoMR databases, and run the helper script that retrieves third-party tools. 
The container bundles Snakemake and the runtime environment, 
but these external assets must live on the host and be bind-mounted at run time.

### Step 1 – Clone the repository

```bash
git clone https://github.com/theLabUpstairs/CoMR.git
cd CoMR
```

### Step 2 – Fetch MitoFates and Mitoprot II

```bash
bash scripts/fetch_third_party.sh
```

### Step 3 – Install TargetP (licensed)

Request a TargetP 2.0 license at [DTU Health Tech](https://services.healthtech.dtu.dk/services/TargetP-2.0/),
download the Linux archive, and extract it somewhere readable by your account and record the path (e.g. to `/your/path/to/targetp-2.0`)

```bash
tar -xzf targetp-2.0.Linux.tar.gz
chmod -R 755 targetp-2.0
```

The directory containing `bin/targetp` will later be mounted read-only inside
the container at `/mnt/software/targetp-2.0`.

### Step 4 – Download the CoMR databases

1. Download and extract (unzip) the CoMR database bundle (alignments, HMM profiles, MitoDB,
   SubtractedDB, UniProt) from the CoMR [Figshare](https://doi.org/10.17044/scilifelab.31361839)(e.g. to `/your/path/to/CoMR_DB_hmm`).

2. Make sure you have a DIAMOND-indexed NR database (`nr.dmnd`) and record the path (e.g. to `/your/path/to/blastdb`). 
   Follow the [DIAMOND documentation](https://github.com/bbuchfink/diamond/wiki) if you
   need to build one. Record whether taxonomy data were included. 
   Note: NR search can be disabled if no diamond indexed NR database is available.

3. If your DIAMOND-indexed NR database was built without taxonomy support, 
   make sure you have access [NCBI taxonomy dumps](https://ftp.ncbi.nih.gov/pub/taxonomy/) 
   (`prot.accession2taxid.gz`, `nodes.dmp`, `names.dmp`).
   CoMR needs them to annotate DIAMOND hits. Record the path (e.g. to `/your/path/to/taxonomy`)

You will bind-mount these directories into the container under `/mnt/databases`, `mnt/blastdb`
and `/mnt/taxonomy` in the docker/singularity command below.


### Step 5 – Prepare the runtime config

Copy the template:

```bash
cp config/config.yaml config/config_runtime.yaml
```
The template is usable as it is, but you may want to specify:

* Diamond options if needed. By default, `taxonomy_enabled` is set on True; if not, set it to False and mount the taxonomy DB inside the container at `/mnt/taxonomy/` in the docker/singularity command below.

* The `misc` to adjust for your system capacity if needed.

* `output_dir` if you want CoMR logs and declared outputs to be written outside
  the repository root.

#### Recommended `misc` values based on available CPUs/RAM

| Node profile | Total cores / RAM | threads | threads DIAMOND | DIAMOND block size | DIAMOND slots | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Workstation | 8 cores / 32 GB | 4-8 | 4-8 | 1 | 1 | Disable NR search |
| Mid HPC node | 32 cores / 128 GB | 24-32 | 8-16 | 1-2 | 1 | Leaves 8 cores free for MAFFT/DeepMito and budgets ~2-3 GB RAM per DIAMOND thread. |
| Large HPC node | ≥64 cores / ≥256 GB + fast scratch | 32-64 | 32-64 | 4+ | 2 | Only increase diamond slots if storage can handle two concurrent DIAMOND runs. Larger block sizes give DIAMOND more RAM to buffer queries. |

Rule of thumb: DIAMOND needs ~2-3 GB RAM per thread plus high I/O, so size
`threads_diamond` accordingly and still leave ≥4 cores for MAFFT/IQ-TREE so tree
building keeps pace once searches finish.

The template leaves `fasta_files` empty on purpose; set it there if you want a
default sample or sample list. You can also provide inputs dynamically at run
time with `--config fasta=sample.pep`, `--config fasta=sample1.pep,sample2.fasta`, or
`--config fasta=/path/to/sample.fasta`.

Keep the docker internal paths (`/mnt/databases/...` and `/mnt/software/...`) aligned
with the bind mounts shown below. Your actual local paths will be specified directly in the docker/singularity command below.

CoMR accepts FASTA inputs in two ways:

* Legacy shorthand: put your input file in `data/` and refer to it by filename,
  for example `--config fasta=proteins.fasta`
* Explicit paths: provide one or more FASTA paths directly, for example
  `--config fasta=/path/to/proteins.fasta` or set `fasta_files` in the runtime
  config to a list of paths. Example:

  ```yaml
  fasta_files:
    - /path/to/sample1.fasta
    - /path/to/sample2.fa
  ```

Accepted input extensions are `.fasta`, `.fa`, `.fas`, `.aa`, and `.pep`.
Internally, CoMR normalizes all accepted inputs to its own `.fasta`
intermediates before running TargetP and MitoFates.

Each input must resolve to a unique sample basename. For example,
`sample.fasta` and `sample.pep` would collide and should not be provided
together unless one is renamed.

CoMR also supports a configurable output root. Set `output_dir` in
`config_runtime.yaml` or pass it dynamically with
`--config output_dir=/path/to/comr_results`. This redirects declared workflow
outputs and logs, including `00_data_format_<FASTA>/`,
`01_analysis_original_<FASTA>/`, `02_analysis_parsed_<FASTA>/`,
`03_alignments_<FASTA>/`, `04_trees_<FASTA>/`, `05_CoMR_<FASTA>/`, and
`logs_<FASTA>/`, under the chosen directory. Internal runtime/cache files such
as `.snakemake/` and `.inline_cache/` remain under the CoMR installation
directory.


### Step 6 – Run CoMR from the container image

#### Docker

1. Pull the docker image:

```bash
docker pull ghcr.io/thelabupstairs/comr:latest
```

2. Declare your host paths once, then launch Snakemake:

```bash

# Declare local paths once 
# Replace "your_storage" placeholders, cores and FASTA name with your actual values

COMR_ROOT=/path/to/CoMR
DB_DIR=/path/to/CoMR_DB_hmm
NR_DMND=/path/to/blastdb/nr.dmnd
TAXONOMY=/path/to/taxonomy
TARGETP=/path/to/targetp-2.0
CORES=32
FASTA_INPUT=proteins.pep
OUTPUT_DIR=/path/to/comr_results
IMAGE=comr:latest

docker run --rm \
  --user "$(id -u)":"$(id -g)" \ # this allows you to run docker as a user
  -e HOME=/opt/CoMR \
  -v "$DB_DIR:/mnt/databases:ro" \  
  -v "$NR_DMND:/mnt/databases/nr.dmnd:ro" \
  -v "$TAXONOMY:/mnt/databases/taxonomy:ro" \
  -v "$TARGETP:/mnt/software/targetp-2.0:ro" \
  -v "$COMR_ROOT:/opt/CoMR" \
  -w /opt/CoMR \
  "$IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_INPUT" output_dir="$OUTPUT_DIR" #enable_customdb=true #enable_nr=false
```

- Use the command `id` to get additional group IDs needed for shared filesystems and
  pass them via repeated `--group-add <gid>`.
- When overriding multiple keys, include them all after a single `--config`
  flag (e.g. `--config fasta=your_sample.pep enable_nr=false` or
  `--config fasta=/path/to/your_sample.fas output_dir=/path/to/comr_results enable_nr=false`) because Snakemake only honors
  the last `--config` option.
- Append `enable_customdb=true` inside the `--config` block to activate the optional
  CustomDB search/parse steps.
- Append `enable_nr=false` to temporarily skip the NR DIAMOND search/parse stage if the database
  is not available.

#### Singularity / Apptainer

On HPC, it might be easier to rely on Singularity/Apptainer containers.

1. Download `CoMR.sif` from the CoMR [Figshare](https://doi.org/10.17044/scilifelab.31361839)

   ```bash
   curl -o CoMR.sif https://doi.org/10.17044/scilifelab.31361839/CoMR.sif
   ```

2. Execute Snakemake:

```bash

# Declare local paths once
# Replace "your_storage" placeholders, cores and FASTA name with your actual values

COMR_ROOT=/path/to/CoMR
COMR_IMAGE=/path/to/CoMR.sif
DB_DIR=/path/to/CoMR_DB_hmm
NR_DMND=/path/to/blastdb/nr.dmnd
TAXONOMY=/path/to/taxonomy
TARGETP=/path/to/targetp-2.0
CORES=32
FASTA_INPUT=proteins.pep
OUTPUT_DIR=/path/to/comr_results

cd "$COMR_ROOT"
mkdir -p "$COMR_ROOT/.inline_cache"

srun singularity exec \
  --bind "$DB_DIR:/mnt/databases:ro" \
  --bind "$NR_DMND:/mnt/blastdb:ro" \
  --bind "TAXONOMY=/mnt/databases/taxonomy:ro" \
  --bind "$TARGETP:/mnt/software/targetp-2.0:ro" \
  --bind "$COMR_ROOT:/opt/CoMR" \
  --bind "$COMR_ROOT/.inline_cache:/opt/software/MitoFates/bin/modules/_Inline" \
  "$COMR_IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_INPUT" output_dir="$OUTPUT_DIR" #enable_customdb=true #enable_nr=false
```

Include `enable_customdb=true` in the `--config` arguments if you need the
CustomDB workflow stages, or `enable_nr=false` if the NR database is not available.
As with Docker, pass multiple overrides in a single `--config` flag so that every
key takes effect.
