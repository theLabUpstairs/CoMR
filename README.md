# CoMR — Comprehensive Mitochondrial Reconstructor

CoMR is a mitochondrial/mitochondrial-related organelles (MRO) proteome
prediction and reconstruction pipeline for both model organisms and eukaryotes
with atypical mitochondrial targeting signals. It combines multiple targeting
predictors, HMM searches, DIAMOND homology searches, and automated downstream
parsing to produce scored candidate lists.

## How to use this repository

| If you need… | See |
| --- | --- |
| Steps to **run** the workflow via prebuilt container images (Docker or Singularity), including database/TargetP setup | Continue below |
| Understand how CoMR works and which steps the pipeline takes, and what outputs to expect | `README_CoMR.md`|
| Instructions to **build** the Docker/Singularity images yourself | `README_BUILD.md` |

---

## Quick start (prebuilt containers, recommended)

Even when you use a prebuilt container image, you must still fetch the licensed
TargetP binary, download the CoMR databases, and run the helper script that
retrieves third-party tools. The container bundles Snakemake and the runtime
environment, but these external assets must live on the host and be bind-mounted
at run time.

### Step 1 – Clone the repository

```bash
git clone https://github.com/theLabUpstairs/CoMR.git
cd CoMR
```

### Step 2 – Fetch MitoFates and Mitoprot II

Populate `third_party/`:

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

1. Download the CoMR database bundle (alignments, HMM profiles, MitoDB,
   SubtractedDB, UniProt) from the [Figshare](10.17044/scilifelab.31361839) and
   extract it (e.g. to `/your/path/to/CoMR_DB_hmm`).

2. Make sure you have a DIAMOND-indexed NR database (`nr.dmnd`) and record the path (e.g. to `/your/path/to/blastdb`). 
   Follow the [DIAMOND documentation](https://github.com/bbuchfink/diamond/wiki) if you
   need to build one. Record whether taxonomy data were included. 
   Note: NR search can be disabled if no diamond indexed NR database is available.

3. If your DIAMOND-indexed NR database was built without taxonomy support, 
   make sure you have access [NCBI taxonomy dumps](https://ftp.ncbi.nih.gov/pub/taxonomy/) 
   (`prot.accession2taxid.gz`, `nodes.dmp`, `names.dmp`).
   CoMR needs them to annotate DIAMOND hits. Record the path (e.g. to `/your/path/to/taxonomy`)

You will bind-mount these directories into the container under `/mnt/databases`, `mnt/blastdb`
and `/mnt/taxonomy` when running Snakemake.


### Step 5 – Prepare the runtime config

Copy the template and adjust any parameters you need:

```bash
cp config/config.yaml config/config_runtime.yaml
```
Note: the template is usable as it is, but you may want to specify:
- Diamond options if needed. By default, `taxonomy_enabled` is set on True; if not, set it to False and mount the taxonomy DB inside the container at `/mnt/taxonomy/` in the docker command below.
- The `misc` to adjust for your system capacity if needed.

Recommended `misc` values based on available CPUs/RAM (NR BLAST runs are HPC workloads; the “workstation” profile below is for debugging only, not full production runs):

| Node profile | Total cores / RAM | `misc. threads` | `misc. threads_diamond` | `diamond_search. block_size` | `misc. diamond_slots` | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Workstation (debug only) | 8 cores / 32 GB | 8 | 8 | 2 | 1 | Good for smoke tests, with NR search disabled.|
| Mid HPC node | 32 cores / 128 GB | 24/32 | 16 | 2-4 | 1 | Leaves 8 cores free for MAFFT/DeepMito and budgets ~2–3 GB RAM per DIAMOND thread. |
| Large HPC node | ≥64 cores / ≥256 GB + fast scratch | 32/48 | 32 | 6–8 | 2 | Only increase `diamond_slots` if storage can handle two concurrent DIAMOND runs. Larger block sizes give DIAMOND more RAM to buffer queries. |

Rule of thumb: DIAMOND needs ~2–3 GB RAM per thread plus high I/O, so size
`threads_diamond` accordingly and still leave ≥4 cores for MAFFT/IQ-TREE so tree
building keeps pace once searches finish.

The template leaves `fasta_files` empty on purpose; set it there if you want a
default sample or sample list, otherwise provide targets dynamically via
`--config fasta=sample or fasta=sample1,sample2`.

Keep the docker internal paths (`/mnt/databases/...` and `/mnt/software/...`) aligned
with the bind mounts shown below. Your actual local paths will be specified directly in the docker command below.

Put your input fasta(s) in the `data/` folder (extension must be .fasta).


### Step 6 – Run CoMR from the container image

#### Docker

1. Pull the docker image:

```bash
docker pull ghcr.io/thelabupstairs/comr:latest
```

2. Declare your host paths once, then launch Snakemake:

```bash

# Declare local paths once so the bind list stays readable
# Replace "your_storage" placeholders, cores and FASTA name with your actual values

COMR_ROOT=/path/to/CoMR
DB_DIR=/path/to/CoMR_DB_hmm
NR_DMND=/path/to/blastdb/nr.dmnd
TAXONOMY=/path/to/taxonomy
TARGETP=/path/to/targetp-2.0
CORES=16
FASTA_BASENAME=fasta
IMAGE=registry.example.org/comr:latest

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
    --config fasta="$FASTA_BASENAME" #enable_customdb=true #enable_nr=false
```

- Use the command `id` to get additional group IDs needed for shared filesystems and
  pass them via repeated `--group-add <gid>`.
- When overriding multiple keys, include them all after a single `--config`
  flag (e.g. `--config fasta=your.fasta enable_nr=false`) because Snakemake only honors
  the last `--config` option.
- Append `enable_customdb=true` inside the `--config` block to activate the optional
  CustomDB search/parse steps.
- Append `enable_nr=false` to temporarily skip the NR DIAMOND search/parse stage if the database
  is not available.

#### Singularity / Apptainer

1. Download `CoMR.sif` from the CoMR [Figshare](10.17044/scilifelab.31361839)

   ```bash
   curl -o CoMR.sif https://storage.example.org/comr/CoMR-1.0.sif
   ```

2. Execute Snakemake:

```bash

# Declare local paths once so the bind list stays readable
# Replace "your_storage" placeholders, cores and FASTA name with your actual values

COMR_ROOT=/path/to/CoMR
COMR_IMAGE=/path/to/CoMR.sif
DB_DIR=/path/to/CoMR_DB_hmm
NR_DMND=/path/to/blastdb/nr.dmnd
TAXONOMY=/path/to/taxonomy
TARGETP=/path/to/targetp-2.0
CORES=32
FASTA_BASENAME=fasta

cd "$COMR_ROOT"
mkdir -p "$COMR_ROOT/.inline_cache"

srun singularity exec \
  --bind "$DB_DIR:/mnt/databases:ro" \
  --bind "$BLAST_DB:/mnt/blastdb:ro" \
  --bind "$TARGETP:/mnt/software/targetp-2.0:ro" \
  --bind "$COMR_ROOT:/opt/CoMR" \
  --bind "$COMR_ROOT/.inline_cache:/opt/software/MitoFates/bin/modules/_Inline" \
  "$COMR_IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_BASENAME" #enable_customdb=true #enable_nr=false
```

Include `enable_customdb=true` in the `--config` arguments if you need the
CustomDB workflow stages, or `enable_nr=false` if the NR database is not available.
As with Docker, pass multiple overrides in a single `--config` flag so that every
key takes effect.
