## Quick start Dardel

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
   SubtractedDB, UniProt) from the [Figshare]() and
   extract it (e.g. to `/your/path/to/CoMR_DB_hmm`).

   For convenience, CoMR database bundle is available at `/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/CoMR_DB_hmm/`

2. The DIAMOND-indexed NR database (including taxonomy) (`nr.dmnd`) is available on dardel at: `/sw/data/diamond_databases/Blast/latest`


You will bind-mount these directories into the container under `/mnt/databases`, `mnt/blastdb` when running Snakemake.


### Step 5 – Prepare the runtime config and input files

Copy the template and adjust any parameters you need:

```bash
cp config/config_hpc.yaml config/config_runtime.yaml
```
Note: the template is usable as it is, but you may want to adapt the `misc` options to adjust performance according to your dataset.
(See recommendations in `README.md > ### Step 5 – Prepare the runtime config`)

The template leaves `fasta_files` empty on purpose; set it there if you want a
default sample or sample list, otherwise provide targets dynamically via
`--config fasta=sample or fasta=sample1,sample2`.

Keep the docker internal paths (`/mnt/databases/...` and `/mnt/software/...`) aligned
with the bind mounts shown below. Your actual local paths will be specified directly in the singularity command below.

Put your input fasta(s) in the `data/` folder (extension must be .fasta).


### Step 6 - Run CoMR

Most HPCs will have their own manual to use singularity containers.
For convenience, CoMR singularity image is available at `/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/CoMR_images/`


1. On dardel HPC, you need to convert a singularity image (.sif) to a sandbox with the build command:

```bash

ml PDC/24.11
ml singularity/4.2.0-cpeGNU-24.11

cd your_storage

singularity build --sandbox comr_latest /cfs/klemming/projects/snic/tango2_lund_storage/nobackup/CoMR_images/comr_latest.sif
```

2. Save the snippet below as `run_comr.sh`, edit the variables once, and
   launch it (`sbatch run_comr.sh`).

```bash
#!/bin/bash -l

#SBATCH -J comr
#SBATCH -t 24:00:00
#SBATCH -A naiss2025-5-253
#SBATCH -p main
#SBATCH -n 1
#SBATCH --mem=440G
#SBATCH --mail-user=your.email@biol.lu.se
#SBATCH --mail-type=END
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.out

ml PDC/24.11
ml singularity/4.2.0-cpeGNU-24.11

# Declare local paths once so the bind list stays readable
# Replace "your_storage" placeholders and FASTA name with your actual values
COMR_ROOT=/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/your_storage/CoMR
COMR_IMAGE=/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/your_storage/comr_latest
DB_DIR=/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/your_storage/CoMR_DB_hmm
BLAST_DB=/sw/data/diamond_databases/Blast/latest
TARGETP=/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/your_storage/targetp-2.0
CORES=128
FASTA_BASENAME=fasta_basename   # matches a .fasta under data/

cd "$COMR_ROOT"
mkdir -p "$COMR_ROOT/.inline_cache"

srun singularity exec -B /cfs/klemming \
  --bind "$DB_DIR:/mnt/databases:ro" \
  --bind "$BLAST_DB:/mnt/blastdb:ro" \
  --bind "$TARGETP:/mnt/software/targetp-2.0:ro" \
  --bind "$COMR_ROOT:/opt/CoMR" \
  --bind "$COMR_ROOT/.inline_cache:/opt/software/MitoFates/bin/modules/_Inline" \
  "$COMR_IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_BASENAME" enable_customdb=true


```