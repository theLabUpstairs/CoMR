## Quick start Perun

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

   ```bash
   for Courtney: you will not be able to git clone from perun without setting up your github account there. Thus, an easiest way is to download the comr.zip directly from the github page, and upload it on perun. Once extracted, it will work the same.
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

   ```bash
   for Courtney: you can find the bundle on dardel: `/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/CoMR_DB_hmm/`
   ```

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


### Step 5 – Prepare the runtime config and input files

Copy the template and adjust any parameters you need:

```bash
cp config/config.yaml config/config_runtime.yaml
```
Note: the template is usable as it is, but you may want to adapt the `misc` options to adjust performance according to your dataset
(See recommendations in `README.md > ### Step 5 – Prepare the runtime config`)

The template leaves `fasta_files` empty on purpose; set it there if you want a
default sample or sample list, otherwise provide targets dynamically via
`--config fasta=sample or fasta=sample1,sample2`.

Keep the docker internal paths (`/mnt/databases/...` and `/mnt/software/...`) aligned
with the bind mounts shown below. Your actual local paths will be specified directly in the singularity command below.

Put your input fasta(s) in the `data/` folder (extension must be .fasta).


### Step 6 – Run CoMR on Perun

Perun uses Sun Grid Engine (SGE). Singularity is already available on the
cluster, so the workflow is to prepare a small submission script and hand it to
`qsub`. 

1. The snippet below assumes you stored the `.sif` image and
the CoMR folder under your project directory.

   ```bash
   for Courtney: you can find the image on dardel: `/cfs/klemming/projects/snic/tango2_lund_storage/nobackup/CoMR_images/`
   ```

2. Save the snippet below as `run_comr.sh`, edit the variables once, and
   submit the script via `qsub -q PerunMaster comr.sh`. 

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -o comr.$JOB_ID.log
#$ -M your.email@example.org
#$ -m be
#$ -pe threaded 48

. /etc/profile

# Declare local paths once so the bind list stays readable
# Replace "your_storage" placeholders and FASTA name with your actual values

COMR_ROOT=/path/to/your/CoMR
COMR_IMAGE=/path/to/comr_latest.sif
DB_DIR=/your/path/to/CoMR_DB_hmm
BLAST_DB=/your/path/to/blastdb
TAXONOMY=/your/path/to/taxonomy
TARGETP=/your/path/to/targetp-2.0
CORES=48
FASTA_BASENAME=fasta_basename   # matches a .fasta living under data/


cd "$COMR_ROOT"
mkdir -p "$COMR_ROOT/.inline_cache"

# for Courtney, you may have to adapt the first line according to perun
singularity exec --cleanenv \
  --bind "$DB_DIR:/mnt/databases:ro" \
  --bind "$BLAST_DB:/mnt/blastdb:ro" \
  --bind "$TARGETP:/mnt/software/targetp-2.0:ro" \
  --bind "$COMR_ROOT:/opt/CoMR" \
  --bind "$COMR_ROOT/.inline_cache:/opt/software/MitoFates/bin/modules/_Inline" \
  "$COMR_IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_BASENAME" enable_customdb=true

# update the placeholder variables with the directories you prepared above
```



