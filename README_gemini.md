## Quick start Gemini

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

Targetp is already available on gemini at `/home/gemini/software/targetp-2.0`

### Step 4 – the CoMR databases

1. The CoMR databases are already available for all members of the lab group at this path: `/media/local_4/CoMR_DB_hmm`

2. The DIAMOND-indexed NR database (including taxonomy) (`nr.dmnd`) is available for all members of the lab group at this path: `/media/local_4/blastdb`

You will bind-mount these directories into the container under `/mnt/databases` and `mnt/blastdb` when running Snakemake.


### Step 5 – Prepare the runtime config and input files

1. Copy the template and adjust any parameters you need:

```bash
cp config/config.yaml config/config_runtime.yaml
```
Note: the template is usable as it is on gemini (as long as nobody else is running heavy job at the same time)
(See recommendations in `README.md > ### Step 5 – Prepare the runtime config`)

2. Put your input fasta(s) in the `data/` folder (extension must be .fasta).


### Step 6 – Run CoMR from the container image

#### Docker

1. The docker container is already available on gemini for all users.

2. Save the snippet below as `run_comr_gemini.sh`, edit the variables once, and
   launch it (`bash run_comr_gemini.sh`). Wrap the command in `nohup … &` or use
   `screen/tmux` only if you need it to survive logout.

```bash
#!/bin/bash
set -euo pipefail

# Declare local paths once so the bind list stays readable
# Replace "your_username" placeholders and FASTA name with your actual values

COMR_ROOT=/media/local_3/users/your_username/CoMR
DB_DIR=/media/local_4/CoMR_DB_hmm
BLAST_DB=/media/local_4/blastdb
TAXONOMY=/media/local_4/taxonomy
TARGETP=/home/gemini/software/targetp-2.0
CORES=24
FASTA_BASENAME=fasta_basename        # .fasta under data/
IMAGE=comr:latest
LOG=$COMR_ROOT/logs_${FASTA_BASENAME}/${FASTA_BASENAME}.log

mkdir -p "$(dirname "$LOG")"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --group-add 1006 \ # this allows docker to access the databases folder
  -e HOME=/opt/CoMR \
  -v "$DB_DIR:/mnt/databases:ro" \
  -v "$BLAST_DB:/mnt/blastdb:ro" \
  -v "$TAXONOMY:/mnt/taxonomy:ro" \
  -v "$TARGETP:/mnt/software/targetp-2.0:ro" \
  -v "$COMR_ROOT:/opt/CoMR" \
  -w /opt/CoMR \
  "$IMAGE" \
  snakemake --cores "$CORES" \
    --configfile config/config_runtime.yaml \
    --config fasta="$FASTA_BASENAME" enable_customdb=true 2>&1 | tee "$LOG"
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
