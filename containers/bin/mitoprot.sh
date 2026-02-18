#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: mitoprot.sh <fasta> <output_dir>" >&2
  exit 1
fi
input_fasta="$1"
output_dir="$2"
mkdir -p "$output_dir"
input_abs="$(readlink -f "$input_fasta")"
output_abs="$(readlink -f "$output_dir")"
base_name="$(basename "$input_fasta")"
sequence_prefix="${base_name%%.*}"
output_file="${output_abs}/${sequence_prefix}.mitoprot"
declare -a CANDIDATE_BINS=(
  "/opt/software/mitoprotII/mitoprot"
  "/opt/software/mitoprotII/mitoprot.lin"
  "/opt/software/mitoprotII/mitoprot.linux"
  "/opt/software/mitoprotII/mitoprot"
)
mitoprot_bin=""
for candidate in "${CANDIDATE_BINS[@]}"; do
  if [[ -x "$candidate" ]]; then
    mitoprot_bin="$candidate"
    break
  fi
done
if [[ -z "$mitoprot_bin" ]]; then
  mitoprot_bin=$(find /opt/software/mitoprotII -maxdepth 2 -type f -executable -name 'mitoprot*' ! -name '*.adb' | head -n 1 || true)
fi
if [[ -z "$mitoprot_bin" ]]; then
  echo "Unable to locate a mitoprot executable under /opt/software/mitoprotII" >&2
  exit 1
fi
(
  cd "$output_abs"
  "$mitoprot_bin" -f a "$input_abs" > "$output_file"
)
echo "Complet results are writed in ${sequence_prefix}.mitoprot"
echo "MitoProt analysis complete. Results saved in $output_dir"
