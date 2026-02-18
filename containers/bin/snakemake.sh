#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="/opt/micromamba/envs/comr/bin/python"
exec "${PYTHON_BIN}" -m snakemake "$@"
