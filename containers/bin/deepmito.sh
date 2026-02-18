#!/usr/bin/env bash
set -euo pipefail
DEEPMITO_ENV="/opt/conda"
python_bin="${DEEPMITO_ENV}/bin/python3"
if [[ ! -x "$python_bin" ]]; then
  echo "DeepMito python interpreter not found at ${python_bin}" >&2
  exit 1
fi
export PATH="$DEEPMITO_ENV/bin:$PATH"
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH="$DEEPMITO_ENV/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$DEEPMITO_ENV"
exec "$python_bin" /usr/src/deepmito/deepmito.py "$@"
