#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: mitofates.sh <fasta> <organism>" >&2
  exit 1
fi
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PERL5LIB="/opt/software/perl5/lib/perl5:/usr/local/share/perl/5.28.1:/usr/local/lib/perl/5.28.1:/usr/local/share/perl5:/usr/local/lib/perl5:/usr/local/lib/x86_64-linux-gnu/perl/5.28.1"
exec /usr/bin/perl /opt/software/MitoFates/MitoFates.pl "$@"
