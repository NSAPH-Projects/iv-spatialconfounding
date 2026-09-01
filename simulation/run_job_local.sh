#!/usr/bin/env bash
set -euo pipefail

if (( $# < 3 || $# > 10 )); then
  echo "Usage: $0 <nsims> <mechanism> <linear|nonlinear> [select_basis] [n_cores] [results_dir] [seed] [methods] [core_fraction] [alpha]" >&2
  exit 2
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

exec Rscript run_simfunc.R "$@"
