#!/bin/bash
set -euo pipefail

# Usage (from simulation/ on the cluster):
#   bash submit_jobs.sh [run_label] [nsims] [options] [tau] [groups] [scenarios] [chunk]
#   run_label: names the results/output/error directories (default "final")
#   nsims:     replicates per scenario x outcome cell (default 100)
#   options:   comma-separated outcome models (default linear,tanh)
#   tau:       exposure R^2 reference for IV-TPS / IV-GL (default 0.5)
#   groups:    comma-separated method groups to submit (default
#              noniv,ivtps,ivgl); only the IV groups depend on tau
#   scenarios: comma-separated (default LS,MS,twoU,coordnoise,stateJump)
#   chunk:     replicates per job (default nsims, i.e. one job per cell). With
#              chunk < nsims, each cell is split into jobs over replicates
#              1..chunk, chunk+1..2*chunk, ...; each writes
#              <scenario>_<option>_<group>_r<first rep>.csv, and the datasets
#              are the same as in an unchunked run.
# Each scenario x outcome cell is split into three jobs, so the slow IV
# methods never hold up the others. Every job writes its own file
# (<scenario>_<option>_<group>.csv); datasets depend only on the replicate
# number, so the files line up. Requires sim.rds (built locally by
# preprocessing.R) in this directory.
#
# Paper run (1000 reps in chunks of 250):
#   bash submit_jobs.sh final 1000 linear,tanh 0.5 noniv,ivtps,ivgl LS,MS,twoU,coordnoise,stateJump 250
# Quick test of one mechanism: bash submit_jobs.sh test 10 linear 0.5 noniv,ivtps,ivgl coordnoise
run_label="${1:-final}"
nsims="${2:-100}"
IFS=',' read -r -a options <<< "${3:-linear,tanh}"
tau="${4:-0.5}"
IFS=',' read -r -a wanted_groups <<< "${5:-noniv,ivtps,ivgl}"
IFS=',' read -r -a scenarios <<< "${6:-LS,MS,twoU,coordnoise,stateJump}"
chunk="${7:-$nsims}"

# group name | methods | cores | memory (MB) | time limit
# noniv jobs peak at ~1.6 GB; IV jobs at ~3.9 GB with 8 cores.
groups=(
  "noniv|oracle,baseline,spatialcoord,trueIV|4|2500|02:00:00"
  "ivtps|IV-TPS|8|6000|12:00:00"
  "ivgl|IV-GL|8|6000|12:00:00"
)

sim_dir="$(pwd)"
results_dir="${sim_dir}/results_${run_label}"
output_dir="${sim_dir}/output_${run_label}"
error_dir="${sim_dir}/error_${run_label}"

if [[ ! -f "${sim_dir}/sim.rds" ]]; then
  echo "sim.rds not found in ${sim_dir}; copy it from your local preprocessing run." >&2
  exit 1
fi
mkdir -p "$results_dir" "$output_dir" "$error_dir"

echo "Submitting ${#scenarios[@]} scenarios x ${#options[@]} options x groups ${wanted_groups[*]}, nsims=${nsims} (chunks of ${chunk}), tau=${tau}, into ${results_dir}"

for scenario in "${scenarios[@]}"; do
  for option in "${options[@]}"; do
    for group in "${groups[@]}"; do
      IFS='|' read -r tag methods cores mem time_limit <<< "$group"
      [[ " ${wanted_groups[*]} " == *" ${tag} "* ]] || continue
      for (( start = 1; start <= nsims; start += chunk )); do
        n=$(( nsims - start + 1 < chunk ? nsims - start + 1 : chunk ))
        file_tag="$tag"
        (( chunk < nsims )) && file_tag="${tag}_r${start}"
        sbatch --job-name="dag_${scenario}_${option}_${file_tag}_${run_label}" \
            -c "$cores" --mem="$mem" -t "$time_limit" \
            --output="${output_dir}/${scenario}_${option}_${file_tag}.out" \
            --error="${error_dir}/${scenario}_${option}_${file_tag}.err" \
            run_job.sh nsims="$n" rep_start="$start" scenarios="$scenario" options="$option" \
            methods="$methods" results="$results_dir" tag="$file_tag" tau="$tau"
        sleep 1 # pause to be kind to the scheduler
      done
    done
  done
done
