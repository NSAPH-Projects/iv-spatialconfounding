#!/bin/bash
set -euo pipefail

# Usage: bash submit_jobs.sh [family_alpha] [run_label] [nsims] [time_limit]
#   family_alpha: Bonferroni-corrected significance level for the basis-
#                 selection departure test (must be in (0, 1); default 0.05)
#   run_label:    used to build a distinct results/output/error directory so
#                 different family_alpha runs never share a results_dir
#                 (default derived from family_alpha, e.g. "alpha0.05")
#   nsims:        replicates per (mechanism, option) cell (default 1000)
#   time_limit:   SLURM --time for each simulation job (default 1-00:00:00)
#
# Pilot:     bash submit_jobs.sh 0.05 pilot 20 02:00:00
# Full run:  bash submit_jobs.sh 0.05
family_alpha="${1:-0.05}"
run_label="${2:-alpha${family_alpha}}"
nsims="${3:-1000}"
time_limit="${4:-1-00:00:00}"

base_dir="/n/home07/swoodward/instrumental_variables_simulation_sep2026"
sim_subdir="october_2026/simulation"
output_subdir="output_${run_label}"
error_subdir="error_${run_label}"
results_subdir="results_stability_${nsims}reps_quadratic_${run_label}"

mkdir -p "${base_dir}/${sim_subdir}/${output_subdir}" "${base_dir}/${sim_subdir}/${error_subdir}" \
  "${base_dir}/${sim_subdir}/${results_subdir}"

# Mechanisms 4 and 8 (scale-reversed) are on hold: with the mismatched basis
# their walk can run hundreds of columns past the boundary, each step
# refactorizing a ~900-column basis. Add them back here when ready.
confounding_mechanisms=(1 2 3 5 6 7)
options=("linear" "nonlinear")
n_cores=8
seed=20260825
r2_reference=0.7
# +spatialcoord IV variants and trueIV-spatialcoord dropped from the study;
# plain spatialcoord is kept. All 6 methods run for every family_alpha value
# (each family_alpha's results directory is fully self-contained), even
# though family_alpha only actually affects IV-TPS/IV-GraphLaplacian --
# oracle/baseline/spatialcoord/trueIV results will be identical (same seed)
# across family_alpha runs.
# Each cell is split into three jobs so no single job runs both slow IV
# methods back to back. Every method writes its own output files, so the
# jobs never touch the same file.
method_groups=("oracle,baseline,spatialcoord,trueIV" "IV-TPS" "IV-GraphLaplacian")
# New directory per run -- funcs.R appends new reps as extra
# columns onto existing output files if the file already exists. Mixing reps
# from different DGP/estimator versions into the same files would silently
# corrupt results.
results_dir="${base_dir}/${sim_subdir}/${results_subdir}/"

echo "Submitting nsims=${nsims}, family_alpha=${family_alpha} (run_label=${run_label}), time limit ${time_limit}, into ${results_dir}"

# Loop through each combination and submit jobs
for cm in "${confounding_mechanisms[@]}"; do
  for option in "${options[@]}"; do
    for methods in "${method_groups[@]}"; do
      group_tag="${methods%%,*}"  # oracle | IV-TPS | IV-GraphLaplacian
      sbatch --job-name="cm${cm}_${option}_${group_tag}_${run_label}" \
          --time="$time_limit" \
          --output="${base_dir}/${sim_subdir}/${output_subdir}/cm${cm}_${option}_${group_tag}.out" \
          --error="${base_dir}/${sim_subdir}/${error_subdir}/cm${cm}_${option}_${group_tag}.err" \
          run_job.sh "$nsims" "$cm" "$option" "TRUE" "$n_cores" \
          "$results_dir" "$seed" "$methods" "$r2_reference" "$family_alpha"
      sleep 1 # pause to be kind to the scheduler
    done
  done
done

# Per-replicate coverage truth (one job per mechanism -- shared across both
# options, see compute_true_psi.R). Submitted concurrently with the main
# jobs above, not dependent on them: it's a pure re-derivation of the same
# confounder field from the same seed/nsims (reusing those exact variables,
# so it can never drift out of sync with the main run), not a function of
# the main jobs' output, so there's no need to wait for them to finish.
for cm in "${confounding_mechanisms[@]}"; do
  sbatch --job-name="truth_cm${cm}_${run_label}" \
      --output="${base_dir}/${sim_subdir}/${output_subdir}/truth_cm${cm}.out" \
      --error="${base_dir}/${sim_subdir}/${error_subdir}/truth_cm${cm}.err" \
      run_job_truth.sh "$cm" "$nsims" "$seed" "$results_dir" "$n_cores"
  sleep 1 # pause to be kind to the scheduler
done
