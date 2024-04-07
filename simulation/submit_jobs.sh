#!/bin/bash

# Define arrays for options and methods
projmatnames=("GFT150" "GFT51100" "nestedstate") # Adjust as needed
options=("linear" "interaction" "nonlinear")
methods=("KennedyERC" "spatialplus" "bobb_exposurepenalized" "spectral_discrete" "keller_szpiro_selectingscale_outcome" "unadjustedOLS" "GPCERF" "GPCERF_nn")

# Loop through each combination and submit spatial linear and oracle jobs
for projmatname in "${projmatnames[@]}"; do
  for option in "${options[@]}"; do
    for method in "${methods[@]}"; do
  sbatch --job-name=${projmatname}_${option}_${method}\
          --output=output/${projmatname}_${option}_${method}.out \
          --error=error/${projmatname}_${option}_${method}.err \
          run_job.sh 1000 "$projmatname" "$option" "$method" "no"
  sleep 1 # pause to be kind to the scheduler
    done
  done
done

sensitivitymethods=('KennedyERC', 'GPCERF', 'GPCERF_nn')
for projmatname in "${projmatnames[@]}"; do
  for option in "${options[@]}"; do
    for method in "${sensitivitymethods[@]}"; do
  sbatch --job-name=${projmatname}_${option}_${method}\
          --output=output/${projmatname}_${option}_${method}.out \
          --error=error/${projmatname}_${option}_${method}.err \
          run_job.sh 1000 "$projmatname" "$option" "$method" "yes"
  sleep 1 # pause to be kind to the scheduler
    done
  done
done