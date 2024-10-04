#!/bin/bash

# Define arrays for options and methods
projmatnames=("GFT150" "GFT51100" "nestedstate" "nestedregion") # Adjust as needed
options=("linear" "interaction" "nonlinear")
methods=("KennedyERC" "spatialcoord" "spatialplus" "bobb_exposurepenalized" "keller_szpiro_selectingscale_outcome" "unadjustedOLS" "oracleU")  # "spectral_discrete" 

# Loop through each combination and submit spatial linear and oracle jobs
for projmatname in "${projmatnames[@]}"; do
  for option in "${options[@]}"; do
    for method in "${methods[@]}"; do
  sbatch --job-name=${projmatname}_${option}_${method}\
    --output=output/${projmatname}_${option}_${method}.out \
    --error=error/${projmatname}_${option}_${method}.err \
    run_job.sh 100 "$projmatname" "$option" "$method" 
  sleep 1 # pause to be kind to the scheduler
    done
  done
done

