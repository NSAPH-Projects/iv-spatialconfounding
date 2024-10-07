#!/bin/bash

# Define arrays for options and methods
rangeus=("smallscale" "tinyscale") # Adjust as needed
options=("linear" "nonlinear")

# Loop through each combination and submit jobs
for rangeu in "${rangeus[@]}"; do
  for option in "${options[@]}"; do
  sbatch --job-name=${rangeu}_${option}\
    --output=output/${rangeu}_${option}.out \
    --error=error/${rangeu}_${option}.err \
    run_job.sh 100 "$rangeu" "$option"  
  sleep 1 # pause to be kind to the scheduler
  done
done

