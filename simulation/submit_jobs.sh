#!/bin/bash

# Define arrays for options and methods
confounding_mechanisms=(1 2 3 4 5 6)
options=("linear" "nonlinear")

# Loop through each combination and submit jobs
for cm in "${confounding_mechanisms[@]}"; do
  for option in "${options[@]}"; do
  sbatch --job-name="cm${cm}_${option}" \
      --output="output/cm${cm}_${option}.out" \
      --error="error/cm${cm}_${option}.err" \
      run_job.sh 1000 "$cm" "$option" "TRUE" "8"
  sleep 1 # pause to be kind to the scheduler
  done
done

