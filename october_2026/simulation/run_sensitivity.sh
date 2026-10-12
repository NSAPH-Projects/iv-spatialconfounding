#!/bin/bash
#SBATCH -c 8 # Request 8 cores for parallel processing
#SBATCH -t 01-00:00:00 # Amount of time needed DD-HH:MM:SS
#SBATCH -p hsph # Partition to submit to
#SBATCH --mem=48000 # Memory for all cores (each mclapply worker is a full R fork)
#SBATCH -o /n/home07/swoodward/instrumental_variables_simulation_sep2026/october_2026/simulation/sensitivity.out #specify where to save the output log
#SBATCH -e /n/home07/swoodward/instrumental_variables_simulation_sep2026/october_2026/simulation/sensitivity.err #specify where to save errors returned by the program
#SBATCH --mail-type=END #notifications for job done
#SBATCH --mail-user=swoodward@g.harvard.edu # send to address
#SBATCH -N 1 #Number of nodes

# Usage (from october_2026/simulation/):  sbatch run_sensitivity.sh
# Output: images/sensitivity_basis_dimension.png (+ .rds with the fits).

my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} Rscript sensitivity_basis_dimension.R
