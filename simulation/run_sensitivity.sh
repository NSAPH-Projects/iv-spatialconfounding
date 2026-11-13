#!/bin/bash
#SBATCH -c 8 # Cores; trajectories run in parallel
#SBATCH -t 06:00:00 # 5 mechanisms take ~1-2 h on 8 cores
#SBATCH -p hsph # Partition to submit to
#SBATCH --mem=8000 # Memory for all cores (peaked at ~3.3 GB)
#SBATCH -o sensitivity.out # Output log (in the submission directory)
#SBATCH -e sensitivity.err # Error log
#SBATCH --mail-type=END,FAIL #notifications for job done
#SBATCH --mail-user=swoodward@g.harvard.edu # send to address
#SBATCH -N 1 #Number of nodes

# Usage (from simulation/):  sbatch run_sensitivity.sh [key=value ...]
#   e.g. sbatch run_sensitivity.sh reps=5 scenarios=LS,MS
# Output: images/sensitivity.rds (trajectories), images/sensitivity.png,
# images/sensitivity_correlation.png.

my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} \
  Rscript sensitivity.R cores="${SLURM_CPUS_PER_TASK:-8}" "$@"
