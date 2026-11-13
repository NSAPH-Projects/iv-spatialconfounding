#!/bin/bash
#SBATCH -c 4 # Cores; submit_jobs.sh overrides -c, --mem, and -t per method group
#SBATCH -t 02:00:00 # Amount of time needed DD-HH:MM:SS
#SBATCH -p hsph # Partition to submit to
#SBATCH --mem=4000 # Memory for all cores
#SBATCH --mail-type=FAIL #notifications for job done
#SBATCH --mail-user=swoodward@g.harvard.edu # send to address
#SBATCH -N 1 #Number of nodes

# Usage: called by submit_jobs.sh; passes key=value arguments to run_sim.R.

my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} \
  Rscript run_sim.R cores="${SLURM_CPUS_PER_TASK:-4}" "$@"
