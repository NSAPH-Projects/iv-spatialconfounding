#!/bin/bash
#SBATCH -c 1 # Number of threads, just 1, not parallelizing
#SBATCH -t 00-03:00:00 # Amount of time needed DD-HH:MM:SS
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=650 # Memory, uses about 550
#SBATCH -o /n/home07/swoodward/instrumental_variables_simulation/error.out #specify where to save errors returned by the program
#SBATCH -e /n/home07/swoodward/instrumental_variables_simulation/log.err #specify where to save the output log
#SBATCH --mail-type=END #notifications for job done
#SBATCH --mail-user=swoodward@g.harvard.edu # send to address
#SBATCH -N 1 #Number of nodes

my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} Rscript run_simfunc.R "$1" "$2" "$3"
