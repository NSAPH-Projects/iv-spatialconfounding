#!/bin/bash
#SBATCH -c 8 # Request 8 cores for parallel processing
#SBATCH -t 00-10:00:00 # Amount of time needed DD-HH:MM:SS
#SBATCH -p hsph # Partition to submit to
#SBATCH --mem=12000 # Memory
#SBATCH -o error.out #specify where to save errors returned by the program
#SBATCH -e log.err #specify where to save the output log
#SBATCH --mail-type=END #notifications for job done
#SBATCH --mail-user=swoodward@g.harvard.edu # send to address
#SBATCH -N 1 #Number of nodes

my_packages=${HOME}/R/ifxrstudio/RELEASE_3_18
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_18.sif"

singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} Rscript run_simfunc.R "$1" "$2" "$3" "$4" "$5"
