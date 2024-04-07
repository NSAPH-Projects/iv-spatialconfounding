#!/bin/bash
#SBATCH -c 8 # Number of threads
#SBATCH -t 02-23:59:00 # Amount of time needed DD-HH:MM:SS
#SBATCH -p shared # Partition to submit to
#SBATCH --mem=10000 #Memory
#SBATCH -o /n/home07/swoodward/simulation/error.out #specify where to save errors returned by the program
#SBATCH -e /n/home07/swoodward/simulation/log.err #specify where to save the output log
#SBATCH --mail-type=END #notifications for job done
#SBATCH --mail-user=swoodward@g.harvard.edu # send to address
#SBATCH -N 1 #Number of nodes

module load R/4.3.3-fasrc01
export R_LIBS_USER=/n/home07/swoodward/apps/R_4.3.3:$R_LIBS_USER  #change this accordingly

# Run the R script
#singularity exec sing_biocond_3.14.sif 
Rscript run_simfunc.R "$1" "$2" "$3" "$4" "$5"