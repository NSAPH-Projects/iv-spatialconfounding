# Simulation

This folder contains the code for the simulation section of the paper. 

- `simulation_preprocessing.R` creates `sim.RData` (the ingredients of the simulation data) and plots one of the simulated datasets.
- `run_job.sh` and `submit_jobs.sh` are the job script for slurm and shell script for job submission on the cluster.
- `run_simfunc.R`, called on by `run_job.sh` runs the simulation with the given parameters.
- The results are saved in the `results_Oct1` folder as csvs.
- `analysis.R` creates plots of the estimated ERFs in `images` and calculates average absolute bias and average RMSE.

## To reproduce

For confounding mechanisms 1-2:

1. Run `simulation_preprocessing.R` to create `sim.RData` and plot one of the simulated datasets.
2. Submit `submit_jobs.sh` to the cluster. This will run `run_simfunc.R` for each combination of parameters and save results as csvs in the `results_Oct1` folder. 
3. Run `analysis.R` to create plots of the estimated ERFs and print tables of average absolute bias and average RMSE. The plots will be saved in the `images` folder.

For confounding mechanism 3:

1. Change `within_state_GP = F` to `within_state_GP = T` in `run_simfunc.R`.
2. Repeat steps 2-3 above.

### Sources of Data 

1. [Census](https://www.census.gov/) TIGER Line 2010 Shapefiles (County,State).
