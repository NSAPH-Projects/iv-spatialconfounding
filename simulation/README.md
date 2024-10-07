# Simulation

This folder contains the code for the simulation section of the paper. 
- `simulation_preprocessing.R` creates `sim.RData`, the ingredients of the simulation data, and plots one of the simulated datasets.
- `run_job.sh` and `submit_jobs.sh` are the job script for slurm and shell script for job submission on the cluster.
- `run_simfunc.R`, called on by `run_job.sh` runs the simulation with the given parameters.
- `analysis.R` creates plots of the estimated ERFs and calculates average absolute bias and average RMSE.

## Sources of Data 
- [ ] 1. [Census](https://www.census.gov/) TIGER Line Shapefiles (County,State)
