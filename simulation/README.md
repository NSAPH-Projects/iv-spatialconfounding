# Simulation

This folder contains the code for the simulation section of the paper. 

- `simulation_preprocessing.R` creates `sim.RData` (input data to simulation) and plots each one of the simulated datasets.
- `run_job.sh` (slurm job script) and `submit_jobs.sh` (submission shell script) run the simulation on the cluster. `submit_jobs.sh` takes an `alpha` value as its first argument (the tuning parameter for the sequential stability criterion, default 1). One call submits a job for every (confounding mechanism 1-8) x (linear, nonlinear) combination for each method.
- `run_simfunc.R`, called by `run_job.sh`, runs the simulation for one (mechanism, option) with the given parameters.
- Each `submit_jobs.sh <alpha>` run saves its results as csvs to its own self-contained directory `results_manuscript_cluster_1000reps_quadratic_alpha<alpha>/`.
- `analysis.R` reads one results directory (set `results_dir` and `results_label` at the top of the script), saves boxplots of estimates to `images`, and prints the tables in the manuscript.

## To reproduce

1. Run `simulation_preprocessing.R` to create `sim.RData` and plot one of the simulated datasets.
2. For each alpha in 0.5, 1, 1.5, 2, submit `bash submit_jobs.sh <alpha>` to the cluster. Each call runs `run_simfunc.R` for every (mechanism 1-8) x (linear, nonlinear) combination and saves csvs to `results_manuscript_cluster_1000reps_quadratic_alpha<alpha>/`.
3. For each alpha, set `results_dir` / `results_label` at the top of `analysis.R` to that alpha's directory and run it to regenerate that alpha's figures (in the `images` folder) and print the table LaTeX. The manuscript reports these results across the four alpha values.

### Sources of Data 

1. [Census](https://www.census.gov/) TIGER Line 2010 Shapefiles (County,State).
