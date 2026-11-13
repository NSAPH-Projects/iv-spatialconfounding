# Simulation

Code for the simulation study (Section 5 and Supplementary Section 4): 5
confounding mechanisms × 2 outcome models, 1000 datasets of n = 1000 each. The
data-generating process is described at the top of `funcs.R`.

Run every command from `simulation/`.

| Step | Command | Output |
|---|---|---|
| **1. Fixed inputs** (spatial surfaces) | `Rscript preprocessing.R` | `sim.rds` |
| **2. Simulation** | `bash submit_jobs.sh final 1000 linear,tanh 0.5 noniv,ivtps,ivgl LS,MS,twoU,coordnoise,stateJump 250` | `results_final/` |
| **3. Supplement tables** | `Rscript analysis.R` | `tables/results_final.tex` |
| **4. Main-text table** | `Rscript make_main_table.R` | `tables/main_rmse_final.tex` |
| **5. Maps of the confounding mechanisms** | `Rscript make_decomp.R` | `images/all-decomp.jpeg` |
| **6. Sensitivity analysis** | `sbatch run_sensitivity.sh` | `images/sensitivity.rds`, `images/sensitivity.png` |
| **7. Sensitivity and selected k** | `Rscript plot_combined.R`; main-text version: `Rscript plot_combined.R images/sensitivity.rds results_final images/sensitivity_selection_main.png LS,MS linear` | `images/sensitivity_selection.png`, `images/sensitivity_selection_main.png` |
