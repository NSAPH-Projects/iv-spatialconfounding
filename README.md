# Harnessing the Scale of Spatial Confounding for Causal Inference

The folder simulation contains the code to run the simulations in Section 4 of the paper. The workflow is:
- 'Preprocessing.R' generates the graph, spatial coordinates, and other data used for the simulation.
- 'simfuncs.R' contains all utility functions for the simulation. The main function is 'simfunc' which will write results of simulations to csv files in a folder called 'results' (hidden).
- 'run_job.sh' and 'submit_jobs.sh' are the job script for slurm and shell script for job submission on the cluster.

## Sources of Data 
- [ ] 1. [Census](https://www.census.gov/) TIGER Line Shapefiles (County,State)

The computations in this paper were run on the FASRC cluster supported by the FAS Division of Science Research Computing Group at Harvard University.
