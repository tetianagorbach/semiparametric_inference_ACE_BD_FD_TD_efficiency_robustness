#!/bin/bash
### SNAC project number, enter your own
#SBATCH -A SNIC2020-5-672 
# Asking for one core
#SBATCH -n 1
#SBATCH --time=60:20:00

# Serial job
# No matter how many processers you request this job will run
# on _only_ one core.

# Load the module first.
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.0.0


R --no-save --no-restore -f sim_study_misspecification/sim_study_misspecification_hpc2n.R