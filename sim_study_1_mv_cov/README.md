Reports contains an R markdown document to construct the Figure and the Tables that represent the results of the simulation study in the paper.

Results contains the results of the simulations after running the simulations on High Performance Computing Center North (HPC2N).

Rscripts contains the following files:

        - if_bounds_utils.R provides R functions for the simulation study
        
        - job.sh provides the parameters and the code used to run the simulation study on the HPC2N. To run the code on HPC2N, first one has to log in to HPC2N and then run the command  sbatch job.sh in Terminal (job.sh should be moved to the root directory).
        
        - run_sim_study_1_mv_cov.R provides the code to run the simulation study
        
        - parameters.R initialises the parameters of the simulation study.




