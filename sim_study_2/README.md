Reports contains an R markdown document to construct the Figure and the Table that represent the results of the simulation study in the paper.

Results contains the results of the simulations after running the simulations on High Performance Computing Center North (HPC2N).

Rscripts contains the following files:
        
        - job_2.sh provides the parameters and the code used to run the simulation study on HPC2N. To run the code on HPC2N, first one has to log in to HPC2N and then run the command  sbatch job_2.sh in Terminal (job_2.sh should be moved to the root directory).
        
        - run_sim_study_2.R provides the code to run the simulation study
        
        - sim_study_2_parameters.R initialises the parameters of the simulation study.



