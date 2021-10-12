# this file  defines the parameters of the simulation study
seed <- 0
number.of.clusters <- 10 # number of clusters for parallel computing
number.of.replicates <- 1000 # number of replicates
sample.size <- c(50, 100, 500, 1000, 5000, 10000, 20000, 30000, 40000, 50000)

parameters <- data.frame(alpha1 = 1, 
                         beta1 = rep(c(0.5, 1.5), each = 4),
                         gamma1 = rep(c(0.5, 0.5, 1.5, 1.5), times = 2),
                         gamma2 = rep (c(0.5, 1.5), times = 4)
                         )
pc <- 0.5
sigma.y <- 1
sigma.z <- 1

output.file.name <- gsub(":| |-", "_", paste0("sim_study_1/results/results_sim_varied_sample_size_8dgp_nrep1000_n50000",  Sys.time(), ".Rdata"))