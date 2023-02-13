# this file  defines the parameters of the simulation study
seed <- 0
number.of.clusters <- 10 # number of clusters for parallel computing
number.of.replicates <- 10 # number of replicates
sample.size <- c(50, 100, 500, 1000)#, 5000, 10000, 20000, 30000, 40000, 50000)


parameters <- data.frame(alpha1 = 1,
                         alpha2 = 1, 
                         beta1 = 0.5,
                         beta2 = 1,
                         beta3 = 1.5,
                         gamma1 = 1.5,
                         gamma2 = 2,
                         gamma3 = 2.5)
pc <- 0.5
sigma.y <- 1
sigma.z <- 1

output.file.name <- gsub(":| |-", "_", paste0("sim_study_1_complex_cov/results/res",  Sys.time(), ".Rdata"))

