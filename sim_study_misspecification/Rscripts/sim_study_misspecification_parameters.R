# this file  defines the parameters of the simulation study
seed <- 0
number.of.clusters <- 10 # number of clusters for parallel computing

number.of.replicates <- 1000 # number of replicates
sample.size <- 1000


parameters <- data.frame(alpha = 1, 
                         beta = c(0.5),
                         gamma1 = 1.5,
                         gamma2 = 1.5
                         )

pc <- 0.5
sigma.y <- 1
sigma.z <- 1

output.file.name <- gsub(":| |-", "_", paste0("sim_study_misspecification/results/sim_study_misspecification",  Sys.time(), ".Rdata"))