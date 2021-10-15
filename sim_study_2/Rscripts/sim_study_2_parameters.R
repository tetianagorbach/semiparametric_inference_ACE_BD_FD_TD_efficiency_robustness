# this file  defines the parameters of the simulation study
# seed <- 2164681
# number.of.clusters <- 10 # number of clusters for parallel computing
# 
# number.of.replicates <- 1000 # number of replicates
# sample.size <- 50000
# 
# 
# parameters <- data.frame(alpha = 1, 
#                          beta = 1.5,
#                          gamma1 = 1.5,
#                          gamma2 = 1.5
#                          )
# 
# pc <- 0.5
# sigma.y <- 1
# sigma.z <- 1
# 

# output.file.name <- gsub(":| |-", "_", paste0("sim_study_2/results/sim_study_misspecification_res",  Sys.time(), ".Rdata"))


seed <- 2164681
number.of.clusters <- 10 # number of clusters for parallel computing

number.of.replicates <- 10 # number of replicates
sample.size <- 100


parameters <- data.frame(alpha = 1, 
                         beta = 1.5,
                         gamma1 = 1.5,
                         gamma2 = 1.5
)

pc <- 0.5
sigma.y <- 1
sigma.z <- 1

output.file.name <- gsub(":| |-", "_", paste0("sim_study_2/results/sim_study_misspecification_res",  Sys.time(), ".Rdata"))