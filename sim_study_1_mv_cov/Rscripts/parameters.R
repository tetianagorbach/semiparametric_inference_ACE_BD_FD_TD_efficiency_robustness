# this file  defines the parameters of the simulation study
seed <- 12345
number.of.clusters <- 2 # number of clusters for parallel computing
number.of.replicates <- 1000 # number of replicates
sample.size <- c(100, 500, 1000, 5000, 10000, 20000, 30000, 40000, 50000)



parameters <- expand.grid(alpha = c(0,1), 
                          beta = c(0, 1.5), 
                          gamma1 = c(0.5, 1.5),
                          gamma2 = c(0.1, 1))
# parameters <- data.frame(alpha1 = , #A|C
#                          alpha2 = 0, 
#                          alpha3 = 0, 
#                          alpha4 = 0, 
#                          alpha5 = 0,
#                          beta1 = rep(c(0, 1.5), each = 8), #Z|A
#                          gamma1 = rep(rep(c(0.5, 1.5), each = 4), times = 2), #Y| Z, C : Z
#                          gamma2 = rep (rep(c(0.5, 1.5), each = 2), times = 4),#Y| Z, C : C
#                          gamma3 = rep(c(0.5, 1.5), times = 8),
#                          gamma4 = 1,
#                          gamma5 = 1,
#                          gamma6 = 1)


pc <- 0.5
sigma.y <- 1
sigma.z <- 1

output.file.name <- gsub(":| |-", "_", paste0("sim_study_1_mv_cov/results/res",  Sys.time(), ".Rdata"))

