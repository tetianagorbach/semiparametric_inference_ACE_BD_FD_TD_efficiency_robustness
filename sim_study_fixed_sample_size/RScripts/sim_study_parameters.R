# this file  defines the parameters of the simulation study
seed <- 11771177
number.of.clusters <- 10 # numebr of clusters for parallel computing

number.of.replicates <- 10 # number of replicates
sample.size <- 100

parameters <- expand.grid(alpha1 = c(-2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2),
                          beta1 = c(-1.5, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 1.5), 
                          gamma1 = c(-2,  -1, -0.5, -0.1, 0, 0.1, 0.5,  1, 2),
                          gamma2 =  c(-2,  -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2))
pc <- 0.5
sigma.y <- 1
sigma.z <- 1
