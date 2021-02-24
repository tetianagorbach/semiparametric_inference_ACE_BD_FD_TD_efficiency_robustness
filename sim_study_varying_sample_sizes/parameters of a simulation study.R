# this file provides the values of the parameters for the simulation study
seed <- 100
number.of.clusters <- 10 # numebr of clusters for parallel computing

number.replicates <- 100 # number of replicates
sample.size <- c(20, 50, 80, 100, 200, 400, 800, 1000, 1500, 2000, 2500, 5000, 10000)

parameters <- list(
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,0.1), par.y.zc=c(1, 2, 0.1),
             sigma.z = 1,sigma.y = 1),# strong relation of Z on Y
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,0.1), par.y.zc=c(1, 0.1, 2),
             sigma.z = 1,sigma.y = 1), # weak relation of Z on Y, but strong C on Y
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,0.1), par.y.zc=c(1, 1, 1),
             sigma.z = 1,sigma.y = 1), # similar effect of C and Z on Y
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,1), par.y.zc=c(1, 2, 0.1),
             sigma.z = 1,sigma.y = 1), 
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,1), par.y.zc=c(1, 0.1, 2),
             sigma.z = 1,sigma.y = 1), 
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,1), par.y.zc=c(1, 1, 1),
             sigma.z = 1,sigma.y = 1), 
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,1.5), par.y.zc=c(1, 2, 0.1),
             sigma.z = 1,sigma.y = 1), 
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,1.5), par.y.zc=c(1, 0.1, 2),
             sigma.z = 1,sigma.y = 1), 
        list(pc=0.5, par.a.c = c(0,0.5), par.z.a = c(0,1.5), par.y.zc=c(1, 1, 1), 
             sigma.z = 1,sigma.y = 1)
)