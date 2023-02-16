# code for the simulation study
require(foreach)
require(doParallel)
require(doRNG) # for independent RNG streams.
require(MASS)

source("sim_study_1_mv_cov/Rscripts/if_bounds_utils_mv_cov.R") # load functions
source("sim_study_1_mv_cov/Rscripts/parameters.R") # load parameters


cl <- makeCluster(number.of.clusters)
registerDoParallel(cl)
set.seed(seed)

sim.results <- matrix(nrow = 0, ncol = 18)
colnames(sim.results) <-  c("dgp",  "n","alpha1",
                            "beta1", 
                            "gamma1", "gamma2",
                            "est.fd.semipar", "est.bd.semipar", "est.td.semipar", "est.eif.td.bd.semipar", "est.eif.fd.td.semipar","est.eif.all.semipar",
                            "est.var.fd.semipar", "est.var.bd.semipar", "est.var.td.semipar", 
                            "est.var.td.bd.semipar","est.var.fd.td.semipar","est.var.all.semipar"
                            )

for (i in 1:nrow(parameters)){
        
        
        print(i)
        pc <- pc
        par.a.c <-  as.numeric(c(0, rep(parameters$alpha[i],5)))
        par.z.a <-  as.numeric(c(0, parameters$beta[i]))
        par.y.zc <-  as.numeric(c(0, parameters$gamma1[i], rep(parameters$gamma2[i], 5)))
        
        # ate <- CalculateMeanPotentialOutcome(astar = 1, pc = pc, par.z.a = par.z.a, par.y.zc = par.y.zc) -
        #         CalculateMeanPotentialOutcome(astar = 0, pc = pc, par.z.a = par.z.a, par.y.zc = par.y.zc)
        
         for (n in sample.size){
                parOut <- foreach(s=1:number.of.replicates, .combine='rbind') %dorng% {
                        # Generate data from DGP i with sample size n
                        data <- GenerateData(n,
                                                           pc = pc,
                                                           par.a.c = par.a.c,
                                                           par.z.a = par.z.a,
                                                           par.y.zc = par.y.zc,
                                                           sigma.z = sigma.z,
                                                           sigma.y = sigma.y)
        
                        # Fit nuisance models
                        fit.z <- lm(z ~ a, data = data)
                        fit.y <- lm(y ~  z  + c1 + c2 + c31 + c32 + c33, data = data)
                        fit.a <- glm(a ~ c1 + c2 + c31 + c32 + c33, data = data, family = binomial)
                        fit.y.az <- lm(y ~ a + z, data = data)
                        fit.y.ac <- lm(y ~ a + c1 + c2 + c31 + c32 + c33, data = data)
        
                        # Semiparametric estimates
                        est.if.bd <- EstimateIFBackDoor(cov.vals.all = data[, c("c1", "c2","c31", "c32","c33")], 
                                                        exposure = data$a, 
                                                        outcome = data$y,
                                                        fit.a.c = fit.a, 
                                                        fit.z.a = fit.z, 
                                                        fit.y.ac = fit.y.ac, 
                                                        astar = 1) -
                                EstimateIFBackDoor(cov.vals.all = data[, c("c1", "c2","c31", "c32","c33")], 
                                                   exposure = data$a, 
                                                   outcome = data$y,
                                                   fit.a.c = fit.a, 
                                                   fit.z.a = fit.z, 
                                                   fit.y.ac = fit.y.ac, 
                                                   astar = 0)
                        
                        est.if.fd <- EstimateIFFrontDoor(exposure = data$a,
                                                         intermediate = data$z,
                                                         outcome = data$y, 
                                                         est.pa = mean(data$a),
                                                         fit.z = fit.z, 
                                                         fit.y.az = fit.y.az, 
                                                         astar = 1) -
                                EstimateIFFrontDoor(exposure = data$a,
                                                    intermediate = data$z,
                                                    outcome = data$y, 
                                                    est.pa = mean(data$a),
                                                    fit.z = fit.z, 
                                                    fit.y.az = fit.y.az,  
                                                    astar = 0)
                        
                        est.if.td <- EstimateIFTwoDoor(cov.vals.all= as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                       exposure = data$a, 
                                                       intermediate = data$z, 
                                                       outcome = data$y,
                                                       fit.a.c = fit.a, 
                                                       fit.z.a = fit.z, 
                                                       fit.y.zc = fit.y, 
                                                       astar = 1) -
                                EstimateIFTwoDoor(cov.vals.all = as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                  exposure = data$a,
                                                  intermediate = data$z,
                                                  outcome = data$y,
                                                  fit.a.c = fit.a,
                                                  fit.z.a = fit.z, 
                                                  fit.y.zc = fit.y, 
                                                  astar = 0)
                        
                        est.eif.bd.td <- EstimateIFBackTwoDoor(cov.vals.all= as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                               exposure = data$a, 
                                                               intermediate = data$z, 
                                                               outcome = data$y,
                                                               fit.a.c = fit.a,
                                                               fit.z.a = fit.z,
                                                               fit.y.zc = fit.y,
                                                               astar = 1) -
                                EstimateIFBackTwoDoor(cov.vals.all = as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                      exposure = data$a,
                                                      intermediate = data$z,
                                                      outcome = data$y,
                                                      fit.a.c = fit.a,
                                                      fit.z.a = fit.z, 
                                                      fit.y.zc = fit.y, 
                                                      astar = 0)
                        
                        
                        est.eif.fd.td <- EstimateIFFrontTwoDoor(cov.vals.all= as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                                   exposure = data$a, 
                                                                   intermediate = data$z, 
                                                                   outcome = data$y, 
                                                                   fit.a.c = fit.a, 
                                                                   fit.z.a = fit.z, 
                                                                   fit.y.zc = fit.y, 
                                                                   astar = 1) -
                                EstimateIFFrontTwoDoor(cov.vals.all = as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                       exposure = data$a,
                                                       intermediate = data$z,
                                                       outcome = data$y,
                                                       fit.a.c = fit.a, 
                                                       fit.z.a = fit.z, 
                                                       fit.y.zc = fit.y,
                                                       astar = 0)
                        est.eif.all <- EstimateIFBackFrontTwoDoor(cov.vals.all = as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                                  exposure = data$a, 
                                                                  intermediate = data$z, 
                                                                  outcome = data$y,
                                                                 fit.a.c = fit.a, 
                                                                 fit.z.a = fit.z, 
                                                                 fit.y.zc = fit.y,
                                                                 astar = 1) -
                                EstimateIFBackFrontTwoDoor(cov.vals.all = as.matrix(data[, c("c1", "c2","c31", "c32","c33")]), 
                                                           exposure = data$a,
                                                           intermediate = data$z,
                                                           outcome = data$y,
                                                           fit.a.c = fit.a, 
                                                           fit.z.a = fit.z, 
                                                           fit.y.zc = fit.y, 
                                                           astar = 0)
        
        
                        # Parametric estimates
                        est.bd.par <- as.numeric(fit.y.ac$coefficients["a"])
                        est.fd.par <- as.numeric(fit.y.az$coefficients["z"]) * as.numeric(fit.z$coefficients["a"]) * (1 - 0)
                        est.td.par <- as.numeric(fit.y$coefficients["z"]) * as.numeric(fit.z$coefficients["a"]) * (1 - 0)
                        est.naive <-  mean(data$y[data$a==1]) - mean(data$y[data$a==0])
        
                        
                        c(i, n, as.numeric(parameters[i, c("alpha")]),
                          parameters[i,"beta"],
                          as.numeric(parameters[i,c("gamma1", "gamma2")]),
                          sapply(list(est.if.fd, est.if.bd, est.if.td, est.eif.bd.td, est.eif.fd.td, est.eif.all), mean),
                          sapply(list(est.if.fd, est.if.bd, est.if.td, est.eif.bd.td, est.eif.fd.td, est.eif.all), var) )
                }
              
                sim.results <- rbind(sim.results, parOut)
        }
        save(sim.results, file = output.file.name)
}


rm(a, h, i, n, S, parOut, ate, sample.sizes, par.a.c, par.y.zc, par.z.a, sample.size,cl )
rm(list = lsf.str())
save(list = ls(), file = output.file.name)

stopCluster(cl)



