# code for the simulation study
require(foreach)
require(doParallel)
require(doRNG) # for independent RNG streams.


source("if_bounds_utils.R") # load functions
source("sim_study_parameters.R") # load parameters

cl <- makeCluster(number.of.clusters) 
registerDoParallel(cl)
set.seed(seed)

sim.results <- matrix(nrow = 0, ncol=29)
colnames(sim.results) <-  c("dgp", "ate", "n","alpha1", "beta1", "gamma1", "gamma2",
                            "est.fd.semipar", "est.bd.semipar", "est.td.semipar", "est.eif.td.bd.semipar", "est.eif.fd.td.semipar","est.eif.all.semipar",
                            "est.fd.par", "est.bd.par", "est.td.par", "est.naive",
                            "est.var.fd.semipar", "est.var.bd.semipar", "est.var.td.semipar", "est.var.td.bd.semipar","est.var.fd.td.semipar","est.var.all.semipar",
                            "bound.fd", "bound.bd", "bound.td", "bound.td.bd", "bound.fd.td","bound.fd.td.bd")


for (i in 1:nrow(parameters)){
        print(i)
        pc <- pc
        par.a.c <-  as.numeric(c(0, parameters$alpha1[i]))
        par.z.a <-  as.numeric(c(0, parameters$beta1[i]))
        par.y.zc <-  as.numeric(c(0, parameters$gamma1[i], parameters$gamma2[i]))
        
        ate <- psi_true_continuous(astar = 1, pc = pc, par.z.a = par.z.a, par.y.zc = par.y.zc)-
                psi_true_continuous(astar = 0, pc = pc, par.z.a = par.z.a, par.y.zc = par.y.zc)
        #### Bounds ###
        bound.bd.value <- BoundBackDoor(astar = 1, a = 0, pc = pc,
                                        par.a.c = par.a.c,
                                        par.y.zc = par.y.zc,
                                        sigma.y = sigma.y,
                                        sigma.z = sigma.z)
        
        bound.td.value <- BoundTwoDoor(astar = 1, a = 0, pc = pc,
                                       par.a.c = par.a.c,
                                       par.z.a  = par.z.a,
                                       par.y.zc = par.y.zc,
                                       sigma.y = sigma.y,
                                       sigma.z = sigma.z)
        
        bound.fd.value <- BoundFrontDoor(astar = 1, a = 0, pc = pc,
                                         par.a.c = par.a.c,
                                         par.z.a = par.z.a,
                                         par.y.zc = par.y.zc,
                                         sigma.y = sigma.y,
                                         sigma.z = sigma.z)
        
        
        bound.eif.td.bd.value <- BoundEIFTwoBackDoor(astar = 1, a = 0, pc = pc,
                                         par.a.c = par.a.c,
                                         par.z.a = par.z.a,
                                         par.y.zc = par.y.zc,
                                         sigma.y = sigma.y,
                                         sigma.z = sigma.z)
        
        bound.eif.fd.td.value <- BoundEIFFrontTwoDoor(astar = 1, a = 0, pc = pc,
                                                     par.a.c = par.a.c,
                                                     par.z.a = par.z.a,
                                                     par.y.zc = par.y.zc,
                                                     sigma.y = sigma.y,
                                                     sigma.z = sigma.z)
        
        bound.eif.fd.td.bd.value <- BoundEIFFrontTwoBackDoor(astar = 1, a = 0, pc = pc,
                                                             par.a.c = par.a.c,
                                                             par.z.a = par.z.a,
                                                             par.y.zc = par.y.zc,
                                                             sigma.y = sigma.y,
                                                             sigma.z = sigma.z)
        n <- sample.size
        parOut <- foreach(s=1:number.of.replicates, .combine='rbind') %dorng% {
                # Generate data from DGP i with sample size n
                data <- gen_data_ca_binary_zy_cont(n,
                                                   pc = pc,
                                                   par.a.c = par.a.c,
                                                   par.z.a = par.z.a,
                                                   par.y.zc = par.y.zc,
                                                   sigma.z = sigma.z,
                                                   sigma.y = sigma.y)

                # Fit nuisance models
                fit.z <- lm(z ~ a , data = data)
                fit.y <- lm(y ~  z  + c, data = data)
                fit.a <- glm(a ~ c, data = data, family = binomial)
                fit.y.az <- lm(y ~ a + z, data = data)
                fit.y.ac <- lm(y ~ a + c, data = data)

                # Semiparametric estimates
                est.if.fd <- EstimateInfluenceFunctionFrontDoor(exposure = data$a,intermediate = data$z,outcome = data$y,
                                                                fit.z = fit.z, fit.y.az = fit.y.az, astar = 1) -
                        EstimateInfluenceFunctionFrontDoor(exposure = data$a,intermediate = data$z,outcome = data$y,
                                                           fit.z = fit.z, fit.y.az = fit.y.az,  astar = 0)
                est.if.bd <- EstimateInfluenceFunctionBackDoor(cov.vals.all = data$c, exposure = data$a, outcome = data$y,
                                                               fit.a = fit.a, fit.z = fit.z, fit.y.ac = fit.y.ac, astar = 1) -
                        EstimateInfluenceFunctionBackDoor(cov.vals.all = data$c, exposure = data$a, outcome = data$y,
                                                          fit.a = fit.a, fit.z = fit.z, fit.y.ac = fit.y.ac, astar = 0)
                est.if.td <- EstimateInfluenceFunctionTwoDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                              fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                        EstimateInfluenceFunctionTwoDoor(cov.vals.all = data$c, exposure = data$a,intermediate = data$z,outcome = data$y,
                                                         fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)
                est.eif.td.bd <- EstimateEIFTwoBackDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                        fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                        EstimateEIFTwoBackDoor(cov.vals.all = data$c, exposure = data$a,intermediate = data$z,outcome = data$y,
                                               fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)
                est.eif.fd.td <- EstimateEIFFrontTwoDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                         fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                        EstimateEIFFrontTwoDoor(cov.vals.all = data$c, exposure = data$a,intermediate = data$z,outcome = data$y,
                                                fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)
                est.eif.all <- EstimateEIFAll(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                              fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                        EstimateEIFAll(cov.vals.all = data$c, exposure = data$a,intermediate = data$z,outcome = data$y,
                                       fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)


                # Parametric estimates
                est.bd.par <- as.numeric(fit.y.ac$coefficients["a"])
                est.fd.par <- as.numeric(fit.y.az$coefficients["z"]) * as.numeric(fit.z$coefficients["a"]) * (1 - 0)
                est.td.par <- as.numeric(fit.y$coefficients["z"]) * as.numeric(fit.z$coefficients["a"]) * (1 - 0)
                est.naive <-  mean(data$y[data$a==1]) - mean(data$y[data$a==0])

                #var.est.bd.par <- summary(fit.y.ac)$coefficients["a", "Std. Error"]^2*n
                c(i,  ate,  n, parameters[i, "alpha1"], parameters[i,"beta1"], parameters[i,"gamma1"], parameters[i,"gamma2"],
                  sapply(list(est.if.fd, est.if.bd, est.if.td, est.eif.td.bd, est.eif.fd.td, est.eif.all), mean),
                  est.fd.par, est.bd.par, est.td.par, est.naive,
                  sapply(list(est.if.fd, est.if.bd, est.if.td, est.eif.td.bd, est.eif.fd.td, est.eif.all), var),
                  bound.fd.value, bound.bd.value, bound.td.value, bound.eif.td.bd.value, bound.eif.fd.td.value, bound.eif.fd.td.bd.value)
        }
        sim.results <- rbind(sim.results, parOut)
        if (i%%100 ==0){
                save(sim.results, file = output.file.name)
        }
}


rm(a, h, i, n, S, parOut, bound.td.value, bound.fd.value, bound.bd.value, ate, sample.sizes, parameters, bound.eif.td.bd.value, cl)
rm(list = lsf.str())
save(list = ls(), file = output.file.name)

stopCluster(cl)



