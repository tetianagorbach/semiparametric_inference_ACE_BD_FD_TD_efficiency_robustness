# code for the simulation study
require(foreach)
require(doParallel)
require(doRNG) # for independent RNG streams.
require(tidyverse)

source("sim_study_misspecification/Rscripts/if_bounds_utils.R") # load functions
source("sim_study_misspecification/Rscripts/sim_study_misspecification_parameters.R") # load parameters
seed <- 1668613

cl <- makeCluster(number.of.clusters) 
registerDoParallel(cl)
set.seed(seed)


sim.results <- matrix(nrow = 0, ncol=8)


estimates <- function(pc.miss, 
                      fit.a.c.miss, pa.miss,
                      fit.z.a.miss, 
                      fit.y.zc.miss,
                      fit.y.ac.miss,
                      fit.y.az.miss ) {
        est.if.bd <- EstimateInfluenceFunctionBackDoor(cov.vals.all = data$c, exposure = data$a, outcome = data$y,
                                                       fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.ac = fit.y.ac.miss, astar = 1) -
                EstimateInfluenceFunctionBackDoor(cov.vals.all = data$c, exposure = data$a, outcome = data$y,
                                                  fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.ac = fit.y.ac.miss, astar = 0)
        
        est.if.fd <- EstimateInfluenceFunctionFrontDoor(exposure = data$a, intermediate = data$z, outcome = data$y, pa = pa.miss,
                                                        fit.z.a = fit.z.a.miss, fit.y.az = fit.y.az.miss, astar = 1) -
                EstimateInfluenceFunctionFrontDoor(exposure = data$a, intermediate = data$z, outcome = data$y, pa = pa.miss,
                                                   fit.z.a = fit.z.a.miss, fit.y.az = fit.y.az.miss,  astar = 0)
        
        
        
        est.if.td <- EstimateInfluenceFunctionTwoDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                      fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 1) -
                EstimateInfluenceFunctionTwoDoor(cov.vals.all = data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                 fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 0)
        
        est.eif.bd.td <- EstimateEIFTwoBackDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 1) -
                EstimateEIFTwoBackDoor(cov.vals.all = data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                       fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 0)
        
        est.eif.fd.td <- EstimateEIFFrontTwoDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                 est.pc = pc.miss, fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 1) -
                EstimateEIFFrontTwoDoor(cov.vals.all = data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                        est.pc = pc.miss, fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 0)
        
        est.eif.all <- EstimateEIFAll(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                      est.pc = pc.miss,  fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 1) -
                EstimateEIFAll(cov.vals.all = data$c, exposure = data$a,intermediate = data$z,outcome = data$y,
                               est.pc = pc.miss, fit.a.c = fit.a.c.miss, fit.z.a = fit.z.a.miss, fit.y.zc = fit.y.zc.miss, astar = 0)
        
        a <- apply(cbind(
                est.if.bd, est.if.fd, est.if.td,
                est.eif.bd.td, est.eif.fd.td, est.eif.all
        ),
        MARGIN = 2, FUN = mean) # calculate estimates
        names(a) <-  c("BD", "FD", "TD", "BD TD", "FD TD", "BD FD TD")
        a
        
}
for (i in 1:nrow(parameters)){
        print(i)
        pc <- pc
        par.a.c <-  as.numeric(c(0, parameters$alpha[i]))
        par.z.a <-  as.numeric(c(0, parameters$beta[i]))
        par.y.zc <-  as.numeric(c(0, parameters$gamma1[i], parameters$gamma2[i]))
        
        ate <- psi_true_continuous(astar = 1, pc = pc, par.z.a = par.z.a, par.y.zc = par.y.zc)-
                psi_true_continuous(astar = 0, pc = pc, par.z.a = par.z.a, par.y.zc = par.y.zc)
        #### Bounds ###
        bound.bd <- BoundBackDoor(astar = 1, a = 0, pc = pc,
                                  par.a.c = par.a.c,
                                  par.y.zc = par.y.zc,
                                  sigma.y = sigma.y,
                                  sigma.z = sigma.z)
        
        bound.fd <- BoundFrontDoor(astar = 1, a = 0, pc = pc,
                                   par.a.c = par.a.c,
                                   par.z.a = par.z.a,
                                   par.y.zc = par.y.zc,
                                   sigma.y = sigma.y,
                                   sigma.z = sigma.z)
        
        bound.td <- BoundTwoDoor(astar = 1, a = 0, pc = pc,
                                 par.a.c = par.a.c,
                                 par.z.a  = par.z.a,
                                 par.y.zc = par.y.zc,
                                 sigma.y = sigma.y,
                                 sigma.z = sigma.z)
        
        
        
        bound.bd.td <- BoundEIFTwoBackDoor(astar = 1, a = 0, pc = pc,
                                           par.a.c = par.a.c,
                                           par.z.a = par.z.a,
                                           par.y.zc = par.y.zc,
                                           sigma.y = sigma.y,
                                           sigma.z = sigma.z)
        
        bound.fd.td <- BoundEIFFrontTwoDoor(astar = 1, a = 0, pc = pc,
                                            par.a.c = par.a.c,
                                            par.z.a = par.z.a,
                                            par.y.zc = par.y.zc,
                                            sigma.y = sigma.y,
                                            sigma.z = sigma.z)
        
        bound.bd.fd.td <- BoundEIFFrontTwoBackDoor(astar = 1, a = 0, pc = pc,
                                                   par.a.c = par.a.c,
                                                   par.z.a = par.z.a,
                                                   par.y.zc = par.y.zc,
                                                   sigma.y = sigma.y,
                                                   sigma.z = sigma.z)
        
        bounds <-  c("bound.BD" = bound.bd, "bound.FD" = bound.fd, "bound.TD" = bound.td,
                     "bound.BDTD" = bound.bd.td, "bound.FDTD" = bound.fd.td, "bound.BDFDTD" = bound.bd.fd.td)
        
        
        
        parOut <- foreach(s=1:number.of.replicates, .combine='rbind') %dorng% {
                # Generate data from DGP i with sample size n
                data <- gen_data_ca_binary_zy_cont(n = sample.size,
                                                   pc = pc,
                                                   par.a.c = par.a.c,
                                                   par.z.a = par.z.a,
                                                   par.y.zc = par.y.zc,
                                                   sigma.z = sigma.z,
                                                   sigma.y = sigma.y)
                
                # Fit nuisance models
                pc <-  mean(data$c)
                fit.a.c <- glm(a ~ c, data = data, family = binomial)
                pa <- mean(data$a)
                fit.z.a <- lm(z ~ a , data = data)
                fit.y.zc <- lm(y ~  z  + c, data = data)
                fit.y.az <- lm(y ~ a + z, data = data)
                fit.y.ac <- lm(y ~ a + c, data = data)
                
                fit.a.c.misspecified <- fit.a.c
                fit.a.c.misspecified$coefficients[ "(Intercept)"] <-  coef(glm(a ~ 1, data = data, family = binomial))
                fit.a.c.misspecified$coefficients[ "c"] <- 0
                
                fit.z.a.misspecified <- fit.z.a
                fit.z.a.misspecified$coefficients[ "(Intercept)"]  <-  lm(z ~ 1 , data = data)$coefficients[ "(Intercept)"] 
                fit.z.a.misspecified$coefficients[ "a"] <- 0
                
                
                fit.y.zc.misspecified <- fit.y.zc 
                fit.y.zc.misspecified$coefficients["z"] <- 0
                fit.y.zc.misspecified$coefficients[c("(Intercept)", "c")] <- lm(y ~ c, data = data)$coefficients[c("(Intercept)", "c")]
                
                
                
                fit.y.az.misspecified <-  fit.y.az
                fit.y.az.misspecified$coefficients[c("(Intercept)", "a")]  <-  lm(y ~ a, data = data)$coefficients[c("(Intercept)", "a")]
                fit.y.az.misspecified$coefficients["z"] <- 0
                
                fit.y.ac.misspecified <-  fit.y.ac
                fit.y.ac.misspecified$coefficients[c("(Intercept)", "a")] <-  lm(y ~ a, data = data)$coefficients[c("(Intercept)", "a")]
                fit.y.ac.misspecified$coefficients["c"] <- 0
                
                
                
                # Semiparametric estimates
                
                
                # Misspecified p(Z|A) - all consistent
                est1 <-  c("misspecification" = 1,   "ate" = ate,
                           estimates(pc.miss = pc, 
                                     fit.a.c.miss = fit.a.c, pa.miss = pa,
                                     fit.z.a.miss = fit.z.a.misspecified, 
                                     fit.y.zc.miss = fit.y.zc,
                                     fit.y.ac.miss = fit.y.ac,
                                     fit.y.az.miss = fit.y.az))
                
                # Only p(Z|A) correctly specified - TD unbiased, but BD&TD is biased
                est2 <- c("misspecification" = 2,   "ate" = ate,
                          estimates(pc.miss = pc/4, 
                                    fit.a.c.miss = fit.a.c.misspecified, pa.miss = pa/4,
                                    fit.z.a.miss = fit.z.a, 
                                    fit.y.zc.miss = fit.y.zc.misspecified,
                                    fit.y.ac.miss = fit.y.ac.misspecified,
                                    fit.y.az.miss = fit.y.az.misspecified))
                # Misspecified p(Z|A) & p(C) - FD, TD, FD&TD - unbiased, but FD&TD, BD&FD&TD  - biased
                est3 <- c("misspecification" = 3,  "ate" = ate,
                          estimates(pc.miss = pc/4, 
                                    fit.a.c.miss = fit.a.c, pa.miss = pa,
                                    fit.z.a.miss = fit.z.a.misspecified, 
                                    fit.y.zc.miss = fit.y.zc,
                                    fit.y.ac.miss = fit.y.ac,
                                    fit.y.az.miss = fit.y.az
                          ))
                # Misspecified p(A|C) & E(Y|Z,C) - BD, TD - unbiased; BD&TD, BD&FD&TD - biased. 
                est4 <- c("misspecification" = 4,  "ate" = ate,
                          estimates(pc.miss = pc, 
                                    fit.a.c.miss = fit.a.c.misspecified, pa.miss = pa,
                                    fit.z.a.miss = fit.z.a, 
                                    fit.y.zc.miss = fit.y.zc.misspecified,
                                    fit.y.ac.miss = fit.y.ac,
                                    fit.y.az.miss = fit.y.az))
                rbind(est1, est2, est3, est4)
        }
        sim.results <- as.data.frame(parOut)%>%
                mutate(misspecification = recode(misspecification,  "1" ="Z", "2" = "C, A, Y", "3" = "C, Z", "4" = "A, Y" ))
        output <-  list(sim.results = sim.results, parameters =  c(unlist(parameters), bounds))
        save (output, file = output.file.name)
}

rm(i,  parOut,  parameters,  cl, bound.td, bound.fd, bound.bd, ate,  bound.bd.td, bound.fd.td, bound.bd.fd.td, bounds, sim.results)
rm(list = lsf.str())
save(list = ls(), file = output.file.name)
