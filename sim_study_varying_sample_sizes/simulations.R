# code for the simulation study
require(foreach)
require(doParallel)
require(doRNG) # for independent RNG streams.

cl <- makeCluster(10) 
registerDoParallel(cl)
source("../sim_study_varying_sample_sizes/if_bounds_utils.R") # load functions
source("../sim_study_varying_sample_sizes/sim_study_parameters.R") # load parameters
set.seed(seed)


sim.results <- matrix(nrow = 0, ncol=22)

h <- 0
for (i in 1:length(parameters)){
        h <- h+1
        print(h)
        ate <- psi_true_continuous(astar = 1, pc = parameters[[i]]$pc, par.z.a = parameters[[i]]$par.z.a,
                                   par.y.zc = parameters[[i]]$par.y.zc)-
                psi_true_continuous(astar = 0, pc = parameters[[i]]$pc, par.z.a = parameters[[i]]$par.z.a,
                                    par.y.zc = parameters[[i]]$par.y.zc)
        #### Bounds ###
        bound.bd.value <- BoundBackDoor(astar = 1, a = 0, pc=parameters[[i]]$pc,
                                        par.a.c = parameters[[i]]$par.a.c,
                                        par.y.zc = parameters[[i]]$par.y.zc,
                                        sigma.y = parameters[[i]]$sigma.y,
                                        sigma.z = parameters[[i]]$sigma.z)
        
        bound.td.value <- BoundTwoDoor(astar = 1, a = 0,  pc = parameters[[i]]$pc,
                                       par.a.c  = parameters[[i]]$par.a.c,
                                       par.z.a  = parameters[[i]]$par.z.a,
                                       par.y.zc = parameters[[i]]$par.y.zc,
                                       sigma.y  = parameters[[i]]$sigma.y,
                                       sigma.z  = parameters[[i]]$sigma.z)
        
        bound.fd.value <- BoundFrontDoor(astar = 1, a = 0, pc=parameters[[i]]$pc,
                                         par.a.c=parameters[[i]]$par.a.c,
                                         par.z.a = parameters[[i]]$par.z.a,
                                         par.y.zc = parameters[[i]]$par.y.zc,
                                         sigma.y = parameters[[i]]$sigma.y,
                                         sigma.z = parameters[[i]]$sigma.z)
        
        
        bound.eif.td.bd.value <- BoundEIFTwoBackDoor(astar = 1, a = 0, pc=parameters[[i]]$pc,
                                         par.a.c=parameters[[i]]$par.a.c,
                                         par.z.a = parameters[[i]]$par.z.a,
                                         par.y.zc = parameters[[i]]$par.y.zc,
                                         sigma.y = parameters[[i]]$sigma.y,
                                         sigma.z = parameters[[i]]$sigma.z)
       for (n in sample.size){
               parOut <- foreach(s=1:number.of.replicates, .combine='rbind') %dorng% {
                       # Generate data from DGP i with sample size n
                       data <- gen_data_ca_binary_zy_cont(n,
                                                          pc = parameters[[i]]$pc,
                                                          par.a.c = parameters[[i]]$par.a.c,
                                                          par.z.a = parameters[[i]]$par.z.a,
                                                          par.y.zc = parameters[[i]]$par.y.zc,
                                                          sigma.z = parameters[[i]]$sigma.z,
                                                          sigma.y = parameters[[i]]$sigma.y)

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
                               EstimateInfluenceFunctionTwoDoor(cov.vals.all = data$c, exposure = data$a,intermediate = data$m,outcome = data$y,
                                                                fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)

                       est.eif.td.bd <- EstimateEIFTwoBackDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                                     fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                               EstimateEIFTwoBackDoor(cov.vals.all = data$c, exposure = data$a,intermediate = data$m,outcome = data$y,
                                                                fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)
                       
                       est.eif.td.fd <- EstimateEIFTwoFrontDoor(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                                     fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                               EstimateEIFTwoFrontDoor(cov.vals.all = data$c, exposure = data$a,intermediate = data$m,outcome = data$y,
                                                                fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)
                       
                       est.eif.all <- EstimateEIFAll(cov.vals.all= data$c, exposure = data$a, intermediate = data$z, outcome = data$y,
                                                               fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 1) -
                               EstimateEIFAll(cov.vals.all = data$c, exposure = data$a,intermediate = data$m,outcome = data$y,
                                                       fit.a = fit.a, fit.z = fit.z, fit.y = fit.y, astar = 0)
                       

                       # Parametric estimates
                       est.bd.par <- as.numeric(fit.y.ac$coefficients["a"])
                       est.fd.par <- as.numeric(fit.y.az$coefficients["z"]) * as.numeric(fit.z$coefficients["a"]) * (1 - 0)
                       est.td.par <- as.numeric(fit.y$coefficients["z"]) * as.numeric(fit.z$coefficients["a"]) * (1 - 0)

                       #var.est.bd.par <- summary(fit.y.ac)$coefficients["a", "Std. Error"]^2*n
                       c(i,  ate,  n, 
                         sapply(list(est.if.fd, est.if.bd, est.if.td, est.eif.td.bd, est.eif.td.fd, est.eif.all), mean),
                         est.fd.par, est.bd.par, est.td.par,
                         sapply(list(est.if.fd, est.if.bd, est.if.td, est.eif.td.bd, est.eif.td.fd, est.eif.all), var),
                         bound.fd.value, bound.bd.value, bound.td.value, bound.eif.td.bd.value)
               }
              sim.results <- rbind(sim.results, parOut)
       }
}
colnames(sim.results) <-  c("dgp", "ate", "n", 
                             "est.fd.semipar", "est.bd.semipar", "est.td.semipar", "est.eif.td.bd.semipar", "est.eif.td.fd.semipar","est.eif.all.semipar",
                              "est.fd.par", "est.bd.par", "est.td.par", 
                              "est.var.fd.semipar", "est.var.bd.semipar", "est.var.td.semipar", "est.var.td.bd.semipar","est.var.td.fd.semipar","est.var.all.semipar",
                              "bound.fd", "bound.bd", "bound.td", "bound.td.bd")

rm(a, h, i, n, S, parOut, bound.td.value, bound.fd.value, bound.bd.value, ate, sample.sizes, parameters)
rm(list = lsf.str())
save(list = ls(), file = paste0("results_sim_study_",  Sys.Date(), ".Rdata"))




