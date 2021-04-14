## Functions for simulations

# data generation and additional functions --------------------------------
expit <- function(x){ exp(x) / (1 + exp(x)) }
p_c <- function(c_value, pc){ #
  # P(C = c_value)
  pc^c_value * (1 - pc)^(1 - c_value)
}
pa_given_c <-  function(a_value, c_value, par.a.c) {  
  # Calculates p(A = a_value|C = c_value)
  p1 <- expit(par.a.c %*% c(1, c_value)) 
  p1^a_value * (1 - p1)^(1 - a_value)
}
pz_given_a <- function(z_value, a_value, par.z.a, sigma.z) {  
  # Calculates p(Z = z_value | A = a_value)
  dnorm(z_value, par.z.a %*% c(1, a_value), sigma.z)
}
# py_given_zc <- function(y_value, z_value, c_value){  
#   # Calculates p(Y = y_value|Z = z_value, C = c_value)
#   dnorm(y_value, par.y.zc %*% c(1, z_value, c_value), sigma.y)
# }

gen_data_ca_binary_zy_cont <- function(n, pc, par.a.c, par.z.a, par.y.zc, sigma.z, sigma.y){
  c <- rbinom(n, 1, pc) # Generate confounder 
  a <- rbinom(n, 1, expit(cbind(rep(1,n),c)%*%par.a.c)) # Generate exposure (given c)
  z <- rnorm(n, cbind(rep(1,n), a )%*% par.z.a, sigma.z) # Generate continuous mediator (given a, c)
  y <- rnorm(n, cbind(rep(1,n),  z, c)%*% par.y.zc, sigma.y) # Generate continuous outcome Y (given a, m, c)
  sim.data <- data.frame(cbind(y, a, z, c))
  return(sim.data)
}

# Calculate true mean of potential outcome --------------------------------
psi_true_continuous <- function(astar, pc,  par.z.a, par.y.zc){
  # The true value of E(a*, Z(a*)) = int_{c,z} E(Y|Z=z,C=c)P(Z=z|A=a*, C=c)P(C=c)dcdz
  # according to the two-door criterion, which is true in our simulated data.
  # E(a*, Z(a*)) = int_{c,z} E(Y|Z=z,C=c)P(Z=z|A=a*, C=c)P(C=c)dcdz = |Z is indep of C given A|
  # = int_{c,z} E(Y|Z=z,C=c)P(Z=z|A=a*)P(C=c)dcdz = par.y.zc[1] + int_{c,z} par.y.zc[2]zP(Z=z|A=a*)P(C=c)dcdz
  # + int_{c,z} par.y.zc[3]cP(Z=z|A=a*)P(C=c)dcdz =
  # par.y.zc[1] + par.y.zc[2]*E(Z=z|A=a*) + par.y.zc[3]*E(C)
  mean.z.given.astar <- c(1, astar) %*% par.z.a
  psi <- c(1, mean.z.given.astar, pc)%*% par.y.zc
  return(psi)
}

# Calculate efficiency bounds ---------------------------------------------
BoundBackDoor <- function(astar, a, pc, par.a.c, par.y.zc, sigma.y, sigma.z){
  (sigma.y^2 + sigma.z^2 * par.y.zc[2]^2 )*
  (p_c(1, pc)/pa_given_c(astar, c_value = 1, par.a.c) + p_c(1, pc)/pa_given_c(a, c_value = 1, par.a.c)+
   p_c(0, pc)/pa_given_c(astar, c_value = 0, par.a.c) + p_c(0, pc)/pa_given_c(a, c_value = 0, par.a.c) )
}
BoundTwoDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z){
  sigma.y^2 * (integrate(function(z){1/sqrt(2*pi)/sigma.z *exp(-2*(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2)+(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2))}, -Inf, Inf)$value-1) * 
    (p_c(c_value=0, pc) * pa_given_c(a_value=a, c_value=0, par.a.c) + p_c(c_value=1, pc) * pa_given_c(a_value=a, c_value=1, par.a.c)) + 
  sigma.y^2 * (integrate(function(z){1/sqrt(2*pi)/sigma.z *exp(-2*(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2)+(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2))}, -Inf, Inf)$value-1)  * 
    (p_c(c_value=0, pc) * pa_given_c(a_value=astar, c_value=0, par.a.c) + p_c(c_value=1, pc) * pa_given_c(a_value=astar, c_value=1, par.a.c)) -
  sigma.y^2 * p_c(c_value=0, pc)*(1/pa_given_c(a_value=astar, c_value=0, par.a.c) + 1/pa_given_c(a_value=a, c_value=0, par.a.c)) -
  sigma.y^2 * p_c(c_value=1, pc)*(1/pa_given_c(a_value=astar, c_value=1, par.a.c) + 1/pa_given_c(a_value=a, c_value=1, par.a.c)) + 
      BoundBackDoor(astar=astar, a=a, pc=pc, par.a.c=par.a.c, par.y.zc=par.y.zc, sigma.y=sigma.y, sigma.z = sigma.z) 
}
BoundFrontDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z){
  fa <- function(a_value){
    pa_given_c(a_value = a_value, c_value = 0, par.a.c = par.a.c) * p_c(c_value=0, pc=pc) +
    pa_given_c(a_value = a_value, c_value = 1, par.a.c = par.a.c) * p_c(c_value=1, pc=pc) }
  
  EY_az <- function(a_value, z_value, par.a.c){
    Ec_given_a <- function(a_value){
     # (expit(c(1,1)%*% par.a.c)*a_value + (1 - expit(c(1,1)%*% par.a.c))*(1-a_value))/fa(a_value)*pc
      pa_given_c(a_value = a_value, c = 1, par.a.c = par.a.c)*pc/fa(a_value = a_value)
    }
    #c(1, z_value, Ec_given_a(a_value))%*%par.y.zc
    par.y.zc[1] + z_value * par.y.zc[2] + Ec_given_a(a_value) * par.y.zc[3]
  }
  
  # EY_az_coef <- function(){
  #   intercept <- par.y.zc[1]  + par.y.zc[3] *(1-expit(c(1,1)%*% par.a.c))/fa(a)*pc
  #   coef_a <- par.y.zc[3] * pc * (expit(c(1,1)%*% par.a.c)/fa(astar) - (1-expit(c(1,1)%*% par.a.c))/fa(a))
  #   return(c(intercept, par.y.zc[2], coef_a))
  # }
  
  fa_varY_az <- function(a_value, z_value){
    fa(a_value) * (sigma.y^2 + (par.y.zc[1] + par.y.zc[2] * z_value)^2) + 
      pa_given_c(a_value = a_value, c_value = 1, par.a.c) * p_c(c_value = 1, pc=pc) * 
      (par.y.zc[3]^2 + 2 * par.y.zc[3]*(par.y.zc[1] + par.y.zc[2]* z_value) ) - fa(a_value)*
      EY_az(a_value = a_value, z_value = z_value, par.a.c=par.a.c)^2
  }
  EZastar <- function(a_value){
    c(1, a_value) %*% par.z.a
  }
  integrate(Vectorize(function(z) {
            fa_varY_az(a_value = astar, z_value = z) * pz_given_a(z_value = z, astar, par.z.a, sigma.z) + # sum_z f(z|a^*)*fa_varY_az(a*) + f(z|a)*fa_varY_az(a)
            fa_varY_az(a_value = a, z_value = z) * pz_given_a(z_value = z, a, par.z.a, sigma.z)
  }), -Inf, Inf)$value -
  2 * integrate(Vectorize(function(z) {
                fa_varY_az(a_value = astar, z_value = z) * pz_given_a(z_value = z, a_value = a, par.z.a, sigma.z) + # sum_z  f(z|a)*fa_varY_az(a*) + f(z|a*)*fa_varY_az(a)
                fa_varY_az(a_value = a, z_value = z) * pz_given_a(z_value = z, a_value = astar, par.z.a, sigma.z)
      }), -Inf, Inf)$value +
  integrate(Vectorize(function(z) {
            fa_varY_az(a_value = astar, z_value = z) /sqrt(2*pi)/sigma.z * exp(-2 * (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) +
            fa_varY_az(a_value = a, z_value = z) /sqrt(2*pi)/sigma.z * exp(-2 * (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))
  }), -Inf, Inf)$value + # f^2(z|a)/f(z|a*)*fa_varY_az(a*) + f^2(z|a*)/f(z|a)*fa_varY_az(a)
  (par.y.zc[1] + par.y.zc[3] * pc)^2 * (1 / fa(a) + 1 / fa(astar)) + # sum_z f(z|a*)()^2/f(a*)) + same for a
  par.y.zc[2]^2 * (
    (sigma.z^2 + EZastar(astar)^2) / fa(astar) + (sigma.z^2 + EZastar(a)^2) / fa(a)
    ) +
  2 * par.y.zc[2] * (par.y.zc[1] + par.y.zc[3] * pc) * (EZastar(astar) / fa(astar) + EZastar(a) / fa(a)) -
  (c(1, EZastar(astar), pc) %*% par.y.zc)^2 / fa(astar) - #  - E^2Y(a*)/f(a*)
  (c(1, EZastar(a), pc) %*% par.y.zc)^2 / fa(a) + #  - E^2Y(a)/f(a)
  par.y.zc[2]^2 * (EZastar(astar) - EZastar(a))^2 - # + sum_a f(a)* (sum_z EY|a,z (f(z|a*)- f(z|a)))^2
  (psi_true_continuous(astar, pc, par.z.a, par.y.zc) - psi_true_continuous(a, pc, par.z.a, par.y.zc))^2 #  - ATE^2
}

BoundEIFTwoBackDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z){
        sigma.y^2 * p_c(c_value=0, pc) * integrate(Vectorize(function(z){
                1/sqrt(2*pi*sigma.z^2) * ( 
                        exp(-(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2)) +
                                exp(-2*(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2) + (z-as.numeric(c(1,1)%*%par.z.a))^2/(2*sigma.z^2)) - 
                                2 * exp(-(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2)))/
                        ( pa_given_c(a_value=astar, c_value=0, par.a.c) +
                                  pa_given_c(a_value=a, c_value=0, par.a.c) * exp(-(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2)+(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2))
                        )
        }), -Inf, Inf)$value +
                sigma.y^2 * p_c(c_value=1, pc) * integrate(Vectorize(function(z){
                        1/sqrt(2*pi*sigma.z^2) * ( 
                                exp(-(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2)) +
                                        exp(-2*(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2) + (z-as.numeric(c(1,1)%*%par.z.a))^2/(2*sigma.z^2)) - 
                                        2 * exp(-(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2)))/
                                ( pa_given_c(a_value=astar, c_value=1, par.a.c) +
                                  pa_given_c(a_value=a, c_value=1, par.a.c) * exp(-(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2)+(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2))
                                )
                }), -Inf, Inf)$value  -
                sigma.y^2 * p_c(c_value=0, pc)*(1/pa_given_c(a_value=astar, c_value=0, par.a.c) + 1/pa_given_c(a_value=a, c_value=0, par.a.c)) -
                sigma.y^2 * p_c(c_value=1, pc)*(1/pa_given_c(a_value=astar, c_value=1, par.a.c) + 1/pa_given_c(a_value=a, c_value=1, par.a.c)) + 
                BoundBackDoor(astar=astar, a=a, pc=pc, par.a.c=par.a.c, par.y.zc=par.y.zc, sigma.y=sigma.y, sigma.z = sigma.z) 
}

# TD ---------------------------------------------------------------------
EstimateInfluenceFunctionTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar){
  n <- length(exposure)
  sigma.z <- summary(fit.z)$sigma
  
  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,"a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  y.par.hat <- summary(fit.y)$coefficients[,1]
  z.par.hat <- summary(fit.z)$coefficients[,1]
  a.par.hat <- summary(fit.a)$coefficients[,1]

  z.mean_astar <- model.matrix.z_astar%*%z.par.hat
  z.mean_ind <- model.matrix.z%*%z.par.hat
  a.mean <- expit(model.matrix.a%*%a.par.hat)
  y.mean <- model.matrix.y%*%y.par.hat
  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- y.mean
  sum.z <- as.matrix(cbind(rep(1,n),z.mean_astar,cov.vals.all))%*%y.par.hat
  sum.az <- as.matrix(cbind(rep(1,n),z.mean_ind,cov.vals.all))%*%y.par.hat
  psi.td.ind <-  (outcome - y.mean)*
                    (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/dnorm(model.matrix.y[,"z"],z.mean_ind,sigma.z))+
                   as.numeric(exposure==astar)/(a.mean*astar + (1 - a.mean)*(1 - astar))*(sum.a - sum.az) +
                  sum.z
  # psi.td.hat <- mean(psi.td.ind)
  
  return(psi.td.ind)
}
# FD ----------------------------------------------------------------------
EstimateInfluenceFunctionFrontDoor <- function(exposure, intermediate, outcome, fit.z, fit.y.az, astar){
  
  n <- length(exposure)
  
  sigma.z <- summary(fit.z)$sigma
  
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.az)))
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,"a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  y.par.hat <- summary(fit.y.az)$coefficients[,1]
  z.par.hat <- summary(fit.z)$coefficients[,1]

  z.mean_astar <- model.matrix.z_astar%*%z.par.hat
  z.mean_ind <- model.matrix.z%*%z.par.hat
  a.mean <- mean(exposure) # E(a) = P(A=1)
  y.mean <- model.matrix.y%*%y.par.hat
  
  
  sum.a <- as.matrix(cbind(rep(1,n),rep(a.mean,n), model.matrix.y[,"z"]))%*%y.par.hat
  sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,"a"],z.mean_astar))%*%y.par.hat
  sum.az <- as.matrix(cbind(rep(1,n),rep(a.mean,n),z.mean_astar))%*%y.par.hat
  
  psi.fd.ind <-  (outcome - y.mean)*
                    (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/dnorm(model.matrix.y[,"z"],z.mean_ind,sigma.z))+
                     as.numeric(exposure==astar)/(a.mean*astar + (1 - a.mean)*(1 - astar))*(sum.a - sum.az)+
                     sum.z
  
  # psi.fd.hat <- mean(psi.fd.ind)
  return(psi.fd.ind)
  
}
# BD ----------------------------------------------------------------------
EstimateInfluenceFunctionBackDoor <- function(cov.vals.all, exposure, outcome, fit.a, fit.z, fit.y.ac, astar){
  
  n <- length(exposure)
  
  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.ac)))
  
  model.matrix.y.astar <- data.frame(model.matrix.y)
  model.matrix.y.astar[, "a"] <- astar
  model.matrix.y.astar <- as.matrix(model.matrix.y.astar)
  
  y.par.hat <- summary(fit.y.ac)$coefficients[,1]
  a.par.hat <- summary(fit.a)$coefficients[,1]

  a.mean <- expit(model.matrix.a%*%a.par.hat)
  y.mean <- model.matrix.y%*%y.par.hat
  
  sum.c <- model.matrix.y.astar%*%y.par.hat
  
  
  psi.bd.ind <-  as.numeric(exposure==astar)/(a.mean*astar + (1 - a.mean)*(1 - astar))*(outcome - sum.c) +
                  sum.c 
  # psi.bd.hat <- mean(psi.bd.ind)
  return(psi.bd.ind)
}
# EIF TD BD ----------------------------------------------------------------------
EstimateEIFTwoBackDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar){
        n <- length(exposure)
        sigma.z <- summary(fit.z)$sigma
        
        model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
        model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
        model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))
        
        model.matrix.z0 <- data.frame(model.matrix.z)
        model.matrix.z0[,"a"] <- 0
        model.matrix.z0 <- as.matrix(model.matrix.z0)
        
        model.matrix.z1 <- data.frame(model.matrix.z)
        model.matrix.z1[,"a"] <- 1
        model.matrix.z1 <- as.matrix(model.matrix.z1)
        
        model.matrix.z_astar <- data.frame(model.matrix.z)
        model.matrix.z_astar[,"a"] <- astar
        model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
        
        y.par.hat <- summary(fit.y)$coefficients[,1]
        z.par.hat <- summary(fit.z)$coefficients[,1]
        a.par.hat <- summary(fit.a)$coefficients[,1]
        
        z.mean_astar <- model.matrix.z_astar%*%z.par.hat
        z.mean_ind <- model.matrix.z%*%z.par.hat
        a.mean <- expit(model.matrix.a%*%a.par.hat)
        y.mean <- model.matrix.y%*%y.par.hat
        # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
        sum.a <- y.mean
        sum.z <- as.matrix(cbind(rep(1,n),z.mean_astar,cov.vals.all))%*%y.par.hat
        sum.az <- as.matrix(cbind(rep(1,n),z.mean_ind,cov.vals.all))%*%y.par.hat
        sum.denom.az <- expit(cbind(1,cov.vals.all) %*%a.par.hat) * dnorm(model.matrix.y[,"z"],model.matrix.z1%*%z.par.hat,sigma.z) + 
                   (1 - expit(cbind(1,cov.vals.all)%*%a.par.hat)) * dnorm(model.matrix.y[,"z"],model.matrix.z0%*%z.par.hat,sigma.z)
       
        psi.td.ind <-  (outcome - y.mean)*
                (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/sum.denom.az) +
                as.numeric(exposure==astar)/(a.mean*astar + (1 - a.mean)*(1 - astar))*(sum.a - sum.az) +
                sum.z
        # psi.td.hat <- mean(psi.td.ind)
        
        return(psi.td.ind)
}

# InfluenceFunctionFrontDoor <- function(exposure, intermediate, outcome, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z, astar){
#   
#   n <- length(exposure)
#   
#   model.matrix.z <- as.matrix(data.frame(Intercept = rep(1,n), a=exposure ))
#   model.matrix.y <- as.matrix(data.frame(Intercept = rep(1,n), a=exposure, z=intermediate ))
#   
#   model.matrix.z_astar <- model.matrix.z
#   model.matrix.z_astar[,"a"] <- astar
#   model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
#   
#   y.par.hat <- c(par.y.zc[1] + par.y.zc[3]*0.5*(1-expit(0.5))/(1/4 +1/2*(1-expit(0.5))),
#                  par.y.zc[3]*0.5*expit(0.5)/(1/4 +1/2*expit(0.5)) - par.y.zc[3]* 0.5*(1-expit(0.5))/(1/4 + 1/2*(1-expit(0.5))),
#                  par.y.zc[2] )
#   z.par.hat <- par.z.a
#   
#   z.mean_astar <- model.matrix.z_astar%*%z.par.hat
#   z.mean_ind <- model.matrix.z%*%z.par.hat
#   a.mean <- 1/4 +1/2*expit(0.5) # E(a) = P(A=1)
#   y.mean <- model.matrix.y%*%y.par.hat
#   
#   
#   sum.a <- as.matrix(cbind(rep(1,n),rep(a.mean,n), model.matrix.y[,"z"]))%*%y.par.hat
#   sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,"a"],z.mean_astar))%*%y.par.hat
#   sum.az <- as.matrix(cbind(rep(1,n),rep(a.mean,n),z.mean_astar))%*%y.par.hat
#   
#   psi.fd.ind <-  (outcome - y.mean)*
#     (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/dnorm(model.matrix.y[,"z"],z.mean_ind,sigma.z))+
#     as.numeric(exposure==astar)/(a.mean*astar + (1 - a.mean)*(1 - astar))*(sum.a - sum.az)+
#     sum.z
#   
#   #psi.fd.hat <- mean(psi.fd.ind)
#   psi.fd.hat <-  psi.fd.ind
#   output <- psi.fd.hat
#   return(output)
#   
# }
# 
# if.fd.true.var <- function(exposure, intermediate, outcome, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z, astar){
#   
#   n <- length(exposure)
#   
#   model.matrix.z <- as.matrix(data.frame(Intercept = rep(1,n), a=exposure ))
#   model.matrix.y <- as.matrix(data.frame(Intercept = rep(1,n), a=exposure, z=intermediate ))
#   
#   model.matrix.z_astar <- model.matrix.z
#   model.matrix.z_astar[,"a"] <- astar
#   model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
#   
#   y.par.hat <- c(par.y.zc[1] + par.y.zc[3]*0.5*(1-expit(0.5))/(1/4 +1/2*(1-expit(0.5))),
#                  par.y.zc[3]*0.5*expit(0.5)/(1/4 +1/2*expit(0.5)) - par.y.zc[3]* 0.5*(1-expit(0.5))/(1/4 + 1/2*(1-expit(0.5))),
#                  par.y.zc[2] )
#   z.par.hat <- par.z.a
#   
#   z.mean_astar <- model.matrix.z_astar%*%z.par.hat
#   z.mean_ind <- model.matrix.z%*%z.par.hat
#   a.mean <- 1/4 +1/2*expit(0.5) # E(a) = P(A=1)
#   y.mean <- model.matrix.y%*%y.par.hat
#   
#   
#   sum.a <- as.matrix(cbind(rep(1,n),rep(a.mean,n), model.matrix.y[,"z"]))%*%y.par.hat
#   sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,"a"],z.mean_astar))%*%y.par.hat
#   sum.az <- as.matrix(cbind(rep(1,n),rep(a.mean,n),z.mean_astar))%*%y.par.hat
#   
#    psi.fd.ind <-  (outcome - y.mean)*#*
#     (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/dnorm(model.matrix.y[,"z"],z.mean_ind,sigma.z)) #+
#     #  as.numeric(exposure==astar)/(a.mean*astar + (1 - a.mean)*(1 - astar))*(sum.a - sum.az)+
#     # sum.z
#   
#   #psi.fd.hat <- mean(psi.fd.ind)
#   
#   output <- psi.fd.ind
#   return(output)
#   
# }
# 


# EIF TD FD ---------------------------------------------------------------
EstimateEIFTwoFrontDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar){
     
        n <- length(exposure)
        sigma.z <- summary(fit.z)$sigma
        
        model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
        model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
        model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))
        
        model.matrix.z_astar <- data.frame(model.matrix.z)
        model.matrix.z_astar[,"a"] <- astar
        model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
        
        y.par.hat <- summary(fit.y)$coefficients[,1]
        z.par.hat <- summary(fit.z)$coefficients[,1]
        a.par.hat <- summary(fit.a)$coefficients[,1]
        
        z.mean_astar <- model.matrix.z_astar%*%z.par.hat
        z.mean_ind <- model.matrix.z%*%z.par.hat
        
        a.mean <- expit(model.matrix.a%*%a.par.hat)
        
        pa <-  t( (expit(c(1,1)%*%a.par.hat)%*%exposure + (1-expit(c(1,1)%*%a.par.hat))%*%(1-exposure))* mean(cov.vals.all) +  
                (expit(c(1,0)%*%a.par.hat)%*%exposure + (1-expit(c(1,0)%*%a.par.hat))%*%(1-exposure))* (1-mean(cov.vals.all)))
                
        y.mean <- model.matrix.y%*%y.par.hat
        
        
        # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
        sum.a <- model.matrix.y[,c("X.Intercept.", "z" )]%*%y.par.hat[c(1,2)] + y.par.hat[3]*mean(cov.vals.all)
        sum.z <- as.matrix(cbind(rep(1,n),z.mean_astar,cov.vals.all))%*%y.par.hat
        sum.az <- as.matrix(cbind(rep(1,n),z.mean_ind,mean(cov.vals.all)))%*%y.par.hat
        psi.td.ind <-  (outcome - y.mean)*
                (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/dnorm(model.matrix.y[,"z"],z.mean_ind,sigma.z))+
                as.numeric(exposure==astar)/pa*(sum.a - sum.az) +
                sum.z
        # psi.td.hat <- mean(psi.td.ind)
        
        return(psi.td.ind)
}


# EIF All -----------------------------------------------------------------
EstimateEIFAll <-  function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar){
        n <- length(exposure)
        sigma.z <- summary(fit.z)$sigma
        
        model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
        model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
        model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))
        
        model.matrix.z0 <- data.frame(model.matrix.z)
        model.matrix.z0[,"a"] <- 0
        model.matrix.z0 <- as.matrix(model.matrix.z0)
        
        model.matrix.z1 <- data.frame(model.matrix.z)
        model.matrix.z1[,"a"] <- 1
        model.matrix.z1 <- as.matrix(model.matrix.z1)
        
        model.matrix.z_astar <- data.frame(model.matrix.z)
        model.matrix.z_astar[,"a"] <- astar
        model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
        
        y.par.hat <- summary(fit.y)$coefficients[,1]
        z.par.hat <- summary(fit.z)$coefficients[,1]
        a.par.hat <- summary(fit.a)$coefficients[,1]
        
        z.mean_astar <- model.matrix.z_astar%*%z.par.hat
        z.mean_ind <- model.matrix.z%*%z.par.hat
        a.mean <- expit(model.matrix.a%*%a.par.hat)
        y.mean <- model.matrix.y%*%y.par.hat
        
        # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
        sum.a <- model.matrix.y[,c("X.Intercept.", "z" )]%*%y.par.hat[c(1,2)] + y.par.hat[3]*mean(cov.vals.all)
        sum.z <- as.matrix(cbind(rep(1,n),z.mean_astar,cov.vals.all))%*%y.par.hat
        sum.az <- as.matrix(cbind(rep(1,n),z.mean_ind, mean(cov.vals.all)))%*%y.par.hat
        sum.denom.az <- expit(cbind(1,cov.vals.all) %*%a.par.hat) * dnorm(model.matrix.y[,"z"],model.matrix.z1%*%z.par.hat,sigma.z) + 
                (1 - expit(cbind(1,cov.vals.all)%*%a.par.hat)) * dnorm(model.matrix.y[,"z"],model.matrix.z0%*%z.par.hat,sigma.z)
        
        pa <-  t( (expit(c(1,1)%*%a.par.hat)%*%exposure + (1-expit(c(1,1)%*%a.par.hat))%*%(1-exposure))* mean(cov.vals.all) +  
                          (expit(c(1,0)%*%a.par.hat)%*%exposure + (1-expit(c(1,0)%*%a.par.hat))%*%(1-exposure))* (1-mean(cov.vals.all)))
        
        psi.td.ind <-  (outcome - y.mean)*
                (dnorm(model.matrix.y[,"z"],z.mean_astar,sigma.z)/sum.denom.az) +
                as.numeric(exposure==astar)/pa*(sum.a - sum.az) +
                sum.z
        # psi.td.hat <- mean(psi.td.ind)
        
        return(psi.td.ind)
}
