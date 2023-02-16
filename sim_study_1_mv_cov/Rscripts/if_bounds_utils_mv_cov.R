## Functions for simulations

# Additional functions --------------------------------
expit <- function(x) {
  exp(x) / (1 + exp(x))
}
#' p_c <- function(c_value, pc) { #
#'   #' Calculates the probability of a Bernoulli(pc) random variable to be equal to c_value
#'   #' @param c_value Value of a random variable
#'   #' @param pc Parameter of the distribution
#'   #'
#'   #' @return Probability of a Bern(pc) random variable to be equal to c_value
#'   pc^c_value * (1 - pc)^(1 - c_value)
#' }
#' pa_given_c <- function(a_value, c_value, par.a.c) {
#'   #' Calculates p(A = a_value|C = c_value)
#'   #' @param a_value Value of A
#'   #' @param c_value Value of C
#'   #' @param par.a.c Parameters of the distribution of A
#'   #'
#'   #' @return Probability p(A = a_value|C = c_value)
#'   p1 <- expit(par.a.c %*% c(1, c_value))
#'   p1^a_value * (1 - p1)^(1 - a_value)
#' }
pz_given_a <- function(z_value, a_value, par.z.a, sigma.z) {
  #' Calculates the density p(Z = z_value | A = a_value)
  #' @param z_value Value of Z
  #' @param a_value Value of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param sigma.z Standard deviation of Z(a_value)
  #'
  #' @return Density p(Z = z_value | A = a_value)
  dnorm(z_value, par.z.a %*% c(1, a_value), sigma.z)
}


# Generate data -----------------------------------------------------------
GenerateData <- function(n, pc, par.a.c, par.z.a, par.y.zc, sigma.z, sigma.y) {
  #' Generates the data for the simulation study
  #' @param n Number of observation to generate
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return  A data frame with generated values of the observed variables C, A, Z, Y
        c1 <- rbinom(n, 1, pc) # Generate a binary confounder
        c2 <- rnorm(n, 0, diag(1)) # Generate a continuous confounder
        c3 <- MASS::mvrnorm(n, c(0,0,0), diag(3)) # Generate a multivariate confounder
        a <- rbinom(n, 1, expit(cbind(rep(1, n), c1, c2, c3) %*% par.a.c)) # Generate the exposure (given c)
        z <- rnorm(n, cbind(rep(1, n), a) %*% par.z.a, sigma.z) # Generate  the mediator (given a and c)
        y <- rnorm(n, cbind(rep(1, n), z, c1, c2, c3) %*% par.y.zc, sigma.y) # Generate the outcome Y (given a, z, c)
        sim.data <- data.frame(cbind(y, a, z, c1 , c2, c3))
        names(sim.data) <-  c("y", "a", "z", "c1", "c2", "c31", "c32", "c33")
  return(sim.data)
}





# Estimate influence functions for EY(astar)------------------------------------
EstimateIFBackDoor <- function(cov.vals.all, exposure, outcome, fit.a.c, fit.z.a, fit.y.ac, astar) {
  #' Estimates the IF for EY(astar) for an observation under the back-door assumption
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param outcome The observed values of the outcomes
  #' @param fit.a.c Fit of the model for the treatment A
  #' @param fit.z.a Fit of the model for the mediator Z on A
  #' @param fit.y.ac Fit of the model for the outcome Y on A and C
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the back-door assumption
  n <- length(exposure)

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a.c)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.ac)))

  model.matrix.y.astar <- data.frame(model.matrix.y)
  model.matrix.y.astar[, "a"] <- astar
  model.matrix.y.astar <- as.matrix(model.matrix.y.astar)

  y.par.hat <- summary(fit.y.ac)$coefficients[, 1]
  a.par.hat <- summary(fit.a.c)$coefficients[, 1]

  a.mean <- expit(model.matrix.a %*% a.par.hat)
  y.mean <- model.matrix.y %*% y.par.hat

  sum.c <- model.matrix.y.astar %*% y.par.hat

  psi.bd.ind <- as.numeric(exposure == astar) / (a.mean * astar + (1 - a.mean) * (1 - astar)) * (outcome - sum.c) +
    sum.c
  return(psi.bd.ind)
}

EstimateIFFrontDoor <- function(exposure, intermediate, outcome, est.pa, fit.z.a, fit.y.az, astar) {
  #' Estimates the IF for EY(astar) for an observation under the front-door assumption
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param est.pa Estimated probability of A being 1
  #' @param fit.z.a Fit of the model for the mediator Z
  #' @param fit.y.az Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the front-door assumption
  n <- length(exposure)

  sigma.z <- summary(fit.z.a)$sigma

  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.az)))

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y.az)$coefficients[, 1]
  z.par.hat <- summary(fit.z.a)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  a.mean <- est.pa # E(a) = P(A=1)
  y.mean <- model.matrix.y %*% y.par.hat


  sum.a <- as.matrix(cbind(rep(1, n), rep(a.mean, n), model.matrix.y[, "z"])) %*% y.par.hat
  sum.z <- as.matrix(cbind(rep(1, n), model.matrix.y[, "a"], z.mean_astar)) %*% y.par.hat
  sum.az <- as.matrix(cbind(rep(1, n), rep(a.mean, n), z.mean_astar)) %*% y.par.hat

  psi.fd.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / dnorm(model.matrix.y[, "z"], z.mean_ind, sigma.z)) +
    as.numeric(exposure == astar) / (a.mean * astar + (1 - a.mean) * (1 - astar)) * (sum.a - sum.az) +
    sum.z

  return(psi.fd.ind)
}

EstimateIFTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a.c, fit.z.a, fit.y.zc, astar) {
  #' Estimates the IF for EY(astar) for an observation under the two-door assumption
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a.c Fit of the model for the treatment A on C
  #' @param fit.z.a Fit of the model for the mediator Z on A
  #' @param fit.y.zc Fit of the model for the outcome Y on Z and C
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the two-door assumption
  n <- length(exposure)
  sigma.z <- summary(fit.z.a)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a.c)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.zc)))

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y.zc)$coefficients[, 1]
  z.par.hat <- summary(fit.z.a)$coefficients[, 1]
  a.par.hat <- summary(fit.a.c)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  a.mean <- expit(model.matrix.a %*% a.par.hat)
  y.mean <- model.matrix.y %*% y.par.hat
  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- y.mean
  sum.z <- as.matrix(cbind(rep(1, n), z.mean_astar, cov.vals.all)) %*% y.par.hat
  sum.az <- as.matrix(cbind(rep(1, n), z.mean_ind, cov.vals.all)) %*% y.par.hat
  psi.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / dnorm(model.matrix.y[, "z"], z.mean_ind, sigma.z)) +
    as.numeric(exposure == astar) / (a.mean * astar + (1 - a.mean) * (1 - astar)) * (sum.a - sum.az) +
    sum.z

  return(psi.td.ind)
}

EstimateIFBackTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a.c, fit.z.a, fit.y.zc, astar) {
  #' Estimates the IF for EY(astar) for an observation under the two- and the back-door assumptions
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a.c Fit of the model for the treatment A on C
  #' @param fit.z.a Fit of the model for the mediator Z
  #' @param fit.y Fit of the model for the outcome Y on Z and C
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the two- and the back-door assumptions
  n <- length(exposure)
  sigma.z <- summary(fit.z.a)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a.c)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.zc)))

  model.matrix.z0 <- data.frame(model.matrix.z)
  model.matrix.z0[, "a"] <- 0
  model.matrix.z0 <- as.matrix(model.matrix.z0)

  model.matrix.z1 <- data.frame(model.matrix.z)
  model.matrix.z1[, "a"] <- 1
  model.matrix.z1 <- as.matrix(model.matrix.z1)

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y.zc)$coefficients[, 1]
  z.par.hat <- summary(fit.z.a)$coefficients[, 1]
  a.par.hat <- summary(fit.a.c)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  a.mean <- expit(model.matrix.a %*% a.par.hat)
  y.mean <- model.matrix.y %*% y.par.hat
  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- y.mean
  sum.z <- as.matrix(cbind(rep(1, n), z.mean_astar, cov.vals.all)) %*% y.par.hat
  sum.az <- as.matrix(cbind(rep(1, n), z.mean_ind, cov.vals.all)) %*% y.par.hat
  sum.denom.az <- expit(cbind(1, cov.vals.all) %*% a.par.hat) * dnorm(model.matrix.y[, "z"], model.matrix.z1 %*% z.par.hat, sigma.z) +
    (1 - expit(cbind(1, cov.vals.all) %*% a.par.hat)) * dnorm(model.matrix.y[, "z"], model.matrix.z0 %*% z.par.hat, sigma.z)

  psi.bd.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / sum.denom.az) +
    as.numeric(exposure == astar) / (a.mean * astar + (1 - a.mean) * (1 - astar)) * (sum.a - sum.az) +
    sum.z

  return(psi.bd.td.ind)
}

EstimateIFFrontTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome,  fit.a.c, fit.z.a, fit.y.zc, astar) {
  #' Estimates the IF for EY(astar) for an observation under the front- and the two-door assumptions
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param est.pc Estimated probability of C being 1
  #' @param fit.a.c Fit of the model for the treatment A on C
  #' @param fit.z.a Fit of the model for the mediator Z
  #' @param fit.y.zc Fit of the model for the outcome Y on Z and C
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the front- and the two-door assumptions
  n <- length(exposure)
  ones <- rep(1, n)
  sigma.z <- summary(fit.z.a)$sigma
  
  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a.c)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.zc)))
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  
  
  y.par.hat <- summary(fit.y.zc)$coefficients[, 1]
  z.par.hat <- summary(fit.z.a)$coefficients[, 1]
  a.par.hat <- summary(fit.a.c)$coefficients[, 1]
  
  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  
  pa <- mean(exposure)
  
  # pa <- t((expit(c(1, 1) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 1) %*% a.par.hat)) %*% (1 - exposure)) * est.pc +
  #                (expit(c(1, 0) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 0) %*% a.par.hat)) %*% (1 - exposure)) * (1 - est.pc))
  # 
  y.mean <- model.matrix.y %*% y.par.hat
  
  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- as.matrix(cbind(ones, intermediate, matrix(colMeans(cov.vals.all), nrow = n, ncol = ncol(cov.vals.all), byrow= T))) %*% y.par.hat
  sum.z <- as.matrix(cbind(1, z.mean_astar, cov.vals.all)) %*% y.par.hat
  sum.az <- as.matrix(cbind(ones, z.mean_ind, matrix(colMeans(cov.vals.all), nrow = n, ncol = ncol(cov.vals.all), byrow= T)))  %*% y.par.hat
  psi.fd.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / dnorm(model.matrix.y[, "z"], z.mean_ind, sigma.z)) +
    as.numeric(exposure == astar) / (pa * astar + (1 - pa) * (1 - astar))  * (sum.a - sum.az) +
    sum.z
  
  
  return(psi.fd.td.ind)
}

EstimateIFBackFrontTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a.c, fit.z.a, fit.y.zc, astar) {
  #' Estimates the IF for EY(astar) for an observation under the back-, front- and the two-door assumptions
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a.c Fit of the model for the treatment A on C
  #' @param fit.z.a Fit of the model for the mediator Z on A
  #' @param fit.y.zc Fit of the model for the outcome Y on Z and C
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the back-, the front- and the two-door assumptions
  n <- length(exposure)
  ones <- rep(1, n)
  sigma.z <- summary(fit.z.a)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a.c)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.zc)))

  model.matrix.z0 <- data.frame(model.matrix.z)
  model.matrix.z0[, "a"] <- 0
  model.matrix.z0 <- as.matrix(model.matrix.z0)

  model.matrix.z1 <- data.frame(model.matrix.z)
  model.matrix.z1[, "a"] <- 1
  model.matrix.z1 <- as.matrix(model.matrix.z1)

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y.zc)$coefficients[, 1]
  z.par.hat <- summary(fit.z.a)$coefficients[, 1]
  a.par.hat <- summary(fit.a.c)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  a.mean <- expit(model.matrix.a %*% a.par.hat)
  y.mean <- model.matrix.y %*% y.par.hat

  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- as.matrix(cbind(ones, intermediate, matrix(colMeans(cov.vals.all), nrow = n, ncol = ncol(cov.vals.all), byrow= T))) %*% y.par.hat
  sum.z <- as.matrix(cbind(rep(1, n), z.mean_astar, cov.vals.all)) %*% y.par.hat
  sum.az <- as.matrix(cbind(ones, z.mean_ind, matrix(colMeans(cov.vals.all), nrow = n, ncol = ncol(cov.vals.all), byrow= T)))  %*% y.par.hat
  sum.denom.az <- expit(cbind(1, cov.vals.all) %*% a.par.hat) * dnorm(model.matrix.y[, "z"], model.matrix.z1 %*% z.par.hat, sigma.z) +
    (1 - expit(cbind(1, cov.vals.all) %*% a.par.hat)) * dnorm(model.matrix.y[, "z"], model.matrix.z0 %*% z.par.hat, sigma.z)

  # pa <- t((expit(c(1, 1) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 1) %*% a.par.hat)) %*% (1 - exposure)) * est.pc + # p(C=1)*p(A|C=1)
  #   (expit(c(1, 0) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 0) %*% a.par.hat)) %*% (1 - exposure)) * (1 - est.pc)) # p(C=0)p(A|C=0)
  
  pa <- mean(exposure)

  psi.bd.fd.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / sum.denom.az) +
    as.numeric(exposure == astar) / (pa * astar + (1 - pa) * (1 - astar)) * (sum.a - sum.az) +
    sum.z

  return(psi.bd.fd.td.ind)
}
