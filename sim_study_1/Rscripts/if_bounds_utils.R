## Functions for simulations

# data generation and additional functions --------------------------------
expit <- function(x) {
  exp(x) / (1 + exp(x))
}
p_c <- function(c_value, pc) { #
  #' Calculates the probability of a Bernoulli(pc) random variable to be equal to c_value
  #' @param c_value Value of a random variable
  #' @param pc Parameter of the distribution
  #'
  #' @return Probability of a Bern(pc) random variable to be equal to c_value
  pc^c_value * (1 - pc)^(1 - c_value)
}
pa_given_c <- function(a_value, c_value, par.a.c) {
  #' Calculates p(A = a_value|C = c_value)
  #' @param a_value Value of A
  #' @param c_value Value of C
  #' @param par.a.c Parameters of the distribution of A
  #'
  #' @return Probability p(A = a_value|C = c_value)
  p1 <- expit(par.a.c %*% c(1, c_value))
  p1^a_value * (1 - p1)^(1 - a_value)
}
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

generate_data_ca_binary_zy_cont <- function(n, pc, par.a.c, par.z.a, par.y.zc, sigma.z, sigma.y) {
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
  c <- rbinom(n, 1, pc) # Generate the confounder
  a <- rbinom(n, 1, expit(cbind(rep(1, n), c) %*% par.a.c)) # Generate the exposure (given c)
  z <- rnorm(n, cbind(rep(1, n), a) %*% par.z.a, sigma.z) # Generate  the mediator (given a and c)
  y <- rnorm(n, cbind(rep(1, n), z, c) %*% par.y.zc, sigma.y) # Generate the outcome Y (given a, m, c)
  sim.data <- data.frame(cbind(y, a, z, c))
  return(sim.data)
}

# Calculate true mean of potential outcome --------------------------------
calculate_mean_potential_outcome <- function(astar, pc, par.z.a, par.y.zc) {
  # The true value of EY(a*, Z(a*)) = int_{c,z} E(Y|Z=z,C=c)P(Z=z|A=a*, C=c)P(C=c)dcdz
  # according to the two-door adjustment, which is true in our simulated data.
  # E(a*, Z(a*)) = int_{c,z} E(Y|Z=z,C=c)P(Z=z|A=a*, C=c)P(C=c)dcdz =
  # = int_{c,z} E(Y|Z=z,C=c)P(Z=z|A=a*)P(C=c)dcdz = par.y.zc[1] + int_{c,z} par.y.zc[2]zP(Z=z|A=a*)P(C=c)dcdz
  # + int_{c,z} par.y.zc[3]cP(Z=z|A=a*)P(C=c)dcdz =
  # par.y.zc[1] + par.y.zc[2]*E(Z=z|A=a*) + par.y.zc[3]*E(C)
  #' Calculates EY(astar)
  #' @param astar The value of the treatment
  #' @param pc Parameter of the distribution of C
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #'
  #' @return  EY(astar)
  mean.z.given.astar <- c(1, astar) %*% par.z.a
  psi <- c(1, mean.z.given.astar, pc) %*% par.y.zc
  return(psi)
}


# Calculate the bounds ----------------------------------------------------
CalculateBoundBackDoor <- function(astar, a, pc, par.a.c, par.y.zc, sigma.y, sigma.z) {
  #' Calculates the back-door bound for the estimation of EY(astar) - EY(a)
  #' @param astar The value of the treatment for treated
  #' @param a The value of the treatment for non-treated
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return The back-door bound for the estimation of EY(astar) - EY(a)
  (sigma.y^2 + sigma.z^2 * par.y.zc[2]^2) *
    (p_c(1, pc) / pa_given_c(astar, c_value = 1, par.a.c) + p_c(1, pc) / pa_given_c(a, c_value = 1, par.a.c) +
      p_c(0, pc) / pa_given_c(astar, c_value = 0, par.a.c) + p_c(0, pc) / pa_given_c(a, c_value = 0, par.a.c))
}

CalculateBoundFrontDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z) {
  #' Calculates the front-door bound for the estimation of EY(astar) - EY(a)
  #' @param astar The value of the treatment for treated
  #' @param a The value of the treatment for non-treated
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return The front-door bound for the estimation of EY(astar) - EY(a)
  fa <- function(a_value) {
    pa_given_c(a_value = a_value, c_value = 0, par.a.c = par.a.c) * p_c(c_value = 0, pc = pc) +
      pa_given_c(a_value = a_value, c_value = 1, par.a.c = par.a.c) * p_c(c_value = 1, pc = pc)
  }

  EY_az <- function(a_value, z_value, par.a.c) {
    Ec_given_a <- function(a_value) {
      pa_given_c(a_value = a_value, c = 1, par.a.c = par.a.c) * pc / fa(a_value = a_value)
    }
    par.y.zc[1] + z_value * par.y.zc[2] + Ec_given_a(a_value) * par.y.zc[3]
  }

  fa_varY_az <- function(a_value, z_value) {
    fa(a_value) * (sigma.y^2 + (par.y.zc[1] + par.y.zc[2] * z_value)^2) +
      pa_given_c(a_value = a_value, c_value = 1, par.a.c) * p_c(c_value = 1, pc = pc) *
        (par.y.zc[3]^2 + 2 * par.y.zc[3] * (par.y.zc[1] + par.y.zc[2] * z_value)) - fa(a_value) *
        EY_az(a_value = a_value, z_value = z_value, par.a.c = par.a.c)^2
  }
  EZastar <- function(a_value) {
    c(1, a_value) %*% par.z.a
  }
  integrate(Vectorize(function(z) {
    fa_varY_az(a_value = astar, z_value = z) * pz_given_a(z_value = z, astar, par.z.a, sigma.z) + # sum_z f(z|a^*)*fa_varY_az(a*) + f(z|a)*fa_varY_az(a)
      fa_varY_az(a_value = a, z_value = z) * pz_given_a(z_value = z, a, par.z.a, sigma.z)
  }), -Inf, Inf, rel.tol = 1e-10)$value -
    2 * integrate(Vectorize(function(z) {
      fa_varY_az(a_value = astar, z_value = z) * pz_given_a(z_value = z, a_value = a, par.z.a, sigma.z) + # sum_z  f(z|a)*fa_varY_az(a*) + f(z|a*)*fa_varY_az(a)
        fa_varY_az(a_value = a, z_value = z) * pz_given_a(z_value = z, a_value = astar, par.z.a, sigma.z)
    }), -Inf, Inf, rel.tol = 1e-10)$value +
    integrate(Vectorize(function(z) {
      fa_varY_az(a_value = astar, z_value = z) / sqrt(2 * pi) / sigma.z * exp(-2 * (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) +
        fa_varY_az(a_value = a, z_value = z) / sqrt(2 * pi) / sigma.z * exp(-2 * (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))
    }), -Inf, Inf, rel.tol = 1e-10)$value + # f^2(z|a)/f(z|a*)*fa_varY_az(a*) + f^2(z|a*)/f(z|a)*fa_varY_az(a)
    (par.y.zc[1] + par.y.zc[3] * pc)^2 * (1 / fa(a) + 1 / fa(astar)) + # sum_z f(z|a*)()^2/f(a*)) + same for a
    par.y.zc[2]^2 * (
      (sigma.z^2 + EZastar(astar)^2) / fa(astar) + (sigma.z^2 + EZastar(a)^2) / fa(a)
    ) +
    2 * par.y.zc[2] * (par.y.zc[1] + par.y.zc[3] * pc) * (EZastar(astar) / fa(astar) + EZastar(a) / fa(a)) -
    (c(1, EZastar(astar), pc) %*% par.y.zc)^2 / fa(astar) - #  - E^2Y(a*)/f(a*)
    (c(1, EZastar(a), pc) %*% par.y.zc)^2 / fa(a) + #  - E^2Y(a)/f(a)
    par.y.zc[2]^2 * (EZastar(astar) - EZastar(a))^2 - # + sum_a f(a)* (sum_z EY|a,z (f(z|a*)- f(z|a)))^2
    (calculate_mean_potential_outcome(astar, pc, par.z.a, par.y.zc) - calculate_mean_potential_outcome(a, pc, par.z.a, par.y.zc))^2 #  - ATE^2
}

CalculateBoundTwoDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z){
  #' Calculates the two-door bound for estimation of EY(astar) - EY(a)
  #' @param astar The value of the treatment for treated
  #' @param a The value of the treatment for non-treated
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return The two-door bound for the estimation of EY(astar) - EY(a)
  sigma.y^2 * (integrate(Vectorize(function(z){1/sqrt(2*pi)/sigma.z * exp(-2*(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2)+(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2))}), -Inf, Inf, rel.tol =1e-10)$value-1) * 
    (p_c(c_value=0, pc) * pa_given_c(a_value=a, c_value=0, par.a.c) + p_c(c_value=1, pc) * pa_given_c(a_value=a, c_value=1, par.a.c)) + 
  sigma.y^2 * (integrate(Vectorize(function(z){1/sqrt(2*pi)/sigma.z * exp(-2*(z-as.numeric(c(1,a)%*%par.z.a))^2/(2*sigma.z^2)+(z-as.numeric(c(1,astar)%*%par.z.a))^2/(2*sigma.z^2))}), -Inf, Inf, rel.tol =1e-10)$value-1)  * 
    (p_c(c_value=0, pc) * pa_given_c(a_value=astar, c_value=0, par.a.c) + p_c(c_value=1, pc) * pa_given_c(a_value=astar, c_value=1, par.a.c)) -
  sigma.y^2 * p_c(c_value=0, pc)*(1/pa_given_c(a_value=astar, c_value=0, par.a.c) + 1/pa_given_c(a_value=a, c_value=0, par.a.c)) -
  sigma.y^2 * p_c(c_value=1, pc)*(1/pa_given_c(a_value=astar, c_value=1, par.a.c) + 1/pa_given_c(a_value=a, c_value=1, par.a.c)) + 
    CalculateBoundBackDoor(astar=astar, a=a, pc=pc, par.a.c=par.a.c, par.y.zc=par.y.zc, sigma.y=sigma.y, sigma.z = sigma.z) 
}

CalculateBoundBackTwoDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z) {
  #' Calculates the bound for the estimation of EY(astar) - EY(a) under the back- and the two-door assumptions
  #' @param astar The value of the treatment for treated
  #' @param a The value of the treatment for non-treated
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return The bound for the estimation of EY(astar) - EY(a) under the back- and the two-door assumptions
  sigma.y^2 * p_c(c_value = 0, pc) * integrate(Vectorize(function(z) {
    1 / sqrt(2 * pi * sigma.z^2) * (
      exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) +
        exp(-2 * (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, 1) %*% par.z.a))^2 / (2 * sigma.z^2)) -
        2 * exp(-(z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))) /
      (pa_given_c(a_value = astar, c_value = 0, par.a.c) +
        pa_given_c(a_value = a, c_value = 0, par.a.c) * exp(-(z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2))
      )
  }), -Inf, Inf, rel.tol = 1e-10)$value +
    sigma.y^2 * p_c(c_value = 1, pc) * integrate(Vectorize(function(z) {
      1 / sqrt(2 * pi * sigma.z^2) * (
        exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) +
          exp(-2 * (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, 1) %*% par.z.a))^2 / (2 * sigma.z^2)) -
          2 * exp(-(z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))) /
        (pa_given_c(a_value = astar, c_value = 1, par.a.c) +
          pa_given_c(a_value = a, c_value = 1, par.a.c) * exp(-(z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2))
        )
    }), -Inf, Inf, rel.tol = 1e-10)$value -
    sigma.y^2 * p_c(c_value = 0, pc) * (1 / pa_given_c(a_value = astar, c_value = 0, par.a.c) + 1 / pa_given_c(a_value = a, c_value = 0, par.a.c)) -
    sigma.y^2 * p_c(c_value = 1, pc) * (1 / pa_given_c(a_value = astar, c_value = 1, par.a.c) + 1 / pa_given_c(a_value = a, c_value = 1, par.a.c)) +
    CalculateBoundBackDoor(astar = astar, a = a, pc = pc, par.a.c = par.a.c, par.y.zc = par.y.zc, sigma.y = sigma.y, sigma.z = sigma.z)
}

CalculateBoundFrontTwoDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z) {
  #' Calculates the bound for the estimation of EY(astar) - EY(a) under the front- and the two-door assumptions
  #' @param astar The value of the treatment for treated
  #' @param a The value of the treatment for non-treated
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return The bound for the estimation of EY(astar) - EY(a) under the front- and the two-door assumptions
  fa <- function(a_value) {
    pa_given_c(a_value = a_value, c_value = 0, par.a.c = par.a.c) * p_c(c_value = 0, pc = pc) +
      pa_given_c(a_value = a_value, c_value = 1, par.a.c = par.a.c) * p_c(c_value = 1, pc = pc)
  }
  EZastar <- function(a_value) {
    c(1, a_value) %*% par.z.a
  }
  sigma.y^2 * (integrate(
    function(z) {
      1 / sqrt(2 * pi) / sigma.z * exp(-2 * (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))
    },
    -Inf, Inf,
    rel.tol = 1e-10
  )$value - 1) *
    (p_c(c_value = 0, pc) * pa_given_c(a_value = a, c_value = 0, par.a.c) + p_c(c_value = 1, pc) * pa_given_c(a_value = a, c_value = 1, par.a.c)) +
    sigma.y^2 * (integrate(
      function(z) {
        1 / sqrt(2 * pi) / sigma.z * exp(-2 * (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) + (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2))
      },
      -Inf, Inf,
      rel.tol = 1e-10
    )$value - 1) *
      (p_c(c_value = 0, pc) * pa_given_c(a_value = astar, c_value = 0, par.a.c) + p_c(c_value = 1, pc) * pa_given_c(a_value = astar, c_value = 1, par.a.c)) + # similar to CalculateBoundTwoDoor
    (par.y.zc[1] + par.y.zc[3] * pc)^2 * (1 / fa(a) + 1 / fa(astar)) + # sum_z f(z|a*)()^2/f(a*)) + same for a, as in CalculateBoundFrontDoor
    par.y.zc[2]^2 * ((sigma.z^2 + EZastar(astar)^2) / fa(astar) + (sigma.z^2 + EZastar(a)^2) / fa(a)
    ) +
    2 * par.y.zc[2] * (par.y.zc[1] + par.y.zc[3] * pc) * (EZastar(astar) / fa(astar) + EZastar(a) / fa(a)) -
    (c(1, EZastar(astar), pc) %*% par.y.zc)^2 / fa(astar) - #  - E^2Y(a*)/f(a*), as in CalculateBoundFrontDoor
    (c(1, EZastar(a), pc) %*% par.y.zc)^2 / fa(a) #  - E^2Y(a)/f(a), as in CalculateBoundFrontDoor
  # the last term is 0 since E[Y|A,z,c] = E[Y|z,c] in the considered data distribution
}

CalculateBoundBackFrontTwoDoor <- function(astar, a, pc, par.a.c, par.z.a, par.y.zc, sigma.y, sigma.z) {
  #' Calculates the bound for the estimation of EY(astar) - EY(a) under the back, the front- and the two-door assumptions
  #' @param astar The value of the treatment for treated
  #' @param a The value of the treatment for non-treated
  #' @param pc Parameter of the distribution of C
  #' @param par.a.c Parameters of the distribution of A
  #' @param par.z.a Parameters of the mean of Z(a_value)|A
  #' @param par.y.zc Parameters of the mean of Y(a_value)|Z, C
  #' @param sigma.z Standard deviation of Z(a_value)|A
  #' @param sigma.y Standard deviation of Y(a_value)|Z,C
  #'
  #' @return The bound for the estimation of EY(astar) - EY(a) under the back, the front- and the two-door assumptions
  fa <- function(a_value) {
    pa_given_c(a_value = a_value, c_value = 0, par.a.c = par.a.c) * p_c(c_value = 0, pc = pc) +
      pa_given_c(a_value = a_value, c_value = 1, par.a.c = par.a.c) * p_c(c_value = 1, pc = pc)
  }
  EZastar <- function(a_value) {
    c(1, a_value) %*% par.z.a
  }

  sigma.y^2 * p_c(c_value = 0, pc) / pa_given_c(a_value = 0, c_value = 0, par.a.c) *
    (2 - 4 * integrate(Vectorize(function(z) {
      1 / sqrt(2 * pi) / sigma.z * exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) /
        (exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2) +
          (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2)) + 1
        )
    }), -Inf, Inf, rel.tol = 1e-10)$value) +
    sigma.y^2 * p_c(c_value = 1, pc) * (
      integrate(Vectorize(function(z) {
        1 / sqrt(2 * pi) / sigma.z * exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) /
          (pa_given_c(a_value = 0, c_value = 1, par.a.c) *
            exp(-(z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2) +
              (z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) +
            pa_given_c(a_value = 1, c_value = 1, par.a.c)
          )
      }), -Inf, Inf, rel.tol = 1e-10)$value +
        integrate(Vectorize(function(z) {
          1 / sqrt(2 * pi) / sigma.z * exp(-(z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2)) /
            (pa_given_c(a_value = 0, c_value = 1, par.a.c) +
              pa_given_c(a_value = 1, c_value = 1, par.a.c) * exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2) +
                (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))
            )
        }), -Inf, Inf, rel.tol = 1e-10)$value -
        2 * integrate(Vectorize(function(z) {
          1 / sqrt(2 * pi) / sigma.z * exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2)) /
            (pa_given_c(a_value = 0, c_value = 1, par.a.c) +
              pa_given_c(a_value = 1, c_value = 1, par.a.c) * exp(-(z - as.numeric(c(1, astar) %*% par.z.a))^2 / (2 * sigma.z^2) +
                (z - as.numeric(c(1, a) %*% par.z.a))^2 / (2 * sigma.z^2))
            )
        }), -Inf, Inf, rel.tol = 1e-10)$value
    ) +
    (par.y.zc[1] + par.y.zc[3] * pc)^2 * (1 / fa(a) + 1 / fa(astar)) + # sum_z f(z|a*)()^2/f(a*)) + same for a, as in CalculateBoundFrontDoor
    par.y.zc[2]^2 * ((sigma.z^2 + EZastar(astar)^2) / fa(astar) + (sigma.z^2 + EZastar(a)^2) / fa(a)
    ) +
    2 * par.y.zc[2] * (par.y.zc[1] + par.y.zc[3] * pc) * (EZastar(astar) / fa(astar) + EZastar(a) / fa(a)) -
    (c(1, EZastar(astar), pc) %*% par.y.zc)^2 / fa(astar) - #  - E^2Y(a*)/f(a*), as in CalculateBoundFrontDoor
    (c(1, EZastar(a), pc) %*% par.y.zc)^2 / fa(a) #  - E^2Y(a)/f(a), as in CalculateBoundFrontDoor
  # the last term is 0 since E[Y|A,z,c] = E[Y|z,c] in the considered data distribution
}


# Estimate influence functions for EY(astar)------------------------------------
EstimateIFBackDoor <- function(cov.vals.all, exposure, outcome, fit.a, fit.z, fit.y.ac, astar) {
  #' Estimates the IF for EY(astar) for an observation under the back-door assumption
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param outcome The observed values of the outcomes
  #' @param fit.a Fit of the model for the treatment A
  #' @param fit.z Fit of the model for the mediator Z
  #' @param fit.y.ac Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the back-door assumption
  n <- length(exposure)

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.ac)))

  model.matrix.y.astar <- data.frame(model.matrix.y)
  model.matrix.y.astar[, "a"] <- astar
  model.matrix.y.astar <- as.matrix(model.matrix.y.astar)

  y.par.hat <- summary(fit.y.ac)$coefficients[, 1]
  a.par.hat <- summary(fit.a)$coefficients[, 1]

  a.mean <- expit(model.matrix.a %*% a.par.hat)
  y.mean <- model.matrix.y %*% y.par.hat

  sum.c <- model.matrix.y.astar %*% y.par.hat

  psi.bd.ind <- as.numeric(exposure == astar) / (a.mean * astar + (1 - a.mean) * (1 - astar)) * (outcome - sum.c) +
    sum.c
  return(psi.bd.ind)
}

EstimateIFFrontDoor <- function(exposure, intermediate, outcome, fit.z, fit.y.az, astar) {
  #' Estimates the IF for EY(astar) for an observation under the front-door assumption
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.z Fit of the model for the mediator Z
  #' @param fit.y.az Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the front-door assumption
  n <- length(exposure)

  sigma.z <- summary(fit.z)$sigma

  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y.az)))

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y.az)$coefficients[, 1]
  z.par.hat <- summary(fit.z)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  a.mean <- mean(exposure) # E(a) = P(A=1)
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

EstimateIFTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar) {
  #' Estimates the IF for EY(astar) for an observation under the two-door assumption
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a Fit of the model for the treatment A
  #' @param fit.z Fit of the model for the mediator Z
  #' @param fit.y Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the two-door assumption
  n <- length(exposure)
  sigma.z <- summary(fit.z)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y)$coefficients[, 1]
  z.par.hat <- summary(fit.z)$coefficients[, 1]
  a.par.hat <- summary(fit.a)$coefficients[, 1]

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

EstimateIFTwoBackDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar) {
  #' Estimates the IF for EY(astar) for an observation under the two- and the back-door assumptions
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a Fit of the model for the treatment A
  #' @param fit.z Fit of the model for the mediator Z
  #' @param fit.y Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the two- and the back-door assumptions
  n <- length(exposure)
  sigma.z <- summary(fit.z)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))

  model.matrix.z0 <- data.frame(model.matrix.z)
  model.matrix.z0[, "a"] <- 0
  model.matrix.z0 <- as.matrix(model.matrix.z0)

  model.matrix.z1 <- data.frame(model.matrix.z)
  model.matrix.z1[, "a"] <- 1
  model.matrix.z1 <- as.matrix(model.matrix.z1)

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y)$coefficients[, 1]
  z.par.hat <- summary(fit.z)$coefficients[, 1]
  a.par.hat <- summary(fit.a)$coefficients[, 1]

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

  psi.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / sum.denom.az) +
    as.numeric(exposure == astar) / (a.mean * astar + (1 - a.mean) * (1 - astar)) * (sum.a - sum.az) +
    sum.z
  # psi.td.hat <- mean(psi.td.ind)

  return(psi.td.ind)
}

EstimateIFFrontTwoDoor <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar) {
  #' Estimates the IF for EY(astar) for an observation under the front- and the two-door assumptions
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a Fit of the model for the treatment A
  #' @param fit.z Fit of the model for the mediator Z
  #' @param fit.y Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the front- and the two-door assumptions
  n <- length(exposure)
  n <- length(exposure)
  ones <- rep(1, n)
  sigma.z <- summary(fit.z)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y)$coefficients[, 1]
  z.par.hat <- summary(fit.z)$coefficients[, 1]
  a.par.hat <- summary(fit.a)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat

  a.mean <- expit(model.matrix.a %*% a.par.hat)

  pa <- t((expit(c(1, 1) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 1) %*% a.par.hat)) %*% (1 - exposure)) * mean(cov.vals.all) +
    (expit(c(1, 0) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 0) %*% a.par.hat)) %*% (1 - exposure)) * (1 - mean(cov.vals.all)))

  y.mean <- model.matrix.y %*% y.par.hat

  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- model.matrix.y[, c("X.Intercept.", "z")] %*% y.par.hat[c("(Intercept)", "z")] +
    y.par.hat["c"] * mean(cov.vals.all)
  sum.z <- as.matrix(cbind(ones, z.mean_astar, cov.vals.all)) %*% y.par.hat
  sum.az <- as.matrix(cbind(ones, z.mean_ind, mean(cov.vals.all))) %*% y.par.hat
  psi.fd.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / dnorm(model.matrix.y[, "z"], z.mean_ind, sigma.z)) +
    as.numeric(exposure == astar) / pa * (sum.a - sum.az) +
    sum.z
  
  return(psi.fd.td.ind)
}

EstimateIFAll <- function(cov.vals.all, exposure, intermediate, outcome, fit.a, fit.z, fit.y, astar) {
  #' Estimates the IF for EY(astar) for an observation under the back-, front- and the two-door assumptions
  #' @param cov.vals.all The observed values of the confounders
  #' @param exposure The observed values of the treatment
  #' @param intermediate The observed values of the mediator
  #' @param outcome The observed values of the outcomes
  #' @param fit.a Fit of the model for the treatment A
  #' @param fit.z Fit of the model for the mediator Z
  #' @param fit.y Fit of the model for the outcome Y
  #' @param astar Treatment value
  #'
  #' @return The estimate of the IF for EY(astar) for an observation under the back-, the front- and the two-door assumptions
  n <- length(exposure)
  sigma.z <- summary(fit.z)$sigma

  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))

  model.matrix.z0 <- data.frame(model.matrix.z)
  model.matrix.z0[, "a"] <- 0
  model.matrix.z0 <- as.matrix(model.matrix.z0)

  model.matrix.z1 <- data.frame(model.matrix.z)
  model.matrix.z1[, "a"] <- 1
  model.matrix.z1 <- as.matrix(model.matrix.z1)

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[, "a"] <- astar
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)

  y.par.hat <- summary(fit.y)$coefficients[, 1]
  z.par.hat <- summary(fit.z)$coefficients[, 1]
  a.par.hat <- summary(fit.a)$coefficients[, 1]

  z.mean_astar <- model.matrix.z_astar %*% z.par.hat
  z.mean_ind <- model.matrix.z %*% z.par.hat
  a.mean <- expit(model.matrix.a %*% a.par.hat)
  y.mean <- model.matrix.y %*% y.par.hat

  # E(Y|a,Z_i,C_i)=E(Y|Z_i, C_i)
  sum.a <- model.matrix.y[, c("X.Intercept.", "z")] %*% y.par.hat[c("(Intercept)", "z")] + mean(cov.vals.all) * y.par.hat["c"]
  sum.z <- as.matrix(cbind(rep(1, n), z.mean_astar, cov.vals.all)) %*% y.par.hat
  sum.az <- as.matrix(cbind(rep(1, n), z.mean_ind, mean(cov.vals.all))) %*% y.par.hat
  sum.denom.az <- expit(cbind(1, cov.vals.all) %*% a.par.hat) * dnorm(model.matrix.y[, "z"], model.matrix.z1 %*% z.par.hat, sigma.z) +
    (1 - expit(cbind(1, cov.vals.all) %*% a.par.hat)) * dnorm(model.matrix.y[, "z"], model.matrix.z0 %*% z.par.hat, sigma.z)

  pa <- t((expit(c(1, 1) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 1) %*% a.par.hat)) %*% (1 - exposure)) * mean(cov.vals.all) + # p(c=1)*p(A|c=1)
    (expit(c(1, 0) %*% a.par.hat) %*% exposure + (1 - expit(c(1, 0) %*% a.par.hat)) %*% (1 - exposure)) * (1 - mean(cov.vals.all))) # p(c=0)p(A|c=0)

  psi.td.ind <- (outcome - y.mean) *
    (dnorm(model.matrix.y[, "z"], z.mean_astar, sigma.z) / sum.denom.az) +
    as.numeric(exposure == astar) / pa * (sum.a - sum.az) +
    sum.z

  return(psi.td.ind)
}
