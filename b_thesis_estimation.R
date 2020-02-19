################################################################################
############                 TU Dortmund University                 ############
############                    Bachelor Thesis                     ############
############    Deconvolution Problems in a Psychometric Context    ############
############                     Robin Grugel                       ############
################################################################################

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tables)
library(xtable)
library(ggpubr)

########## important functions for estimation procedure ##########

# N(mean, sd^2)
cfnorm <- function(t, mean = 0, sd = 1) {
  exp(1i * mean * t - (sd^2 * t^2)/2)
}

# Log-N(meanlog, sdlog^2)
cflnorm <- function(t, meanlog = 1, sdlog = 0.5) {
  if (identical(all.equal(t, 0), TRUE)) { 
    return(1 + (0+0i))
  }
  else {
    t <- t * exp(meanlog)
    R <- integrate(function(y) exp(-log(y/t)^2/2/sdlog^2) *
                       cos(y)/y, lower = 0, upper = t)$value + 
      integrate(function(y) exp(-log(y * t)^2/2/sdlog^2) * 
                  cos(1/y)/y, lower = 0, upper = 1/t)$value
    I <- integrate(function(y) exp(-log(y/t)^2/2/sdlog^2) *
                       sin(y)/y, lower = 0, upper = t)$value +
      integrate(function(y) exp(-log(y * t)^2/2/sdlog^2) *
                       sin(1/y)/y, lower = 0, upper = 1/t)$value
    return((R + 1i * (I))/(sqrt(2 * pi) * sdlog))
  } 
}

# function: ecf() calculates the empirical characteristic function
ecf <- function(t, x) {
  sapply(t, function(t) mean(exp(1i * t * x)))
}

# function: f.K() is the kernel often used, because of the useful compact form
#           in its fourier transform
f.K <- function(x) {
  first <- ((48 * cos(x)) / (pi * x^4)) * (1 - (15/x^2))
  second <- ((144 * sin(x)) / (pi * x^5)) * (2 - (5/x^2))
  res <- first - second
  ifelse(x == 0, 0, res) # ensures, that 
}

# function: phi.k() is the fourier transform of the kernel f.K()
phi.K <- function(t) (1 - t^2)^3 * ifelse(t < -1 | t > 1, 0, 1)

# function: cf.tau() estimates the characteristic function of the random 
#           variable tau of the measurement error model X = tau + epsilon
# input: t - values of the frequency space
#        x - observed values of the variable space
#        phi.K - characteristic function (fourier transform) of the kernel
#        phi.eps - characterstic function of the error variable
# output: estimated function values of the characteristic function of tau
phi.tau <- function(t, y, b, kern = phi.K, error = cfnorm, ...) {
  sapply(t, function(t) 
    (mean(exp(1i * t * y)) * kern(t * b))  / error(t, mean = 0, ...))
}

# function: lnorm.var() calculates the variance of the log-normal distribution
lnorm.var <- function(logmean, logsd) {
  (exp(logsd^2) - 1) * exp(2 * logmean + logsd^2)
}

# function: lnorm.exp() calculates the expected value of the log-normal 
#           distribution
lnorm.exp <- function(logmean, logsd) {
  exp(logmean + (logsd^2)/2)
}

# function: boot.bw() selects the optimal bandwidth by the algorithm 
#           of delaigle (2004a, 2004b)
# input: x - contaminated observations
#        sig - square root of the error variance
#        K - characteristic function of the kernel
#        error - characteristic function of the measurement error
#        rot - should the rule of thumb be applied to estimate the pilot bw
# output: bw - estimate for the optimal bandwidth, that minimizes the MISE
boot.bw <- function(x, sig, error = cfnorm, rot = FALSE, subdivs = 1000, upperh = 10) {
  n <- length(x)
  
  # using rule of thumb as the pilot bandwidth
  if (rot) {
    g2 <- sqrt(2) * sig * ((log(n) / 2) ^ (-0.5))
  } else {
  # calculate second moment for the kernel function
  # mu2 <- integrate(function(x) (x^2) * f.K(x), -Inf, Inf, 
  #                  subdivisions = 100000000)$value
  mu2 <- 6
  
  # asymptotic bias for R(f.hat_x^(r)(.;g_r))
  ABias <- function(g, R, r, error) {
    abias.int.fun <- 
      function(t) (t^(2*r)) * (abs(phi.K(t))^2) * (abs(error(t / g, sd = sig))^(-2))
    abias.int <-  integrate(abias.int.fun, -1, 1, subdivisions = subdivs)$value
    -(g^2) * mu2 * R + ((2 * pi * n * (g^(2*r + 1)))^(-1)) * abias.int
  }
  
  # normal reference method (assuming a parametric model for the true distribution)
  # R4 <- (factorial(8) * sqrt(var(x) - sig^2)^(-9)) / ((2^9) * factorial(4) * sqrt(pi))
  R4 <- (40320 * sqrt(var(x) - sig^2)^(-9)) / (12288 * sqrt(pi))
  
  ecf <- function(t) sapply(t, function(t) mean(exp(1i * t * x)))
  
  # estimation of R(f.hat_x^(r)(.;g_r)) like in delaigle (2004b)
  R.hat <- function(x, r, g) {
    R.hat.int.fun <- function(t) (t^(2*r)) * (abs(ecf(t / g))^2) * (abs(phi.K(t))^2) * (abs(error(t / g, sd = sig))^(-2))
    R.hat.int <- integrate(R.hat.int.fun, -1, 1, subdivisions = subdivs)$value
    (1 / (2 * pi * g^(2*r + 1))) * R.hat.int
  }
  
  # step 1
  g3 <- optimize(function(g) abs(ABias(g, R4, 3, cfnorm)), c(0, upperh))$minimum
  R3 <- (1 / (2 * pi * g3^(2*3 + 1))) * integrate(function(t) (t^(2*3)) * (abs(ecf(t / g3))^2) * (abs(phi.K(t))^2) * (abs(error(t / g3, sd = sig))^(-2)), lower = -1, upper = 1)$value
  R3 <- R.hat(x, 3, g3)
  
  # step 2
  g2 <- optimize(function(g) abs(ABias(g, R3, 2, cfnorm)), c(0, upperh))$minimum
  }
  
  # estimation of the MISE_2(h) from delaigle (2004a)
  MISE2 <- function(h) {
    est_cf <- function(t) phi.tau(t, y = x, b = g2, phi.K, sd = sig)
    
    fac1 <- 1 / (2 * pi * n * h)
    fac2 <- (1 - n ^ (-1)) * ((2 * pi) ^ (-1))
    fac3 <- 1 / pi
    
    int1.fun <- function(t) (abs(phi.K(t))^2) / (abs(error(t / h, sd = sig))^2)
    int1 <- integrate(int1.fun, -1, 1, subdivisions = subdivs)$value
    int2.fun <- function(t) (abs(est_cf(t))^2) * (abs(phi.K(h * t))^2)
    int2 <- integrate(int2.fun, -1, 1, subdivisions = subdivs)$value
    int3.fun <- function(t) (abs(est_cf(t))^2) * phi.K(h * t)
    int3 <- integrate(int3.fun, -1, 1, subdivisions = subdivs)$value
    
    fac1 * int1 + fac2 * int2 - fac3 * int3
  }
  
  bw <- optimize(MISE2,  c(0, upperh))$minimum
  
  return(bw)
}

# function: f.tau() is the deconvoluting kernel density estimation
f.tau <- function(X, sig, rot.pilot = FALSE, error = cfnorm) {
  h <- boot.bw(X, sig, rot.pilot, error = error)
  out <- deconvolve(FUN = function(t) 
    phi.tau(t = t, y = X, b = h, kern = phi.K, error = error, sd = sig),
    a = -1 / h,
    b = 1 / h,
    c = min(X) - sd(X),
    d = max(X) + sd(X)
  )
  out$eval[which(Re(out$eval) < 0)] <- 0 + 0i
  
  out
}

# function: deconvolve() calculates a fourier integral numerically by 
#           approximation through the DFT of the integral scaled and adjusted
#           to concern asymmetric integration and varying range in the
#           transform variable space as described by Inverarity (2002)
# input: FUN - function of which
#        a - lower integration limit
#        b - upper integration limit
#        c - lower limit for evaluated values in transform variable space
#        d - upper limit for evaluated values in transform variable space
#        m - resolution of the transform (ideally 2^k for k > 8)
# output: w - values in transform variable space
#         eval - evaluated transform in w
deconvolve <- function(FUN, a, b, c, d, m = 2^12) {
  res <- list()
  
  beta <- (b - a) / m
  gamma <- (d - c) / m
  delta <- beta * gamma / 2
  j1 <- seq(0, m - 1, length.out = m)
  j2 <- seq(m, 2 * m - 1, length.out = m)
  t <- a + beta * j1
  f_t <- FUN(t)
  res$w <- w <- c + gamma * 0:(m - 1)
  
  y <- complex(2 * m)
  y[1:m] <- f_t * exp(-1i * j1 * (beta * c + delta * j1))
  
  z <- complex(2 * m)
  z[1:m] <- exp(delta * (j1 ^ 2) * 1i)
  z[(m + 1):(2 * m)] <- exp(delta * ((j2 - 2 * m) ^ 2) * 1i)
  
  # another algorithm for the fft might be an improvement
  fft_k <- (fft(fft(y) * fft(z), inverse = TRUE) / length(y))[1:m]
  
  # factor to secure, that the inverse transform deconvolves to a valid
  # density
  fac <- (2 * pi) ^ (-1)
  
  res$eval <- fac * beta * exp(-1i * (a * w + delta * j1 ^ 2)) * fft_k
  
  return(res)
}

# function: cronbach() calculates cronbach's alpha
cronbach <- function(x) {
  k <- ncol(x)
  sum.var <- sum(apply(x, 2, var, na.rm = TRUE))
  x.var <- var(rowSums(x))
  (k / (k - 1)) * (1 - sum.var / x.var)
}

# function: split.var.weight()
split.var.weight <- function(m = 4, big = 4/5){
  res.weight <- rep(big, m - 1) * (1 - big)^(0:(m-2))
  res.weight <- c(res.weight, 1 - sum(res.weight))
  res.weight
}