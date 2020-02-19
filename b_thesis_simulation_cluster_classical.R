################################################################################
############                 TU Dortmund University                 ############
############                    Bachelor Thesis                     ############
############    Deconvolution Problems in a Psychometric Context    ############
############                     Robin Grugel                       ############
################################################################################

# Caution: this version of the code already is modified to obtain the 
#          distribution of the simulated ISE to also evaluate standard errors,
#          which isn't part of the actual thesis.

source("b_thesis_estimation.R")
library(parallel)

# function: calc.ISE() calculates the integrated scared error
# input: n - sample size
#        no_var - number of indicators
#        true_rel - true reliability
#        essentially - specifies, if models describes essential variables
#        tau_congeneric - specifies, if model is tau-congeneric 
#        tau_equivalent - specifies, if model is tau-equivalent 
calc.ISE <- function(n, no_var, true_rel, dis, true_var, essentially,
                     tau_congeneric, tau_equivalent, ...) {
  parts <- 1 / rep(no_var, no_var)
  lambda <- rep(1, no_var)
  translation <- rep(0, no_var)
  
  if (tau_equivalent) {
    parts <- split.var.weight(m = no_var)
  }
  
  if (tau_congeneric) {
    lambda <- switch(
      log(no_var) / log(2),
      c(0.91, 0.85),
      c(0.95, 0.91, 0.86, 0.89),
      c(0.89, 0.93, 0.87, 0.99, 0.85, 0.91, 0.92, 0.89)
    )
  }
  
  if (essentially) {
    translation <- switch(
      log(no_var) / log(2),
      c(0.27,-0.27), 
      c(0.31, -0.09, -0.12, -0.1),
      c(0.1, -0.68, 0.57, -0.2, 0.05, 0.3, 0.13, -0.27)
    )
  }
  
  tau <- do.call(paste0("r", dis), list(n = n, ...))
  true.f <- function(x)
    do.call(paste0("d", dis), list(x = x, ...))
  
  sum_sigma_eps <- 
    true_var * (mean(lambda) ^ 2) * (no_var ^ 2) * (((1 - true_rel) / true_rel))
  sigma_eps <- parts * sum_sigma_eps
  add.trans <-
    matrix(translation,
           nrow = n,
           ncol = no_var,
           byrow = TRUE)
  multiply.load <-
    matrix(lambda,
           nrow = n,
           ncol = no_var,
           byrow = TRUE)
  add.eps <- sapply(sigma_eps, function(x) rnorm(n, mean = 0, sd = sqrt(x)))
  x <- ((replicate(no_var, tau) + add.trans) * multiply.load) + add.eps
  
  x.mean <- rowMeans(x)
  out.kde <- density(x.mean, bw = "nrd")
  
  # evaluate ISE by approximated functions
  est.f <- approxfun(out.kde$x, out.kde$y)
  SE <- function(x) return((est.f(x) - true.f(x)) ^ 2)
  ISE <- integrate(SE, min(out.kde$x), max(out.kde$x), 
                   subdivisions = 100000)$value
  
  ISE
}

# function: sim.MISE() simulates distribution of ISE to estimate MISE for N(1,4)
#           logN(1,0.25) distributed true scores
# input: seq_n - sequence of sample sizes
#        seq_no_var - sequence of number of indicators
#        seq_rel - sequence of true reliability
#        seq_model - sequence of desired models
sim.MISE <- function(cl, no_sim, seq_n, seq_no_var, seq_rel, seq_model) {
  
  par.grid <- expand.grid(list(n = seq_n, no_var = seq_no_var, 
                               true_rel = seq_rel, tau = seq_model)) 
  

  clusterExport(cl, list("calc.ISE", "cronbach", "f.tau", "boot.bw", 
                         "phi.K", "cfnorm", "phi.tau", "deconvolve",
                         "lnorm.var", "split.var.weight"))
  
  sim.norm <- parApply(cl, par.grid, 1, function(x)
    replicate(
      no_sim,
      calc.ISE(
        n = x[["n"]],
        no_var = x[["no_var"]],
        true_rel = x[["true_rel"]],
        dis = "norm",
        true_var = 4,
        mean = 1,
        sd = sqrt(4),
        essentially = FALSE,
        tau_congeneric = FALSE,
        tau_equivalent = x[["tau"]]
      )
    ))
  sim.norm <- cbind(as.vector(sim.norm), 
                    as.data.frame(lapply(par.grid, rep, each = no_sim)),
                    rep("norm", no_sim))
  names(sim.norm) <- c("MISE", "n", "no_var", "true_rel", "tau", "dis")
  

  sim.lnorm <- parApply(cl,par.grid, 1, function(x)
    replicate(
      no_sim,
      calc.ISE(
        n = x[["n"]],
        no_var = x[["no_var"]],
        true_rel = x[["true_rel"]],
        dis = "lnorm",
        true_var = lnorm.var(1, 0.5),
        meanlog = 1,
        sdlog = 0.5,
        essentially = FALSE,
        tau_congeneric = FALSE,
        tau_equivalent = x[["tau"]]
      )
    ))
  
  sim.lnorm <- cbind(as.vector(sim.lnorm), 
                    as.data.frame(lapply(par.grid, rep, each = no_sim)),
                    rep("lnorm", no_sim))
  names(sim.lnorm) <- c("MISE", "n", "no_var", "true_rel", "tau", "dis")
  
  result <- rbind(sim.norm, sim.lnorm)
  return(result)
}

# cl <- makeCluster(4)
# sim_results_1000 <- sim.MISE(
#   cl = cl,
#   no_sim = 1000,
#   seq_n = c(250, 2500, 10000),
#   seq_no_var = c(2, 4, 8),
#   seq_rel = c(0.80, 0.90, 0.95, 0.99),
#   seq_model = c(TRUE, FALSE)
# )
# write.csv(sim_results_1000, "sim_results_1000_classical.csv")
# stopCluster(cl)
