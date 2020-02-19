################################################################################
############                 TU Dortmund University                 ############
############                    Bachelor Thesis                     ############
############    Deconvolution Problems in a Psychometric Context    ############
############                     Robin Grugel                       ############
################################################################################

source("b_thesis_estimation.R")
library(parallel)

sim.eps <- function(cl, no_sim, seq_n, seq_no_var, seq_rel, seq_model) {
  
  par.grid <- expand.grid(list(n = seq_n, no_var = seq_no_var, 
                               true_rel = seq_rel, tau = seq_model)) 
  
  
  calc.eps <- function(n, no_var, true_rel, tau_equivalent) {
  tau <- rnorm(n, mean = 1, sd = 2)
  true_var <- 4
  parts <- 1 / rep(no_var, no_var)
  
  sum_sigma_eps <- true_var * (no_var ^ 2) * (((1 - true_rel) / true_rel))
  
  if(tau_equivalent) {
    parts <- split.var.weight(no_var)
  }
  
  sigma_eps <- parts * sum_sigma_eps
  add.eps <- sapply(sigma_eps, function(x) rnorm(n, mean = 0, sd = sqrt(x)))
  x <- replicate(no_var, tau) + add.eps
  
  a <- cronbach(x)
  x.mean <- rowMeans(x)
  est_sig <- sqrt(var(x.mean) * (1 - a))
  
  est_sig
  }

  sim.eps.norm <- apply(par.grid, 1, function(x)
    replicate(
      no_sim,
      calc.eps(
        n = x[["n"]],
        no_var = x[["no_var"]],
        true_rel = x[["true_rel"]],
        tau_equivalent = x[["tau"]]
      )
    ))
  sim.eps.norm <- cbind(as.vector(sim.eps.norm), 
                    as.data.frame(lapply(par.grid, rep, each = no_sim)),
                    rep("norm", no_sim))
  names(sim.eps.norm) <- c("var_eps", "n", "no_var", "true_rel", "tau", "dis")
  
  
  result <- sim.eps.norm
  return(result)
}

# cl <- makeCluster(4)
# sim_eps_1000 <- sim.eps(
#   cl = cl,
#   no_sim = 1000,
#   seq_n = c(250, 2500, 10000),
#   seq_no_var = c(2, 4, 8),
#   seq_rel = c(0.80, 0.90, 0.95, 0.99),
#   seq_model = c(TRUE, FALSE)
# )
# write.csv(sim_eps_1000, "sim_eps_1000.csv")
# stopCluster(cl)

# plot.sim.var <- function(appendix) {
#   dat <- read.csv("sim_eps_1000.csv")
#   dat$X <- NULL
#   dat$n <- as.factor(dat$n)
#   dat$no_var <- as.factor(dat$no_var)
#   dat$tau <- as.factor(dat$tau)
#   levels(dat$n) <- c("$n = 250$", "$n = 2500$", "$n = 10000$")
#   levels(dat$dis) <- c("$\\text{log}N(1,0.25)$", "$N(1,4)$")
#   levels(dat$tau) <- c("parallel", "$\\tau$-equivalent")
#   
#   dat1 <- dat
#   dat <- aggregate(var_eps ~ ., FUN = mean, data = dat1)
#   dat2 <- aggregate(var_eps ~ ., FUN = sd, data = dat1)
#   dat$sd <- dat2$var_eps
#   
# 
#   g <- ggplot(dat, aes(x = no_var, y = var_eps, col = tau)) 
#   g <- g + geom_point(size = 0.95) + theme_bw()
#   if(appendix) {
#   g <- g +  facet_grid(as.factor(true_rel) ~ n, scales = "free")
#   } else {
#     g <- g +  facet_grid(. ~ n, scales = "free")
#   }
#   g <- g + geom_errorbar(aes(ymin=var_eps-sd, ymax=var_eps+sd), 
#                          width=.3, position=position_dodge(0.09))
#   g <- g + xlab("m")
#   g <- g + ylab("$\\widehat{\\sigma^2_{\\overline{\\epsilon}}}$")
#   g <- g + scale_color_manual(name = "CTT model", values=brewer.pal(3, "Dark2"))
# 
#   g
# }
# plot.sim.var(TRUE)

plot.sim.var2 <- function() {
  dat <- read.csv("sim_eps_1000.csv")
  dat$X <- NULL
  dat$n <- as.factor(dat$n)
  dat$no_var <- as.factor(dat$no_var)
  dat$true_rel <- as.factor(dat$true_rel)
  dat$tau <- as.factor(dat$tau)
  levels(dat$n) <- c("$n = 250$", "$n = 2500$", "$n = 10000$")
  levels(dat$dis) <- c("$\\text{log}N(1,0.25)$", "$N(1,4)$")
  levels(dat$tau) <- c("parallel", "$\\tau$-equivalent")
  
  dat1 <- dat
  dat <- aggregate(var_eps ~ ., FUN = mean, data = dat1)
  dat2 <- aggregate(var_eps ~ ., FUN = sd, data = dat1)
  dat$sd <- dat2$var_eps
  
  g <- ggplot(dat, aes(x = no_var, y = sd, col = tau, group = tau)) 
  g <- g + geom_point(size = 0.35) + geom_line(size = 0.55) + theme_bw()
  g <- g + facet_grid(n ~ true_rel)
  g <- g + xlab("m")
  g <- g + ylab("$\\widehat{\\sigma}_{\\widehat{\\sigma^2_{\\overline{\\epsilon}}}}$")
  g <- g + scale_color_manual(name = "CTT model", 
                              values=brewer.pal(3, "Dark2"))

  g
}
plot.sim.var2()

get.var.table <- function(appendix) {
  dat <- read.csv("sim_eps_1000.csv")
  dat$X <- NULL
  dat$n <- as.factor(dat$n)
  dat$no_var <- as.factor(dat$no_var)
  dat$true_rel <- as.factor(dat$true_rel)
  dat$tau <- as.factor(dat$tau)
  levels(dat$n) <- c("$n = 250$", "$n = 2500$", "$n = 10000$")
  levels(dat$tau) <- c("parallel", "$\\tau$-equivalent")
  
  dat1 <- dat
  dat <- aggregate(var_eps ~ ., FUN = mean, data = dat1)
  dat2 <- aggregate(var_eps ~ ., FUN = sd, data = dat1)
  dat$sd <- dat2$var_eps
  
  dat$sd <- round(dat$sd, digits = 5)
  
  jo <-
    tabular(
      Heading() * n * Heading() * no_var * Heading() * 
        true_rel ~ Heading() *
        tau * Heading() * identity * (sd),
      data = dat
    )
  latex(jo)
}

get.table <- function(logN) {
  dat.new <- read.csv("sim_results_1000_deconvolvefinal.csv")
  dat.new$fun <- "\\widehat{f}_{\\tau}(\\cdot)"
  dat.kde <- read.csv("sim_results_1000_classicalfinal.csv")
  dat.kde$fun <- "$\\widehat{f}_{X}(\\cdot)$"
  dat <- rbind(dat.new, dat.kde)
  dat$X <- NULL
  dat$n <- as.factor(dat$n)
  dat$no_var <- as.factor(dat$no_var)
  dat$fun <- as.factor(dat$fun)
  dat$tau <- as.factor(dat$tau)
  dat$true_rel <- as.factor(dat$true_rel)
  levels(dat$n) <- c("$n = 250$", "$n = 2500$", "$n = 10000$")
  levels(dat$dis) <- c("$\\text{log}N(1,0.5)$", "$N(1,4)$")
  levels(dat$no_var) <- c("$m = 2$", "$m = 4$", "$m = 8$")
  levels(dat$tau) <- c("parallel", "$\\tau$-equivalent")
  
  if(logN) {
    dat <- subset(dat, dis == "$\\text{log}N(1,0.5)$")
  } else {
    dat <- subset(dat, dis == "$N(1,4)$")
  }
  dat$MISE <- round(dat$MISE, digits = 5)
  
  jo <-
    tabular(
      Heading() * n * Heading() * no_var * Heading() * true_rel ~ Heading() *
        tau * Heading() * fun * Heading() * identity * (MISE),
      data = dat
    )
  latex(jo)
}


plot.sim <- function(logN) {
  dat.new <- read.csv("sim_results_1000_deconvolvefinal.csv")
  dat.new$fun <- "$\\widehat{f}_{\\tau}(\\cdot)$"
  dat.kde <- read.csv("sim_results_1000_classicalfinal.csv")
  dat.kde$fun <- "$\\widehat{f}_{\\overline{X}}(\\cdot)$"
  dat <- rbind(dat.new, dat.kde)
  dat$X <- NULL
  dat$n <- as.factor(dat$n)
  dat$no_var <- as.factor(dat$no_var)
  dat$fun <- as.factor(dat$fun)
  dat$tau <- as.factor(dat$tau)
  levels(dat$n) <- c("$n = 250$", "$n = 2500$", "$n = 10000$")
  levels(dat$dis) <- c("$\\text{log}N(1,0.25)$", "$N(1,4)$")
  levels(dat$no_var) <- c("$m = 2$", "$m = 4$", "$m = 8$")
  levels(dat$tau) <- c("parallel", "$\\tau$-equivalent")
  
  dat1 <- dat
  dat <- aggregate(MISE ~ ., FUN = mean, data = dat1)
  dat2 <- aggregate(MISE ~ ., FUN = sd, data = dat1)
  dat$sd <- dat2$MISE
  
  dat <- subset(dat, dis == "$N(1,4)$" & fun == "$\\widehat{f}_{\\tau}(\\cdot)$")
  
  
  ggplot(dat, aes(x = true_rel, y = MISE, col = tau, linetype = fun)) +
    geom_point(size = 0.55) + geom_line(size = 0.35) + theme_bw() +
    facet_grid(no_var ~ n) + 
    geom_errorbar(aes(ymin=MISE-sd, ymax=MISE+sd), width=.01,
                  position=position_dodge(0.001)) + 
    xlab("$\\text{Rel}(\\overline{X})$") +
    ylab("$\\widehat{\\text{MISE}}\\{\\widehat{g}(\\cdot)\\}$") +
    scale_color_manual(name = "CTT model", values=brewer.pal(3, "Dark2")) +
    scale_linetype_discrete(name = "Estimator $\\widehat{g}(\\cdot)$")
}

plot.years <- function() {
  bl <- read.csv2("FE_1933_2006_data.csv")
  bl2 <- read.csv2("FE_1933_2006_variables_and_coding.csv")
  names(bl) <- bl2$Variables
  bl[bl == -9] <- NA
  neu <- subset(bl, select = c("group", "IDNR", "sex", "grade", "age"))
  neu$a1 <- rowSums(bl[, 6:21], na.rm = TRUE) * 2
  neu$a2 <- rowSums(bl[, 22:41], na.rm = TRUE) * 2
  neu$a3 <- rowSums(bl[, 42:65], na.rm = TRUE)
  neu$a4 <- rowSums(bl[, 66:105], na.rm = TRUE) - rowSums(!bl[, 66:105], 
                                                          na.rm = TRUE)
  neu$a5 <- bl[, 106] * 0.3
  neu$b1 <- rowSums(bl[, 109:130], na.rm = TRUE) * 2
  neu$b2 <- rowSums(bl[, 131:170], na.rm = TRUE)
  neu$b3 <- rowSums(bl[, 171:210], na.rm = TRUE) - rowSums(!bl[, 171:210], 
                                                           na.rm = TRUE)
  neu$b4 <- rowSums(bl[, 211:242], na.rm = TRUE)
  neu$b5 <- rowSums(bl[, 243:292], na.rm = TRUE) - rowSums(!bl[, 243:292], 
                                                           na.rm = TRUE)
  X <- subset(neu, group == 2,
              select = c("a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4", "b5"))
  a <- cronbach(X)
  X.bar <- rowMeans(X)
  siggi <- sqrt(var(X.bar) * (1 - a))
  out1 <- as.data.frame(f.tau(X.bar, siggi))
  X <- subset(neu, group == 1,
              select = c("a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4", "b5"))
  a <- cronbach(X)
  X.bar <- rowMeans(X)
  siggi <- sqrt(var(X.bar) * (1 - a))
  out2 <- as.data.frame(f.tau(X.bar, siggi))
  out3 <- rbind(out1, out2)
  out3$year <-
    as.factor(rep(c("2006", "1933/36"), each = length(out1$w)))
  ggplot(out3, aes(x = w, y = Re(eval), col = year)) +
    geom_line(size = 1) + theme_bw() + xlab("$\\omega$") +
    ylab("$\\widehat{f}_{\\tau_{\\text{time}}}(\\omega)$") +
    scale_color_manual(name = "Year", 
                       values = brewer.pal(3, "Dark2"))
}

plot.comparison <- function(n, no_var, seq_rel, dis, true_var, essentially,
                            tau_congeneric, tau_equivalent,
                            selected.seed, ...) {
  res <- list()
  
  parts <- 1 / rep(no_var, no_var)
  lambda <- rep(1, no_var)
  translation <- rep(0, no_var)
  
  if (tau_equivalent) {
    parts <- split.var.weight(no_var)
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
      c(0.27,-0.27), # with additional constraint
      # c(0.87, -0.17),
      c(0.31, -0.09, -0.12, -0.1),
      # c(1.1, 1.09, 0.12, -0.1), # with additional constraint
      c(0.1, -0.68, 0.57, -0.2, 0.05, 0.3, 0.13, -0.27)
      # c(0.1, 6.8, -0.57, -0.2, 0.05, 2.3, 0.13, -0.27) # with additional constraint
    )
  }
  
  tau <- do.call(paste0("r", dis), list(n = n, ...))
  dat <- data.frame()
  
  for (r in seq_rel) {
    set.seed(selected.seed)
    sum_sigma_eps <-
      true_var * (mean(lambda) ^ 2) * (no_var ^ 2) * (((1 - r) / r))
    sigma_eps <- parts * sum_sigma_eps
    multiply.load <-
      matrix(lambda,
             nrow = n,
             ncol = no_var,
             byrow = TRUE)
    add.trans <-
      matrix(translation,
             nrow = n,
             ncol = no_var,
             byrow = TRUE)
    add.eps <-
      sapply(sigma_eps, function(x)
        rnorm(n, mean = 0, sd = sqrt(x)))
    x <-
      ((replicate(no_var, tau)  + add.trans) * multiply.load) + add.eps
    
    a <- cronbach(x)
    x.bar <- rowMeans(x)
    siggi <- sqrt(var(x.bar) * (1 - a))
    out <- f.tau(x.bar, siggi)
    dat <- rbind(dat, out)
  }
  
  dat$reliability <- as.factor(rep(
    paste0(
      "$\\alpha$ = ",
      seq_rel
    ),
    each = length(out$w)
  ))
  
  tit <- ifelse(tau_equivalent, "$\\tau$-equivalent", "parallel")
  
  res$g <- ggplot(dat, aes(x = w, y = Re(eval))) +
    facet_wrap(. ~ reliability) + geom_line() +
    geom_line(aes(x = w, y = do.call(paste0("d", dis), list(x = w, ...)))
              , col = "red") +
    xlab("$\\omega$") + ylab("$\\widehat{f}_{\\tau}(\\omega)$") + theme_bw() + ggtitle(tit)
  
  return(res$g)
}

