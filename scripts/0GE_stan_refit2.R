library(pertInv)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(1)


mm <- list()
fit_mc <- list()
fit_vb <- list()
ee <- list()
# ---------------------------


dat <- list(n_c = 64, n_g = 32, n_r = 16);

ii <- 1
mm[[ii]] <-  stan_model_builder("stan_lib/1_simulate.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, chains=1, iter = 1000, warmup=0, algorithm="Fixed_param")
ee[[ii]] <- extract(fit_mc[[ii]])

dat  <-  c(dat, index_sample(ee[[ii]], i=1))

ii <- 2
mm[[ii]] <-  stan_model_builder("stan_lib/1_fit.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
ee[[ii]] <- extract(fit_mc[[ii]])

dat <- index_sample(e[[ii]], 1)
dat$n_c  <-  ncol(dat$X) # number of cells
dat$n_g  <-  nrow(dat$X) # number of genes
dat$n_r  <-  ncol(dat$D) # number of sgRNAs
n <- dat$n_c * dat$n_g
dat$ii_test <-  sample(n, ceiling(n/5))
dat$n_test <- length(dat$ii_test)
