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
cat(fit_mc[[ii]]@stanmodel@model_code, file="stan_lib/to_read.stan")
ee[[ii]] <- extract(fit_mc[[ii]])
list2env(dat,globalenv())

# check predictions
for (i in names(ee[[ii]])) {
  tmp = ee[[ii]][[i]]
  dim(tmp) = c(dim(tmp)[1],prod(dim(tmp)[-1]))
  hist(tmp[,1], main=i)
  abline(v=dat[[i]][1],col="red")
}


ii <- 3
mm[[ii]] <-  stan_model_builder("stan_lib/1_fit_vec.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
ee[[ii]] <- extract(fit_mc[[ii]])



ii <- 4
mm[[ii]] <-  mm[[ii-1]]
inits = index_samples(ee[[1]], is=2:5)
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, init = inits)
ee[[ii]] <- extract(fit_mc[[ii]])

dt <- rbindlist(lapply(fit_mc[[ii]]@sim$samples,function(x) as.data.table(c(x,attr(x, "sampler_params")))[,chain__:=stringi::stri_rand_strings(1, 3)]))
dt[,i__:=seq_len(.N),by=chain__]

scalar_pars <- stringr::str_subset(colnames(dt), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
dt <- dt[,scalar_pars, with=FALSE]
ggplot(dt[i__>1000],aes(y=energy__,x=mu_log_p_R_r))+geom_point()
ggplot(melt(dt[i__>1000],id.vars=c("energy__","chain__","i__")),aes(y=energy__,x=value))+geom_point()+facet_wrap("variable", scales="free_x")+
 geom_smooth(method="lm")

scalar_pars <- stringr::str_subset(names(fit_mc[[ii]]@sim$samples[[1]]), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
pairs(fit_mc[[ii]],pars=scalar_pars[1:ceiling(length(scalar_pars)/2)])
pairs(fit_mc[[ii]],pars=scalar_pars[-(1:ceiling(length(scalar_pars)/2))])
pairs(fit_mc[[ii]],pars=c("E_c[1]","sd_E","lp__"))
pairs(fit_mc[[ii]],pars=c("sd_mu_X","mu_X_g[1]","mu_X_g[2]","mu_X","lp__"))

sstan = as.shinystan(fit_mc[[ii]], pars=fit_mc[[ii]]@sim$pars_oi[!(fit_mc[[ii]]@sim$pars_oi %in% c("X", "X_train", "X_test"))])

shinystan::launch_shinystan(sstan)






dat <- index_sample(e[[ii]], 1)
dat$n_c  <-  ncol(dat$X) # number of cells
dat$n_g  <-  nrow(dat$X) # number of genes
dat$n_r  <-  ncol(dat$D) # number of sgRNAs
n <- dat$n_c * dat$n_g
dat$ii_test <-  sample(n, ceiling(n/5))
dat$n_test <- length(dat$ii_test)
