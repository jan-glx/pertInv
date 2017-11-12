# -----------------------
library(pertInv)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")
load(file="results/batch.RData")

count_matrix <- count_matrix[1:100,1:100]#count_matrix#count_matrix[1:1000,1:100]
guide_matrix <- guide_matrix[1:100,]#count_matrix#count_matrix[1:1000,1:100]
batch <- model.matrix(~ batch + 0, as.data.table(batch[1:100]))#count_matrix#count_matrix[1:1000,1:100]
n_g <- ncol(count_matrix)
n_c <- nrow(count_matrix)
n_r <- ncol(guide_matrix)
n <- n_g*n_c


thecounts <- as.vector(t(count_matrix))
set.seed(1)
ii_test = sample(n,ceiling(n/5))

dat <- list(n_c = n_c, n_g = n_g, Y = thecounts, n_test=length(ii_test), ii_test=ii_test);
mm <- list()
fit_mc <- list()
fit_vb <- list()
e <- list()
# ---------------------------
# ii <- 1
# mm[[ii]] <- readRDS("stan_lib/model_simplest_hirachical_A.rds") #stan_model_builder("stan_lib/model_simplest_hirachical_B.stan")
# fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
# e[[ii]] <- extract(fit_mc[[ii]])

# dat <- lapply(e[[ii]], mean)


dat <-  list(
  n_c = 64, # number of cells
  n_g = 16, # number of genes
  n_r =  4, # number of gRNAs
  # hyper parameters
  sd_mu_X = 1,
  sd_E = 0.2,
  mu_log_var_X_g = log(1),
  sd_gRNA_effects = 0.3,
  mu_X = 1,
  p_D_given_R = c(0.9, 0.01)
)
dat$sd_log_var_X_g <- log(1.3)
dat$mu_log_p_R_r <- log(1/dat$n_r)
dat$sd_log_p_R_r <- log(1.3)


ii <- 2
mm[[ii]] <-  stan_model_builder("stan_lib/sim_params_with_guide.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, warmup=1000, chains=4)
e[[ii]] <- extract(fit_mc[[ii]])

dat  <-  c(dat, index_sample(e[[ii]], i=1))

ii <- 3
mm[[ii]] <-  stan_model_builder("stan_lib/sim_data_with_guide.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])

dat <- index_sample(e[[ii]], 1)
dat$n_c  <-  ncol(dat$X) # number of cells
dat$n_g  <-  nrow(dat$X) # number of genes
dat$n_r  <-  ncol(dat$D) # number of sgRNAs
n <- dat$n_c * dat$n_g
dat$ii_test <-  sample(n, ceiling(n/5))
dat$n_test <- length(dat$ii_test)

ii <- 4
mm[[ii]] <-  stan_model_builder("stan_lib/model_with_guide_A.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, iter = 4000, control = list(adapt_delta = 0.999, stepsize=0.01)) # does not fit correctly
e[[ii]] <- extract(fit_mc[[ii]])

scalar_pars <- stringr::str_subset(names(fit_mc[[ii]]@sim$samples[[1]]), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
pairs(fit_mc[[ii]],pars=scalar_pars[1:ceiling(length(scalar_pars)/2)])
pairs(fit_mc[[ii]],pars=scalar_pars[-(1:ceiling(length(scalar_pars)/2))])
pairs(fit_mc[[ii]],pars=c("E_c[1]","sd_E","lp__"))
pairs(fit_mc[[ii]],pars=c("sd_mu_X","mu_X_g[1]","mu_X_g[2]","mu_X","lp__"))

sstan = as.shinystan(fit_mc[[ii]], pars=fit_mc[[ii]]@sim$pars_oi[!(fit_mc[[ii]]@sim$pars_oi %in% c("X", "X_train", "X_test"))])

shinystan::launch_shinystan(sstan)


plot(get_sampler_params(fit_mc[[ii]], inc_warmup = FALSE)[[4]][,"energy__"])



ii <- 6
mm[[ii]] <-  stan_model_builder("stan_lib/model_with_guide_C.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test) # -466.0178  -445.5888


dat$n_test = ceiling(dat$n_c/5)
dat$ii_test = sample(dat$n_c, dat$n_test)

ii <- 7
mm[[ii]] <-  stan_model_builder("stan_lib/model_with_guide_D.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])



inf.dt <-  data.table(infered=as.vector(colMeans(e[[ii]]$D_mm)[,dat$ii_test]),
                      label=as.vector(t(dat$D[dat$ii_test,])),
                      true=as.vector(dat$K[,dat$ii_test]))

ggplot(inf.dt,aes(x=true==1,y=infered))+geom_jitter()+stat_summary(fun.data=mean_se,size=2,color="red")+xlab("sgRNA present (R)")

mu_X_g <- lapply(e[[2]], drop)$mu_X_g
p_R_r <- lapply(e[[2]], drop)$p_R_r
sd_X_g <- lapply(e[[2]], drop)$sd_X_g


sstan = as.shinystan(fit_mc[[ii]], pars=fit_mc[[ii]]@sim$pars_oi[!(fit_mc[[ii]]@sim$pars_oi %in% c("X", "X_train", "X_test"))])

shinystan::launch_shinystan(sstan)




K <-  dat$K
D <-  dat$D

jitter(K, index_sample(e[[ii]],1000)$K)

scalar_pars <- stringr::str_subset(names(fit_mc[[ii]]@sim$samples[[1]]), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
sstan <- as.shinystan(fit_mc[[ii]], pars=scalar_pars)
shinystan::launch_shinystan(fit_mc[[ii]])

mean(e[[ii]]$ll_test)
t.test(e[[ii]]$ll_test,e[[ii-1]]$ll_test)


plot((colMeans(e[[ii]]$E_c)), log(rowSums(count_matrix)))
plot((colMeans(e[[ii]]$mu_X_g)), log(colSums(count_matrix)))


