# -----------------------
library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(pertInv)

data_set = "GSM2396858_k562_tfs_7"
# "GSM2396861_k562_ccycle"
# "GSM2396858_k562_tfs_7"
# "GSM2396859_k562_tfs_13"
# "GSM2396860_k562_tfs_highmoi"
# "GSM2396856_dc_3hr"
# "GSM2396857_dc_0hr"
data_folder <- paste0('data_processed/', data_set)

load(file = file.path(data_folder, "batch_matrix.RData"))
load(file = file.path(data_folder, "count_matrix.RData"))
n_g <- ncol(count_matrix)
n_c <- nrow(count_matrix)
load(file = file.path(data_folder, "guide_matrix.RData"))
n_r <- ncol(guide_matrix)
covariates.dt <- fread(file.path(data_folder, "covariates.dt.csv"))

n_g <- 100
n_c <- 100

seletec_guides <- colnames(guide_matrix) %in% c("p_sgGABPA_9","p_INTERGENIC393453","p_sgYY1_3")
selected_cells <- which(rowSums(guide_matrix[,seletec_guides])>0)
n_c <- min(length(selected_cells), n_c)
selected_cells <- selected_cells[seq_len(n_c)]

batch_matrix <- batch_matrix[selected_cells,,drop=F]
batch_matrix <- batch_matrix[,colSums(batch_matrix)>0,drop=F]>0
guide_matrix <- guide_matrix[selected_cells,seletec_guides,drop=F]

count_matrix <- count_matrix[selected_cells,][,1:n_g]

n_g <- ncol(count_matrix)
n_c <- nrow(count_matrix)
n_r <- ncol(guide_matrix)
n <- n_g*n_c


set.seed(1)
ii_test = sample(n,ceiling(n/5))

dat <- list(n_c = n_c, n_g = n_g, Y = count_matrix, R=guide_matrix*1L, D= guide_matrix*1L, n_test=length(ii_test), ii_test=ii_test);
mm <- list()
fit_mc <- list()
fit_vb <- list()
ee <- list()
# ---------------------------
ii <- 1
mm[[ii]] <- stan_model_builder("stan_lib/1C_fit_vec_nb_disc.stan")
fit_vb[[ii]] <- vb(mm[[ii]], data = dat, adapt_engaged=FALSE, eta=0.1, tol_rel_obj=0.001)
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
ee[[ii]] <- extract(fit_mc[[ii]])






ii <- 2
fit_vb[[ii]] <- fit_vb[[1]]

inits <- index_samples(extract(fit_vb[[ii]]) , sample(dim((extract(fit_vb[[ii]]))[[1]])[1],4))
inits <- lapply(inits, function(x) {x$K <-  NULL;x})

mm[[ii]] <- mm[[1]]
fit_vb[[ii]] <- vb(mm[[1]], data = dat, adapt_engaged=FALSE, init=inits[[1]], eta=0.01, tol_rel_obj=0.001)
fit_vb[[ii]] <- vb(mm[[1]], data = dat, adapt_engaged=FALSE, init=inits[[1]], eta=0.01, tol_rel_obj=0.001, algorithm="fullrank")

fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, init=inits)

ee[[ii]] <- extract(fit_vb[[ii]])




sstan <- as.shinystan(fit_mc[[ii]])
sstan = as.shinystan(fit_mc[[ii]], pars=fit_mc[[ii]]@sim$pars_oi[!(fit_mc[[ii]]@sim$pars_oi %in% c("gRNA_effects", "K"))])

shinystan::launch_shinystan(sstan)


scalar_pars <- stringr::str_subset(names(fit_mc[[ii]]@sim$samples[[1]]), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
pairs(fit_mc[[ii]],pars=scalar_pars[1:ceiling(length(scalar_pars)/2)])
pairs(fit_mc[[ii]],pars=scalar_pars[-(1:ceiling(length(scalar_pars)/2))])
pairs(fit_mc[[ii]],pars=c("E_c[1]","sd_E","lp__"))
pairs(fit_mc[[ii]],pars=c("sd_mu_X","mu_X_g[1]","mu_X_g[2]","mu_X","lp__"))

dat <- lapply(ee[[ii]], mean)
lapply(ee[[ii]], mean)
lapply(ee[[ii]], sd)
lapply(ee[[ii]], dim)
lapply(ee[[ii]], function(x) prod(quantile(x,c(0.025,0.975)))>0)

wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))
MUC <- rowMeans(count_matrix)
MUC <- mean(MUC)/MUC
MUC <- exp(log(MUC)-mean(log(MUC)))
res.dt <- data.table(MUC=1/MUC, wMUC=(1/wMUC[,1]), bayes=colMeans(e[[ii]]$E_c))
figure("Correllation of size factors", GGally::ggpairs(res.dt),width=5,height=5)


res.dt <- melt(setnames(data.table(e[[ii]]$sd_X_g),colnames(count_matrix)),id.vars=integer(0),value.name = "dispersion")
res.dt <- cbind(melt(setnames(data.table(e[[ii]]$mu_X_g),colnames(count_matrix)),id.vars=integer(0),value.name = "mean"),res.dt)
res.dt <-
  res.dt[,.(mean_mean=mean(mean),mean_median=median(mean),mean_upper=quantile(mean,0.975),mean_lower=quantile(mean,0.025),
            dispersion_mean=mean(dispersion),dispersion_median=median(dispersion),dispersion_upper=quantile(dispersion,0.975),dispersion_lower=quantile(dispersion,0.025)),by=variable]

figure("dispersion estimates over mean estimates",
ggplot(res.dt,aes(x=mean_mean,y=dispersion_mean,
                  xmin=mean_lower,ymin=dispersion_lower,
                  xmax=mean_upper,ymax=dispersion_upper,
                  group=variable)) +
  geom_errorbar( width = 0.1,color="gray") +
  geom_errorbarh(height=0.01,color="gray") +
  geom_point() +
  xlab("log mean expression of gene") +
  ylab("dispersion") +
  scale_y_log10()
)

Y = sweep(count_matrix, MARGIN=1, wMUC, '*')
res.dt2 <- data.table(gene=colnames(count_matrix),dispersion=(matrixStats::colVars(Y)-colMeans(Y))/colMeans(Y)^2, log_mean=log(colMeans(Y)/mean(colMeans(Y))))
ggplot(res.dt2,aes(x=log_mean,y=dispersion)) +
  geom_point() +
  xlab("log mean expression of gene") +
  ylab("dispersion")

GGally::ggpairs(res.dt[res.dt2,.('moment matching'=dispersion, 'bayes'=dispersion_mean), on=c("variable"="gene")])

GGally::ggpairs(data.table(freq=colMeans(guide_matrix), p_R_r=colMeans(e[[ii]]$p_R_r)))


dim(apply(e[[ii]]$K, MARGIN=2:3, FUN=function(x) (quantile(x,c(0.025,0.5,0.975)))))

ggplot(data.table(detected=as.vector(guide_matrix), knockout= as.vector(t(colMeans(e[[ii]]$K)))), aes(x=detected,y=knockout))+geom_violin()
ggplot(data.table(detected=as.vector(guide_matrix), knockout= as.vector(t(colMeans(e[[ii]]$K)))), aes(x=knockout))+facet_wrap("detected",scales="free")+geom_density()



inits <- index_samples(e[[ii]] , sample(dim((e[[ii]])[[1]])[1],4))
inits <- lapply(inits, function(x) {x$K <-  NULL;x})


ii <- 2
mm[[ii]] <- stan_model_builder("stan_lib/1C_fit_vec_nb_disc.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, init=inits, iter=40,control=list(stepsize=0.06))
e[[ii]] <- extract(fit_mc[[ii]])


dat <- list(n_c = n_c, n_g = n_g, Y = count_matrix, R=guide_matrix*1L, D= guide_matrix*1L, n_test=length(ii_test), ii_test=ii_test);

R_0 <- dat$R
R_1 <- R_0
i <- sample(length(R_1),1)
R_1[i] <- !R_1[i]
dat$R <- R_1

inits <- index_samples(e[[ii]] , sample(dim((e[[ii]])[[1]])[1],4))

mm[[ii]] <- stan_model_builder("stan_lib/1C_fit_vec_nb_disc.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, init=inits, iter=40)
e[[ii]] <- extract(fit_mc[[ii]])





1)

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


ii <- 3
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


