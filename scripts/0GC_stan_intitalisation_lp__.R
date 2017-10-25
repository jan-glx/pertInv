# ---------------------------
library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix[1:100,1:100]#count_matrix#count_matrix[1:1000,1:100]
guide_matrix <- guide_matrix[1:100,]#count_matrix#count_matrix[1:1000,1:100]
n_g <- ncol(count_matrix)
n_c <- nrow(count_matrix)
n <- n_g*n_c


thecounts <- as.vector(t(count_matrix))
set.seed(1)
ii_test = sample(n,ceiling(n/5))

dat <- list(n_c = n_c, n_g = n_g, Y = thecounts, n_test=length(ii_test), ii_test=ii_test);
m1 <- stan_model_builder("stan_lib/model_simplest.stan")

library(rstan)

stime <- proc.time()
fit_vb <- rstan::vb(m1, data = dat, sample_file = 'norm2.csv', eta=0.1, adapt_engaged=FALSE, output_samples=4)
inits = index_samples(extract(fit_vb), seq_len(4))
fit_mc_init <- sampling(m1, data = dat, init=inits)
time_init <- proc.time()-stime

stime <- proc.time()
fit_mc_no_init <- sampling(m1, data = dat)
time_no_init <- proc.time()-stime

time_init
time_no_init

my_sso <- as.shinystan(fit3, pars = c("gene_variace_mode","gene_variance_scale","mean_ln_mean_expression", "sd_ln_mean_expression","sd_sf"))



t_counts <- t(count_matrix)
launch_shinystan(my_sso)


e_vb <- extract(fit_vb)
e_mc <- extract(fit_mc)


tmp = e_vb[!(names(e_vb) %in% "lp__")]


fit_vb <-  fit_vb2
e_vb <- extract(fit_vb)
lp___vb <- sapply(seq_len(dim(e_vb[[1]])[1]), function(i) {log_prob(fit_vb, unconstrain_pars(fit_vb, index_sample(e_vb, i)))})

lp___mc <- e_mc$lp__

ggplot(rbind(data.table(method="NUTS", log_posterior=lp___mc),
             data.table(method="VB",   log_posterior=lp___vb)),
       aes(x=method,y=log_posterior))+geom_boxplot()
