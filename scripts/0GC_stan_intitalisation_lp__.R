
index_first_dim <- function(x, i) {
  sx <- quote(x[i])
  for (j in seq_len(length(dim(x))-1)) {
    sx[[j+3]] <- quote(x[,])[[3]] #add empty symbol
  }
  eval(sx)
}

index_sample <-  function(samples, i) {
  lapply(samples, index_first_dim, i=i)
}

index_samples <-  function(samples, is) {
  lapply(is, function(i) index_sample(samples, i))
}

# ---------------------------
library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix[1:100,1:100]#count_matrix#count_matrix[1:1000,1:100]
guide_matrix <- guide_matrix[1:100,]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)


library(rstan)
options(mc.cores = parallel::detectCores())

library(shinystan)
thecounts <- as.vector(t(count_matrix))

dat <- list(N_cells = nrow(count_matrix), N_genes = ncol(count_matrix),  counts = thecounts);

# ---------------------------
ii <- 1
mm=list()
fit_vb=list()
mm[[ii]] <-  stan_model("scripts/model_poisonly.stan")
# ======

fit_vb[[ii]] <- vb(mm[[ii]], data = dat, sample_file = 'norm2.csv', tol_rel_obj=0.001)#, iter=100000, adapt_engaged=FALSE, eta=0.0001)
# ---------------------------------------
e <- extract(fit2)
inits = index_samples(e, seq_len(4)) #[!(names(e)%in% c("sim_counts", "lp__"))]

fit3 <- sampling(m2, data = dat, init=inits)

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
