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
ii <- 1
mm[[ii]] <-  stan_model_builder("stan_lib/model_simplest.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test)

fit_vb[[ii]] <- vb(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_vb[[ii]])
mean(e[[ii]]$ll_test)

ii <- 2
mm[[ii]] <-  stan_model_builder("stan_lib/model_simplest_hirachical.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test)

t.test(e[[1]]$ll_test, e[[2]]$ll_test)
mean(e[[1]]$ll_test) #-3014.313
mean(e[[2]]$ll_test) #-2987.932 # make per gene mean hirachical



ii <- 3
mm[[ii]] <-  stan_model_builder("stan_lib/model_simplest_hirachical_A.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test) # -2834.916 # making additionally per gene variance hirachical

ii <- 4
mm[[ii]] <-  stan_model_builder("stan_lib/model_simplest_hirachical_B.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test) # -2827.271 # making additionally size factor hirachical
t.test(e[[ii]]$ll_test,e[[ii-1]]$ll_test)


dat <- list(n_c = n_c, n_g = n_g, Y = thecounts, n_test=length(ii_test), ii_test=ii_test, n_r = n_r,
            D = 1*t(guide_matrix));

ii <- 5
mm[[ii]] <-  stan_model_builder("stan_lib/model_with_guide_A.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test) # -2827.579 # with guide info
t.test(e[[ii]]$ll_test,e[[ii-1]]$ll_test) # not significant


plot((colMeans(e[[ii]]$E_c)), log(rowSums(count_matrix)))
plot((colMeans(e[[ii]]$mu_X_g)), log(colSums(count_matrix)))


dat <- list(n_c = n_c, n_g = n_g, Y = thecounts, n_test=length(ii_test), ii_test=ii_test, n_r = n_r,
            D = 1*t(guide_matrix));
ii <- 6
mm[[ii]] <-  stan_model_builder("stan_lib/model_with_guide_B.stan")
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
e[[ii]] <- extract(fit_mc[[ii]])
mean(e[[ii]]$ll_test) # -2781.355 # with guide info
t.test(e[[ii]]$ll_test,e[[ii-2]]$ll_test) # significant


