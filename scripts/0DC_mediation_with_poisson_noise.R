n <- 1000
A <- rnorm(n)
B <- (rnorm(n)+A)/2
C <- (rnorm(n)+B)/2

mA <- rpois(n,lambda=exp(A))
mB <- rpois(n,lambda=exp(B))
mC <- rpois(n,lambda=exp(C))

sd_noise <- 0.3 # if you reduce this further, the second test below quickly looses power
c_measured <- C + rnorm(n, sd = sd_noise)
ppcor::pcor.test(A,B,C) # not significant
ppcor::pcor.test(A,B,c_measured) # significant
# I am looking for some thing like that:
ppcor::pcor.test(A,B,c_measured, z_noise_sd=sd_noise) # not significant


dat <- list(n=n, mA=mA, mB=mB, mC=mC)
stan_model <- stan_model_builder("stan_lib/2_mediation_poisson.stan")
#fit_vb <- vb(mm[[ii]], data = dat, adapt_engaged=FALSE, eta=0.1, tol_rel_obj=0.001)
fit <- sampling(stan_model, data = dat)
ee <- extract(fit)
summarize_sample(ee)[variable %in% c("beta_AB","beta_AC","beta_BC")]
# variable        mean      median      lower     upper
# 1:  beta_AB  0.49516446  0.49437076  0.4002955 0.5969427
# 2:  beta_AC -0.05567849 -0.04537204 -0.2745368 0.1256840
# 3:  beta_BC  0.54754365  0.53296445  0.2813668 0.9120366
