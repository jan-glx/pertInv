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

