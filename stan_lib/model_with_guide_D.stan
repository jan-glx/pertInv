functions{
#include index_set_diff.stan
}

data {
#include _data_size_with_guide_A.stan
#include _dicrete_data_with_guide_A.stan

int<lower=0> n_test; // number of holdout
int<lower=1,upper=n_c> ii_test[n_test]; // indices of holdout
}

transformed data {
int<lower=0> n = n_c;
int<lower=0> n_train = n-n_test; // number of not holdout
int<lower=1,upper=n> ii_train[n_train] = index_set_diff(ii_test, n_test, n); // indices of holdout
matrix<lower=0,upper=1>[n_r,n_c] D_m = to_matrix(to_array_1d(D), n_r, n_c); // sgRNA detected
int<lower=0> Y_flat[n_g*n_c] = to_array_1d(Y);
}

parameters {
#include _hira_prior_params_with_guide_A.stan
#include _cont_hira_params_with_guide_A.stan
matrix[n_g,n_c] X; // log expression
matrix<lower=0,upper=1>[n_r,n_test] D_m_latent; // sgRNA detected

// non-hirachical parameters
real mu_X;
}

transformed parameters {
matrix[n_g,n_train] X_train = X[,ii_train]; // log expression
matrix[n_g,n_test]  X_test  = X[,ii_test]; // log expression
matrix<lower=0,upper=1>[n_r,n_c] D_mm = D_m; // sgRNA detected
D_mm[,ii_test] = D_m_latent;
}

model {
#include _cont_hira_dists_with_guide_A.stan

// prior for hirachical paramters
target += normal_lpdf(mu_log_p_R_r | -log(n_r), log(n_r)/2);
target += normal_lpdf(mu_log_var_X_g | -1, 2);
target += cauchy_lpdf(sd_gRNA_effects | 0, 1);
target += cauchy_lpdf(sd_E | 0, 1);
target += cauchy_lpdf(sd_mu_X | 0, 1);
target += cauchy_lpdf(sd_log_p_R_r | 0, 1);
target += cauchy_lpdf(sd_log_var_X_g | 0, 1);


target += normal_lpdf(to_vector(X) | to_vector(rep_matrix(mu_X+mu_X_g, n_c) + gRNA_effects*D_mm), to_vector(rep_matrix(sd_X_g, n_c)));
target += poisson_lpmf(Y_flat  | to_vector(exp(X+rep_matrix(E_c, n_g))));
}

generated quantities{
real ll_test;
ll_test = poisson_lpmf(Y_flat[ii_test]  | to_vector(exp(X+rep_matrix(E_c, n_g)))[ii_test]); // conditional predictive ordinate
}
