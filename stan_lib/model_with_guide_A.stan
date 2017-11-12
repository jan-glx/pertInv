functions{
#include index_set_diff.stan
}

data {
#include _data_size_with_guide_A.stan
#include _dicrete_data_with_guide_A.stan

int<lower=0> n_test; // number of holdout
int<lower=1,upper=n_c*n_g> ii_test[n_test]; // indices of holdout
}

transformed data {
int<lower=0> n = n_c*n_g;
int<lower=0> n_train = n-n_test; // number of not holdout
int<lower=1,upper=n> ii_train[n_train] = index_set_diff(ii_test, n_test, n); // indices of holdout
matrix<lower=0,upper=1>[n_r,n_c] D_m = to_matrix(to_array_1d(D), n_r, n_c); // sgRNA detected
int<lower=0> Y_flat[n_g*n_c] = to_array_1d(Y);
}

parameters {
#include _hira_prior_params_with_guide_A.stan
#include _cont_hira_params_with_guide_A.stan
matrix[n_g,n_c] X; // log expression

// non-hirachical parameters
real mu_X;
}

transformed parameters {
}

model {
#include _hira_prior_params_dists_with_guide_A.stan
#include _cont_hira_dists_with_guide_A.stan
target += normal_lpdf(mu_X | 0, 1);
target += normal_lpdf(to_vector(X) | to_vector(rep_matrix(mu_X+mu_X_g, n_c) + gRNA_effects*D_m), to_vector(rep_matrix(sd_X_g, n_c)));
target += poisson_lpmf(Y_flat[ii_train]  | to_vector(exp(X+rep_matrix(E_c, n_g)))[ii_train]);
}

generated quantities{
real ll_test;
ll_test = poisson_lpmf(Y_flat[ii_test]  | to_vector(exp(X+rep_matrix(E_c, n_g)))[ii_test]);
}
