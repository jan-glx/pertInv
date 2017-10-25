functions{
#include index_set_diff.stan
}

data {
int<lower=0> n_c; // number of cells
int<lower=0> n_g; // number of genes
int<lower=0> Y[n_g*n_c]; // counts
int<lower=0> n_test; // number of holdout
int<lower=1,upper=n_c*n_g> ii_test[n_test]; // indices of holdout
}

transformed data {
int<lower=0> n = n_c*n_g;
int<lower=0> n_train = n-n_test; // number of not holdout
int<lower=1,upper=n> ii_train[n_train]; // indices of holdout
ii_train = index_set_diff(ii_test, n_test, n);
}

parameters {
matrix[n_g,n_c] X; // log expression
row_vector[n_c] E; // log exposure 'size factor' or 'library size'
vector<lower=0>[n_g] var_X_g;
vector[n_g] mu_X_g;
real<lower=0> var_mu_X;
real<lower=0> sd_E;
//real mode_var_gene;
//real<lower=0> scale_vare_gene;
real mu_X;
}

transformed parameters {
}

model {
//target += normal_lpdf(mode_var_gene | 0, 3);
//target += normal_lpdf(scale_vare_gene | 1, 10);
//target += normal_lpdf(var_gene | mode_var_gene, scale_vare_gene); //cauchy_lpdf
target += normal_lpdf(to_vector(X) | to_vector(rep_matrix(mu_X+mu_X_g, n_c)), to_vector(rep_matrix(sqrt(var_X_g), n_c)));

//target += inv_gamma_lpdf(sd_E | 40, 10);
//target += normal_lpdf(E | 0, sd_E);

//target += inv_gamma_lpdf(var_gene_ln_expression | 5, 5);
target += normal_lpdf(mu_X_g | 0, sqrt(var_mu_X));

//target += normal_lpdf(mu_X | 0, 10);

target += poisson_lpmf(Y[ii_train]  | to_vector(exp(X+rep_matrix(E, n_g)))[ii_train]);
}

generated quantities{
/*int sim_counts[n_g, n_c];
for (g in 1:n_g) {
  for (c in 1:n_c) {
    sim_counts[g,c] = poisson_rng(expression[g,c]);
  }
}*/
real ll_test;
ll_test = poisson_lpmf(Y[ii_test]  | to_vector(exp(X+rep_matrix(E, n_g)))[ii_test]);
}
