functions{
#include index_set_diff.stan
#include common.stan
}

data {
int<lower=0> n_c; // number of cells
int<lower=0> n_g; // number of genes
int<lower=0> n_r; // number of gRNAs
int<lower=0> Y[n_g*n_c]; // counts
matrix<lower=0,upper=1>[n_r,n_c] D; // gRNA detected
int<lower=0> n_test; // number of holdout
int<lower=1,upper=n_c*n_g> ii_test[n_test]; // indices of holdout

}

transformed data {
int<lower=0> n = n_c*n_g;
int<lower=0> n_train = n-n_test; // number of not holdout
int<lower=1,upper=n> ii_train[n_train]  = index_set_diff(ii_test, n_test, n); // indices of holdout
int<lower=1,upper=2> D_int[n_c,n_r]; // gRNA detected 1=True 2=False [n_r,n_c]
for (c in 1:n_c) {
  for (r in 1:n_r) {
    D_int[c,r] = (D[r,c]>0.5)? 1 : 2;
  }
}
}

parameters {
matrix[n_g,n_c] X; // log expression
row_vector[n_c] E_c; // log exposure 'size factor' or 'library size'
vector<lower=0>[n_g] sd_X_g; // std. dev.s of genes' expression
real<lower=0> alpha_var_X_g;
real<lower=0> beta_var_X_g;
vector[n_g] mu_X_g; // means of genes' expression
real<lower=0> sd_mu_X;
real<lower=0> sd_E;
//real mode_var_gene;
//real<lower=0> scale_vare_gene;
real mu_X;
matrix[n_g,n_r] gRNA_effects; // total effects of gRNAs
real<lower=0> sd_gRNA_effects;
vector<lower=0,upper=1>[2] p_D_given_R;
matrix[n_r,n_c] K; // gRNA detected
row_vector<lower=0,upper=1>[n_r] p_R_r; // gRNA_prior_probs
real<lower=0> alpha_p_R_r;
real<lower=0> beta_p_R_r;
}

transformed parameters {
matrix<lower=0,upper=1>[2,n_r] p_R_r_given_D_r; // [[R1_given_D_1, R2_given_D_1, ...],
                                                //  [R1_given_D_0, R2_given_D_0, ...]]

p_R_r_given_D_r = inv(1 - p_D_given_R * (1 - inv(p_R_r)));
}

model {
//target += normal_lpdf(mode_var_gene | 0, 3);
//target += normal_lpdf(scale_vare_gene | 1, 10);
target += inv_gamma_lpdf(sd_X_g | alpha_var_X_g, beta_var_X_g); //cauchy_lpdf
target += normal_lpdf(to_vector(X) | to_vector(rep_matrix(mu_X+mu_X_g, n_c) + gRNA_effects*K), to_vector(rep_matrix(sd_X_g, n_c)));

target += normal_lpdf(to_vector(gRNA_effects) | 0, sd_gRNA_effects);

for (c in 1:n_c) {
  for (r in 1:n_r) {
    target += log_mix(p_R_r_given_D_r[D_int[c,r],r], normal_lpdf(K[r,c] | 1, 0.1), normal_lpdf(K[r,c] | 0, 0.1));
  }
}

target += beta_lpdf(p_R_r | alpha_p_R_r, beta_p_R_r);


//target += inv_gamma_lpdf(sd_E | 40, 10);
target += normal_lpdf(E_c | 0, sd_E);

//target += inv_gamma_lpdf(var_gene_ln_expression | 5, 5);
target += normal_lpdf(mu_X_g | 0, sd_mu_X);

//target += normal_lpdf(mu_X | 0, 10);

target += poisson_lpmf(Y[ii_train]  | to_vector(exp(X+rep_matrix(E_c, n_g)))[ii_train]);
}

generated quantities{
/*int sim_counts[n_g, n_c];
for (g in 1:n_g) {
  for (c in 1:n_c) {
    sim_counts[g,c] = poisson_rng(expression[g,c]);
  }
}*/
real ll_test;
ll_test = poisson_lpmf(Y[ii_test]  | to_vector(exp(X+rep_matrix(E_c, n_g)))[ii_test]);
}
