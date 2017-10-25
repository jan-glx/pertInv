functions{
#include index_set_diff.stan
#include common.stan
}

data {
#include _data_size_with_guide_A.stan

#include _cont_hira_params_with_guide_A.stan
#include _disc_hira_params_with_guide_A.stan

// non-hirachical parameters
real mu_X;
vector<lower=0,upper=1>[2] p_D_given_R;
}

transformed data {
}

parameters {
// continous data
matrix[n_g,n_c] X; // log expression
}

transformed parameters {
matrix[n_r,n_c] K; // knockout strengths of sgRNAs and cells
for (c in 1:n_c) {
  for (r in 1:n_r) {
    K[r,c] = R[c,r];
  }
}
}

model {
// continous distributions
target += normal_lpdf(to_vector(X) | to_vector(rep_matrix(mu_X+mu_X_g, n_c) + gRNA_effects*K), to_vector(rep_matrix(sd_X_g, n_c)));
}

generated quantities{
#include _dicrete_data_with_guide_A.stan

// distributions of discrete data
for (c in 1:n_c) {
  for (g in 1:n_g) {
    Y[c,g] = poisson_rng(exp(X[g,c]+E_c[c]));
  }
}
for (c in 1:n_c) {
  for (r in 1:n_r) {
    D[c,r] = bernoulli_rng(p_D_given_R[R[c,r]?1:2]);
  }
}
}
