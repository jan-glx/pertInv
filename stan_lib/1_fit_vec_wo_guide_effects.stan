functions{
  #include index_set_diff.stan
  #include common.stan
}

data {
  #include 1_size.stan
  #include 1_observed.stan

  #include 1_latent_disc.stan
}

transformed data {
  int R_2m[n_c,n_r];
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      R_2m[c,r] =2-R[c,r];
    }
  }
}

parameters {
  #include 1_latent_cont.stan
}

transformed parameters {
  row_vector<lower=0,upper=1>[n_r] p_R_r; // gRNA_prior_probs
  matrix[n_g,n_c] X; // log expression

  for (r in 1:n_r) {
    p_R_r[r] = 1-exp(-p_R_r_helper[r]); // distribution of abundance of gRNAs
  }

  X = rep_matrix(mu_X+mu_X_g, n_c) + gRNA_effects*K + X_noise .* rep_matrix(sd_X_g, n_c);
}

model {
  // hirachical prior parameters
  target += normal_lpdf(mu_log_p_R_r | log(1.0 / n_r), log(1.3)); // prior of mean of log abundance of sgRNAs
  target += lognormal_lpdf(sd_log_p_R_r | -1, 1); // prior of sd of log abundance of sgRNAs
  target += lognormal_lpdf(sd_mu_X | -1, 1);  // prior of sd of mean log expression of genes
  target += lognormal_lpdf(sd_E | -1, 1); // prior of sd of log size factors
  target += normal_lpdf(mu_log_sd_X_g | -1, 1); // prior of mean log variance of log expression of genes
  target += lognormal_lpdf(sd_log_sd_X_g | -1, 1);  // prior of sd of log variance of log expression of genes
  target += lognormal_lpdf(sd_gRNA_effects | -3, 1); // prior of sd of sgRNA effects

  // non-hirachical prior parameters
  target += normal_lpdf(mu_X | 0, 10);
  target += beta_lpdf(p_D_given_R[1] | 10, 0.1);
  target += beta_lpdf(p_D_given_R[2] | 1, 100);

  // continous hirachical parameters
  target += normal_lpdf(to_vector(gRNA_effects) | 0, sd_gRNA_effects); // distribution of effects of gRNAs on genes
  target += lognormal_lpdf(sd_X_g | mu_log_sd_X_g, sd_log_sd_X_g); // distribution of expression variance of genes
  target += lognormal_lpdf(p_R_r_helper | mu_log_p_R_r, sd_log_p_R_r); // distribution of abundance of gRNAs
  target += normal_lpdf(E_c | 0, sd_E); // distribution of size factors of cells
  target += normal_lpdf(mu_X_g | 0, sd_mu_X); // distribution of mean expression of genes

  // discrete latent variables
  target += bernoulli_lpmf(to_array_1d(R) | to_vector(rep_matrix(p_R_r, n_c)));

  // further continous latent variables
  target += normal_lpdf(to_vector(K) | to_vector(to_array_1d(R)), 0.1);
  target += normal_lpdf(to_vector(X_noise) | 0, 1);

  // data
  target += poisson_lpmf(to_array_1d(Y)  | to_vector(exp(X + rep_matrix(E_c, n_g))));
  target += bernoulli_lpmf(to_array_1d(D) | p_D_given_R[to_array_1d(R_2m)]);
}

generated quantities{
}
