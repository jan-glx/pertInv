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
}

parameters {
  #include 1_latent_cont.stan
}

transformed parameters {
  row_vector<lower=0,upper=1>[n_r] p_R_r; // gRNA_prior_probs

  for (r in 1:n_r) {
    p_R_r[r] = 1-exp(-p_R_r_helper[r]); // distribution of abundance of gRNAs
  }
}

model {
  // hirachical prior parameters
  mu_log_p_R_r ~ normal(log(1.0 / n_r), log(1.3)); // prior of mean of log abundance of sgRNAs
  sd_log_p_R_r ~ lognormal(-1, 0.1); // prior of sd of log abundance of sgRNAs
  sd_mu_X ~ lognormal(-1, 0.1);  // prior of sd of mean log expression of genes
  sd_E ~ lognormal(-1, 0.1); // prior of sd of log size factors
  mu_log_var_X_g ~ normal(-1, 0.1); // prior of mean log variance of log expression of genes
  sd_log_var_X_g ~ lognormal(-1, 0.1);  // prior of sd of log variance of log expression of genes
  sd_gRNA_effects ~ lognormal(-3, 0.1); // prior of sd of sgRNA effects

  // non-hirachical prior parameters
  mu_X ~ normal(0, 1);
  p_D_given_R[1] ~ beta(100, 1);
  p_D_given_R[2] ~ beta(10, 1000);

  // continous hirachical parameters
  for (r in 1:n_r) {
    for (g in 1:n_g) {
      gRNA_effects[g,r] ~ normal(0, sd_gRNA_effects); // distribution of effects of gRNAs on genes
    }
  }
  for (g in 1:n_g) {
    sd_X_g[g] ~ lognormal(mu_log_var_X_g, sd_log_var_X_g); // distribution of expression variance of genes
  }
  for (r in 1:n_r) {
    p_R_r_helper[r] ~ lognormal(mu_log_p_R_r, sd_log_p_R_r); // distribution of abundance of gRNAs
  }
  for (c in 1:n_c) {
    E_c[c] ~ normal(0, sd_E); // distribution of size factors of cells
  }
  for (g in 1:n_g) {
    mu_X_g[g] ~ normal(0, sd_mu_X); // distribution of mean expression of genes
  }

  // discrete latent variables
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      R[c,r] ~ bernoulli(p_R_r[r]);
    }
  }

  // further continous latent variables
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      K[r,c] ~ normal(R[c,r], 0.1);
    }
  }

  for (c in 1:n_c) {
    for (g in 1:n_g) {
      X[g,c] ~ normal(mu_X + mu_X_g[g] + gRNA_effects[g]*col(K,c), sd_X_g[g]);
    }
  }

  // data
  for (c in 1:n_c) {
    for (g in 1:n_g) {
      Y[c,g] ~ poisson(exp(X[g,c]+E_c[c]));
    }
  }
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      D[c,r] ~ bernoulli(p_D_given_R[R[c,r]?1:2]);
    }
  }
}

generated quantities{
}
