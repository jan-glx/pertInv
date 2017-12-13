functions{
  #include index_set_diff.stan
  #include common.stan
}

data {
  #include 1_size.stan
}

transformed data {
}

parameters {
}

transformed parameters {
}

model {
}

generated quantities{
  #include 1B_latent_cont.stan
  #include 1_latent_disc.stan
  #include 1_observed.stan
  row_vector<lower=0,upper=1>[n_r] p_R_r; // gRNA_prior_probs

  // hirachical prior parameters
  mu_log_p_R_r = normal_rng(log(1.0 / n_r), log(1.3)); // prior of mean of log abundance of sgRNAs
  sd_log_p_R_r = lognormal_rng(-1, 0.1); // prior of sd of log abundance of sgRNAs
  sd_mu_X = lognormal_rng(-1, 0.1);  // prior of sd of mean log expression of genes
  sd_E = lognormal_rng(-1, 0.1); // prior of sd of log size factors
  mu_log_sd_X_g = normal_rng(-1, 0.1); // prior of mean log variance of log expression of genes
  sd_log_sd_X_g = lognormal_rng(-1, 0.1);  // prior of sd of log variance of log expression of genes
  sd_gRNA_effects = lognormal_rng(-2, 0.1); // prior of sd of sgRNA effects

  // non-hirachical prior parameters
  mu_X = normal_rng(0, 1);
  p_D_given_R[1] = beta_rng(100, 1);
  p_D_given_R[2] = beta_rng(10, 1000);

  // continous hirachical parameters
  for (r in 1:n_r) {
    for (g in 1:n_g) {
      gRNA_effects[g,r] = normal_rng(0, sd_gRNA_effects); // distribution of effects of gRNAs on genes
    }
  }
  for (g in 1:n_g) {
    sd_X_g[g] = lognormal_rng(mu_log_sd_X_g, sd_log_sd_X_g); // distribution of expression variance of genes
  }
  for (r in 1:n_r) {
    p_R_r_helper[r] = lognormal_rng(mu_log_p_R_r, sd_log_p_R_r); // distribution of abundance of gRNAs
  }
  for (r in 1:n_r) {
    p_R_r[r] = 1-exp(-p_R_r_helper[r]); // distribution of abundance of gRNAs
  }
  for (c in 1:n_c) {
    E_c[c] = normal_rng(0, sd_E); // distribution of size factors of cells
  }
  for (g in 1:n_g) {
    mu_X_g[g] = normal_rng(0, sd_mu_X); // distribution of mean expression of genes
  }

  // discrete latent variables
  for (c in 1:n_c) {
    int total = 0;
    while (total<=0) {
      for (r in 1:n_r) {
        R[c,r] = bernoulli_rng(p_R_r[r]);
        total = total+R[c,r];
      }
    }
  }

  // further continous latent variables
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      K[r,c] = normal_rng(R[c,r], 0.1);
    }
  }


  // data
  for (c in 1:n_c) {
    for (g in 1:n_g) {
      Y[c,g] = neg_binomial_2_log_rng(mu_X + mu_X_g[g] + gRNA_effects[g] * col(K,c)+E_c[c], 1/sd_X_g[g]);
    }
  }
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      D[c,r] = bernoulli_rng(p_D_given_R[R[c,r]?1:2]);
    }
  }
}
