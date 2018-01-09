functions{
  #include /functions/index_set_diff.stan
  #include /functions/common.stan
}

data {
  #include /chunks/1_size.stan
}

transformed data {
  #include /chunks/1C_prior_params.stan
}

parameters {
}

transformed parameters {
}

model {
}

generated quantities{
  #include /chunks/1C_latent_cont.stan
  #include /chunks/1_latent_disc.stan
  #include /chunks/1_observed.stan
  row_vector[n_r] logit_p_R_r; // gRNA_prior_probs
  vector<lower=0,upper=1>[2] p_D_given_R;

  // hirachical prior parameters
  mu_log_p_R_r = normal_rng(mu_mu_log_p_R, sd_mu_log_p_R); // prior of mean of log abundance of sgRNAs
  sd_log_p_R_r = lognormal_rng(mu_log_sd_log_p_R_r, sd_log_sd_log_p_R_r); // prior of sd of log abundance of sgRNAs
  sd_mu_X = lognormal_rng(mu_log_sd_mu_X, sd_log_sd_mu_X);  // prior of sd of mean log expression of genes
  sd_E = lognormal_rng(mu_log_sd_E, sd_log_sd_E); // prior of sd of log size factors
  mu_log_sd_X_g = normal_rng(mu_mu_log_sd_X_g, sd_mu_log_sd_X_g); // prior of mean log variance of log expression of genes
  sd_log_sd_X_g = lognormal_rng(mu_log_sd_log_sd_X_g, sd_log_sd_log_sd_X_g);  // prior of sd of log variance of log expression of genes
  sd_gRNA_effects = lognormal_rng(mu_log_sd_gRNA_effects, sd_log_sd_gRNA_effects); // prior of sd of sgRNA effects

  // non-hirachical prior parameters
  mu_X = normal_rng(mu_muX, sd_muX);
  logit_p_D_given_R[1] = normal_rng(mu_logit_p_D_given_R1, sd_logit_p_D_given_R1);
  logit_p_D_given_R[2] = normal_rng(mu_logit_p_D_given_R2, sd_logit_p_D_given_R2); // not present, detected  +  present,no knockout, detected
  p_D_given_R = inv_logit(logit_p_D_given_R);

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
    logit_p_R_r[r] = log1m_exp(-p_R_r_helper[r])+p_R_r_helper[r]; // distribution of abundance of gRNAs
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
        R[c,r] = bernoulli_logit_rng(logit_p_R_r[r]);
        total = total+R[c,r];
      }
    }
  }


  // data
  for (c in 1:n_c) {
    for (g in 1:n_g) {
      Y[c,g] = neg_binomial_2_log_rng(mu_X + mu_X_g[g] + gRNA_effects[g] * to_vector(R[c]) + E_c[c],inv(sd_X_g[g]));
    }
  }
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      D[c,r] = bernoulli_logit_rng(logit_p_D_given_R[R[c,r]?1:2]);
    }
  }
}

