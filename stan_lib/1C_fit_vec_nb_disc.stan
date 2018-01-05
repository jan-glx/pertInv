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
  #include 1C_prior_params.stan
  int R_2m[n_c,n_r];
  for (c in 1:n_c) {
    for (r in 1:n_r) {
      R_2m[c,r] =2-R[c,r];
    }
  }
}

parameters {
  #include 1C_latent_cont.stan
}

transformed parameters {
  row_vector[n_r] logit_p_R_r; // gRNA_prior_probs

  for (r in 1:n_r) {
    logit_p_R_r[r] =log1m_exp(-p_R_r_helper[r])+p_R_r_helper[r];// log(x/(1-x))=log(x)-log(1-x)//logit(1-exp(-p_R_r_helper[r])); // distribution of abundance of gRNAs
  }
}

model {
  // hirachical prior parameters
  target += normal_lpdf(mu_log_p_R_r | mu_mu_log_p_R, sd_mu_log_p_R); // prior of mean of log abundance of sgRNAs
  target += lognormal_lpdf(sd_log_p_R_r | mu_log_sd_log_p_R_r, sd_log_sd_log_p_R_r); // prior of sd of log abundance of sgRNAs
  target += lognormal_lpdf(sd_mu_X | mu_log_sd_mu_X, sd_log_sd_mu_X);  // prior of sd of mean log expression of genes
  target += lognormal_lpdf(sd_E | mu_log_sd_E, sd_log_sd_E); // prior of sd of log size factors
  target += normal_lpdf(mu_log_sd_X_g | mu_mu_log_sd_X_g, sd_mu_log_sd_X_g); // prior of mean log variance of log expression of genes
  target += lognormal_lpdf(sd_log_sd_X_g | mu_log_sd_log_sd_X_g, sd_log_sd_log_sd_X_g);  // prior of sd of log variance of log expression of genes
  target += lognormal_lpdf(sd_gRNA_effects | mu_log_sd_gRNA_effects, sd_log_sd_gRNA_effects); // prior of sd of sgRNA effects

  // non-hirachical prior parameters
  target += normal_lpdf(mu_X | mu_muX, sd_muX);
  target += normal_lpdf(logit_p_D_given_R[1] | mu_logit_p_D_given_R1, sd_logit_p_D_given_R1) ;
  target += normal_lpdf(logit_p_D_given_R[2] | mu_logit_p_D_given_R2, sd_logit_p_D_given_R2) ;

  // continous hirachical parameters
  target += normal_lpdf(to_vector(gRNA_effects) | 0, sd_gRNA_effects); // distribution of effects of gRNAs on genes
  target += lognormal_lpdf(sd_X_g | mu_log_sd_X_g, sd_log_sd_X_g); // distribution of expression variance of genes
  target += lognormal_lpdf(p_R_r_helper | mu_log_p_R_r, sd_log_p_R_r); // distribution of abundance of gRNAs
  target += normal_lpdf(E_c | 0, sd_E); // distribution of size factors of cells
  target += normal_lpdf(mu_X_g | 0, sd_mu_X); // distribution of mean expression of genes

  // discrete latent variables
  target += bernoulli_logit_lpmf(to_array_1d(R) | to_vector(rep_matrix(logit_p_R_r, n_c))); // distribution of guide presence
  // TODO: selection bias against Sum(R)==0 cells

  // further continous latent variables
  //target += normal_lpdf(to_vector(K) | to_vector(to_array_1d(R)), 0.1); // Knockout strength

  // data
  target += neg_binomial_2_log_lpmf(to_array_1d(Y)  | to_vector(rep_matrix(mu_X+mu_X_g, n_c) + gRNA_effects*to_matrix(to_vector(to_array_1d(R)), n_r, n_c) + rep_matrix(E_c, n_g)), inv(to_vector(rep_matrix(sd_X_g, n_c))));
  target += bernoulli_logit_lpmf(to_array_1d(D) | logit_p_D_given_R[to_array_1d(R_2m)]);
}

generated quantities{
  // discrete latent variables
  matrix[n_r,n_c] lp_R_terms;
  vector<lower=0,upper=1>[2] p_D_given_R;

  for (c in 1:n_c) {
    for (r in 1:n_r) {
      lp_R_terms[r,c] = bernoulli_logit_lpmf(R[c,r] | logit_p_R_r[r]); // distribution of guide presence
      // TODO: selection bias against Sum(R)==0 cells
    }
  }
  p_D_given_R = inv_logit(logit_p_D_given_R);
}
