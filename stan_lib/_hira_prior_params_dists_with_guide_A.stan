
// continous hirachical prior parameters distributions
// all priors are waekly informative
target += normal_lpdf(mu_log_p_R_r | log(1.0 / n_r), log(1.3)); // prior of mean of log abundance of sgRNAs
target += cauchy_lpdf(sd_log_p_R_r | 0, log(1.3)); // prior of sd of log abundance of sgRNAs
target += cauchy_lpdf(sd_mu_X | 0, 1);  // prior of sd of mean log expression of genes
target += cauchy_lpdf(sd_E | 0, 0.2); // prior of sd of log size factors
target += normal_lpdf(mu_log_var_X_g | log(1.0), 1); // prior of mean log variance of log expression of genes
target += cauchy_lpdf(sd_log_var_X_g | 0, log(1.3));  // prior of sd of log variance of log expression of genes
target += cauchy_lpdf(sd_gRNA_effects | 0, 0.3); // prior of sd of sgRNA effects
