// continous hirachical parameters distributions
target += normal_lpdf(to_vector(gRNA_effects) | 0, sd_gRNA_effects); // distribution of effects of gRNAs on genes
target += lognormal_lpdf(sd_X_g | mu_log_var_X_g, sd_log_var_X_g); // distribution of expression variance of genes
target += lognormal_lpdf(p_R_r | mu_log_p_R_r, sd_log_p_R_r); // distribution of abundance of gRNAs
target += normal_lpdf(E_c | 0, sd_E); // distribution of size factors of cells
target += normal_lpdf(mu_X_g | 0, sd_mu_X); // distribution of mean expression of genes
