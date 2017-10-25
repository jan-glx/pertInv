// continous hirachical parameters
matrix[n_g,n_r] gRNA_effects; // total effects of gRNAs
vector<lower=0>[n_g] sd_X_g; // std. dev.s of genes' expression
row_vector<lower=0,upper=1>[n_r] p_R_r; // gRNA_prior_probs
row_vector[n_c] E_c; // log exposure 'size factor' or 'library size'
vector[n_g] mu_X_g; // means of genes' expression
