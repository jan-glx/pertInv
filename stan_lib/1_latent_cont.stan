
// hirachical prior parameters
real mu_log_p_R_r;
real<lower=0> sd_log_p_R_r;
real<lower=0> sd_mu_X;
real<lower=0> sd_E;
real mu_log_sd_X_g;
real<lower=0> sd_log_sd_X_g;
real<lower=0> sd_gRNA_effects;

// non-hirachical prior parameters
real mu_X;
vector<lower=0,upper=1>[2] p_D_given_R;

// continous hirachical parameters
matrix[n_g,n_r] gRNA_effects; // total effects of gRNAs
vector<lower=0>[n_g] sd_X_g; // std. dev.s of genes' expression
row_vector<lower=0>[n_r] p_R_r_helper; // gRNA_prior_probs
row_vector[n_c] E_c; // log exposure 'size factor' or 'library size'
vector[n_g] mu_X_g; // means of genes' expression

// further continous latent variables
matrix[n_g,n_c] X_noise; // log expression noise
matrix[n_r,n_c] K; // knockout strengths of sgRNAs and cells
