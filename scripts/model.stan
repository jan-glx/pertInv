data {
int<lower=0> N_cells;
int<lower=0> N_genes;
int<lower=0> N_guides;
int<lower=0> counts[N_genes*N_cells];
int<lower=0,upper=1> guide_detected[N_cells,N_guides];
}

parameters {
vector[N_genes*N_cells] ln_expression_noise;
matrix[N_guides,N_cells] knockout;
matrix[N_genes,N_guides] knockout_effect;
vector<lower=0>[N_guides] lambda_guide;
row_vector[N_cells] ln_sf; // 'size factor' or 'library size'
vector<lower=0>[N_genes] gene_variance;
vector[N_genes] ln_mean_expression;
real<lower=0> alpha_sf;
real<lower=0> sd_sf;
//real<lower=0> alpha_gene_variance;
real<lower=0> beta_gene_variance;
real<lower=0> sd_ln_mean_expression;
real mean_ln_mean_expression;
vector<lower=0,upper=1>[N_cells] remaining_time;
vector<lower=0,upper=1>[N_guides] detection_prob;
real<lower=0,upper=1> false_detection_prob;
vector[N_guides] mean_knockout;
vector<lower=0>[N_guides] sd_knockout;
}

transformed parameters {
real<lower=0> T;
T = sum(lambda_guide);
}

model {

vector[N_guides] sum_term;
vector[N_guides] sum_term_norm;
vector[N_guides] sum_term_present_true;
real lpdf_detected_given_present_true;
real lpdf_ko_detected_given_present_true;
real lpdf_detected_given_present_false;
real lpdf_ko_detected_given_present_false;
real tmp;

//alpha_gene_variance = 1.0;
target += gamma_lpdf(beta_gene_variance | 4, 10);
target += inv_gamma_lpdf(gene_variance | 1.0, beta_gene_variance);
target += normal_lpdf(ln_expression_noise | 0.0, to_vector(rep_matrix(gene_variance, N_cells)));

target += inv_gamma_lpdf(sd_sf | 40, 10);
target += normal_lpdf(ln_sf | 0, sd_sf);

target += inv_gamma_lpdf(sd_ln_mean_expression | 5, 5);
target += normal_lpdf(ln_mean_expression | 0, sd_ln_mean_expression);

target += normal_lpdf(mean_ln_mean_expression | 0, 10);


target += normal_lpdf(mean_knockout | 0, 10);
target += inv_gamma_lpdf(sd_knockout | 5, 5);

target += exponential_lpdf(1-remaining_time| T);

target += inv_gamma_lpdf(lambda_guide | 10, (10.0-1)/N_guides);
target += beta_lpdf(detection_prob | 6, 3);

target += beta_lpdf(false_detection_prob | 1, 100*N_guides);


for (i in 1:N_cells) {
  for (k in 1:N_guides) {

    lpdf_detected_given_present_true = bernoulli_lpmf(guide_detected[i,k] | detection_prob[k]);
    lpdf_ko_detected_given_present_true = normal_lpdf(knockout[k,i] | mean_knockout[k], sd_knockout[k]) + lpdf_detected_given_present_true;

    lpdf_detected_given_present_false = bernoulli_lpmf(guide_detected[i,k]| false_detection_prob);
    lpdf_ko_detected_given_present_false = normal_lpdf(knockout[k,i] | 0, 0.001) +lpdf_detected_given_present_false;
    tmp = log_mix(exp(-lambda_guide[k] * remaining_time[i]),
                      lpdf_ko_detected_given_present_true,
                      lpdf_ko_detected_given_present_false);
    target += tmp;
    sum_term[k] = log(lambda_guide[k])-log(T)+lpdf_ko_detected_given_present_true-tmp;

    tmp = log_mix(exp(-lambda_guide[k] * remaining_time[i]),
                  lpdf_detected_given_present_true,
                  lpdf_detected_given_present_false);
    target += -tmp; // for normalization
    sum_term_norm[k] = log(lambda_guide[k])-log(T)+lpdf_detected_given_present_true-tmp;
  }
  target += log_sum_exp(sum_term);
  target += -log_sum_exp(sum_term_norm);
}



target += poisson_lpmf(counts  | exp(
  to_vector(rep_matrix(mean_ln_mean_expression+ln_mean_expression, N_cells) +
  rep_matrix(ln_sf, N_genes)) +
  ln_expression_noise +
  to_vector(knockout_effect * knockout)
));
}


generated quantities {
matrix[N_guides, N_cells] llr_guide_present;
real p_g;

for (k in 1:N_guides) {
  for (i in 1:N_cells) {
    p_g = lambda_guide[k]/T+(1-lambda_guide[k]/T)*exp(-lambda_guide[k] * remaining_time[i]);
    llr_guide_present[k,i] = bernoulli_lpmf(guide_detected[i,k]| 0.00001) + normal_lpdf(knockout[k,i] | 0, 0.001) + log(p_g)-
      (log1m(p_g) + bernoulli_lpmf(guide_detected[i,k] | detection_prob[k]) + normal_lpdf(knockout[k,i] | mean_knockout[k], sd_knockout[k]));
  }
}
}
