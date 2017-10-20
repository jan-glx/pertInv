data {
int<lower=0> N_cells;
int<lower=0> N_genes;
int<lower=0> N_gRNAs;
int<lower=0> counts[N_genes*N_cells];
int<lower=0,upper=1> guide_detected[N_cells,N_gRNAs];
}

parameters {
vector[N_genes*N_cells] ln_expression_noise;
matrix[N_gRNAs,N_cells] knockout;
matrix[N_genes,N_gRNAs] knockout_effect;
vector<lower=0>[N_gRNAs] lambda_guide;
row_vector[N_cells] ln_sf; // 'size factor' or 'library size'
vector<lower=0>[N_genes] gene_variance;
vector[N_genes] ln_mean_expression;
real<lower=0> alpha_sf;
real<lower=0> sd_sf;
//real<lower=0> alpha_gene_variance;
real<lower=0> beta_gene_variance;
real<lower=0> sd_ln_mean_expression;
real mean_ln_mean_expression;
vector<lower=0,upper=1>[N_gRNAs] detection_prob;
vector<lower=0,upper=1>[N_gRNAs] knockout_rate;
real<lower=0,upper=1> false_detection_prob;
real<lower=0,upper=1> pos_selection_prob;
real<lower=0,upper=1> neg_selection_prob;
vector[N_gRNAs] mean_knockout;
vector<lower=0>[N_gRNAs] sd_knockout;
}

transformed parameters {
}

model {

real lpdf_detected_present_true;
real lpdf_ko_detected_present_true;
real lpdf_detected_present_false;
real lpdf_ko_detected_present_false;
real prod_term;
real prod_of_sums_term;
real prod_term_norm;
real prod_of_sums_term_norm;
matrix[N_gRNAs, N_cells] llr_guide_present;

//alpha_gene_variance = 1.0;
target += gamma_lpdf(beta_gene_variance | 4, 10);
target += inv_gamma_lpdf(gene_variance | 1.0, beta_gene_variance);
target += normal_lpdf(ln_expression_noise | 0.0, to_vector(rep_matrix(gene_variance, N_cells)));

target += inv_gamma_lpdf(sd_sf | 40, 10);
target += normal_lpdf(ln_sf | 0, sd_sf);

target += inv_gamma_lpdf(sd_ln_mean_expression | 5, 5);
target += normal_lpdf(ln_mean_expression | 0, sd_ln_mean_expression);

target += normal_lpdf(mean_ln_mean_expression | 0, 10);


target += normal_lpdf(mean_knockout | 0, 2);
target += inv_gamma_lpdf(sd_knockout | 5, 5);

target += inv_gamma_lpdf(lambda_guide | 10, (10.0-1)/N_gRNAs);
target += beta_lpdf(detection_prob | 7, 3);
target += beta_lpdf(knockout_rate | 7, 3);

target += beta_lpdf(false_detection_prob | 1, 100*N_gRNAs);
target += beta_lpdf(neg_selection_prob | 1, 100);
target += beta_lpdf(pos_selection_prob | 100, 10);
target += normal_lpdf(to_vector(knockout_effect) | 0, 1);


for (c in 1:N_cells) {
  prod_of_sums_term = log(pos_selection_prob);
  prod_term = log(neg_selection_prob-pos_selection_prob);
  prod_of_sums_term_norm = prod_of_sums_term;
  prod_term_norm = prod_term;
  for (r in 1:N_gRNAs) {

    lpdf_detected_present_true = bernoulli_lpmf(guide_detected[c,r] | detection_prob[r]) +
      log1m(exp(-lambda_guide[r]));

    lpdf_ko_detected_present_true =
      log_mix(knockout_rate[r],
              normal_lpdf(knockout[r,c] | mean_knockout[r], sd_knockout[r]),
              normal_lpdf(knockout[r,c] | 0, 0.001)
      ) + lpdf_detected_present_true;

    lpdf_detected_present_false =
      bernoulli_lpmf(guide_detected[c,r]| false_detection_prob) - lambda_guide[r];
    lpdf_ko_detected_present_false =
      normal_lpdf(knockout[r,c] | 0, 0.001) + lpdf_detected_present_false;

    prod_of_sums_term =
      prod_of_sums_term + log_sum_exp(lpdf_ko_detected_present_false, lpdf_ko_detected_present_true);
    prod_term = prod_term + lpdf_ko_detected_present_false;

    prod_of_sums_term_norm =
      prod_of_sums_term_norm + log_sum_exp(lpdf_detected_present_false, lpdf_detected_present_true);
    prod_term_norm = prod_term_norm + lpdf_detected_present_false;

    llr_guide_present[r,c] = lpdf_ko_detected_present_true + log(pos_selection_prob) -
      (lpdf_ko_detected_present_false  +
       log_mix(exp(-(sum(lambda_guide)-lambda_guide[r])), neg_selection_prob, pos_selection_prob)
      );
  }
  target += log_sum_exp(prod_of_sums_term, prod_of_sums_term);
  target += -log_sum_exp(prod_of_sums_term_norm, prod_term_norm);
}



target += poisson_lpmf(counts  | exp(
  to_vector(rep_matrix(mean_ln_mean_expression+ln_mean_expression, N_cells) +
  rep_matrix(ln_sf, N_genes)) +
  ln_expression_noise +
  to_vector(knockout_effect * knockout)
));
}

generated quantities{
matrix[N_gRNAs, N_cells] llr_guide_present;
{
  real lpdf_detected_present_true;
  real lpdf_ko_detected_present_true;
  real lpdf_detected_present_false;
  real lpdf_ko_detected_present_false;
  real prod_term;
  real prod_of_sums_term;
  real prod_term_norm;
  real prod_of_sums_term_norm;

  for (c in 1:N_cells) {
    prod_of_sums_term = log(pos_selection_prob);
    prod_term = log(neg_selection_prob-pos_selection_prob);
    prod_of_sums_term_norm = prod_of_sums_term;
    prod_term_norm = prod_term;
    for (r in 1:N_gRNAs) {

      lpdf_detected_present_true = bernoulli_lpmf(guide_detected[c,r] | detection_prob[r]) + log1m(exp(-lambda_guide[r]));
      lpdf_ko_detected_present_true = log_mix(knockout_rate[r], normal_lpdf(knockout[r,c] | mean_knockout[r], sd_knockout[r]),normal_lpdf(knockout[r,c] | 0, 0.001)) + lpdf_detected_present_true;


      lpdf_detected_present_false = bernoulli_lpmf(guide_detected[c,r]| false_detection_prob) -lambda_guide[r];
      lpdf_ko_detected_present_false = normal_lpdf(knockout[r,c] | 0, 0.01) + lpdf_detected_present_false;

      prod_of_sums_term = prod_of_sums_term + log_sum_exp(lpdf_ko_detected_present_false, lpdf_ko_detected_present_true);
      prod_term = prod_term + lpdf_ko_detected_present_false;

      prod_of_sums_term_norm = prod_of_sums_term_norm + log_sum_exp(lpdf_detected_present_false, lpdf_detected_present_true);
      prod_term_norm = prod_term_norm + lpdf_detected_present_false;

      llr_guide_present[r,c] = lpdf_ko_detected_present_true + log(pos_selection_prob) -
        (lpdf_ko_detected_present_false  +  log_mix(exp(-sum(lambda_guide)), neg_selection_prob, pos_selection_prob));
    }
  }
}
}
