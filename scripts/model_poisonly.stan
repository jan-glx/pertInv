data {
int<lower=0> N_cells;
int<lower=0> N_genes;
int<lower=0> counts[N_genes*N_cells];
}

parameters {
matrix[N_genes,N_cells] ln_expression_noise;
row_vector[N_cells] ln_sf; // 'size factor' or 'library size'
vector<lower=0>[N_genes] gene_variance;
vector[N_genes] ln_mean_expression;
real<lower=0> sd_sf;
real gene_variace_mode;
real<lower=0> gene_variance_scale;
real<lower=0> sd_ln_mean_expression;
real mean_ln_mean_expression;
}

transformed parameters {

matrix[N_genes,N_cells] expression =
  exp(
    rep_matrix(mean_ln_mean_expression+ln_mean_expression, N_cells) +
    rep_matrix(ln_sf, N_genes) +
    ln_expression_noise
  );
}

model {

//target += normal_lpdf(gene_variace_mode | 0, 3);
//target += normal_lpdf(gene_variance_scale | 1, 10);
//target += normal_lpdf(gene_variance | gene_variace_mode, gene_variance_scale); //cauchy_lpdf
target += normal_lpdf(to_vector(ln_expression_noise) | 0.0, to_vector(rep_matrix(gene_variance, N_cells)));

//target += inv_gamma_lpdf(sd_sf | 40, 10);
//target += normal_lpdf(ln_sf | 0, sd_sf);

//target += inv_gamma_lpdf(sd_ln_mean_expression | 5, 5);
target += normal_lpdf(ln_mean_expression | 0, sd_ln_mean_expression);

target += normal_lpdf(mean_ln_mean_expression | 0, 10);

target += poisson_lpmf(counts  | to_vector(expression));
}

generated quantities{
/*int sim_counts[N_genes, N_cells];
for (g in 1:N_genes) {
  for (c in 1:N_cells) {
    sim_counts[g,c] = poisson_rng(expression[g,c]);
  }
}*/
}
