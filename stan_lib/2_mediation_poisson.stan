data {
  int<lower=0> n;
  int<lower=0> mA[n];
  int<lower=0> mB[n];
  int<lower=0> mC[n];
}

parameters {
  vector[n] A;
  vector[n] B;
  vector[n] C;
  real mu_A;
  real<lower=0> sd_A;
  real mu_B;
  real<lower=0> sd_B;
  real mu_C;
  real<lower=0> sd_C;
  real beta_AB;
  real beta_AC;
  real beta_BC;
}

model {
  target += normal_lpdf(A | mu_A, sd_A);
  target += normal_lpdf(B | mu_B + beta_AB * A, sd_B);
  target += normal_lpdf(C | mu_C + beta_AC * A + beta_BC * B, sd_C);
  target += poisson_log_lpmf(mA | A);
  target += poisson_log_lpmf(mB | B);
  target += poisson_log_lpmf(mC | C);
}
