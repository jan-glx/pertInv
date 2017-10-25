functions{
#include index_set_diff.stan
#include common.stan
}

data {
#include _data_size_with_guide_A.stan

#include _hira_prior_params_with_guide_A.stan

}

transformed data {
}

parameters {
#include _cont_hira_params_with_guide_A.stan
}

transformed parameters {
}

model {
#include _cont_hira_dists_with_guide_A.stan
}

generated quantities{
#include _disc_hira_params_with_guide_A.stan
// discrete hirachical parameters distributions
for (c in 1:n_c) {
  int total = 0;
  while (total<=0) {
    for (r in 1:n_r) {
      R[c,r] = bernoulli_rng(p_R_r[r]);
      total = total+R[c,r];
    }
  }
}
}
