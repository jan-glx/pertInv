
# ---------------------------
library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix[1:100,1:100]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)


library(rstan)
options(mc.cores = parallel::detectCores())
# ---------------------------
stanmodelcode <- "
data {
int<lower=0> N_cells;
int<lower=0> N_genes;
int<lower=0> counts[N_genes*N_cells];
}

parameters {
vector[N_genes*N_cells] ln_expression_noise;
row_vector[N_cells] ln_sf; // 'size factor' or 'library size'
vector<lower=0>[N_genes] gene_variance;
vector[N_genes] ln_mean_expression;
real<lower=0> alpha_sf;
real<lower=0> sd_sf;
//real<lower=0> alpha_gene_variance;
real<lower=0> beta_gene_variance;
real<lower=0> sd_ln_mean_expression;
real mean_ln_mean_expression;
real intercept;
}

model {
//alpha_gene_variance = 1.0;
target += gamma_lpdf(beta_gene_variance | 4, 10);
target += inv_gamma_lpdf(gene_variance | 1.0, beta_gene_variance);
target += normal_lpdf(ln_expression_noise | 0.0, to_vector(rep_matrix(gene_variance, N_cells)));

target += inv_gamma_lpdf(sd_sf | 40, 10);
target += normal_lpdf(ln_sf | 0, sd_sf);

target += inv_gamma_lpdf(sd_ln_mean_expression | 5, 5);
target += normal_lpdf(ln_mean_expression | 0, sd_ln_mean_expression);

target += normal_lpdf(intercept | -1, 3);
target += poisson_lpmf(counts  | exp(to_vector(rep_matrix(intercept+ln_mean_expression, N_cells) + rep_matrix(ln_sf, N_genes)) + ln_expression_noise));
}
"

dat <- list(N_cells = nrow(count_matrix), N_genes = ncol(count_matrix), counts = as.vector(t(count_matrix)));
fit <- stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 2012, chains = 3, sample_file = 'norm.csv',
         #   control = list(adapt_delta = 0.99),
            verbose = TRUE)

# ---------------------------------------
print(fit)
traceplot(fit)

# extract samples
e <- extract(fit, permuted = TRUE) # return a list of arrays
mu <- e$mu

m <- extract(fit, permuted = FALSE, inc_warmup = FALSE) # return an array
print(dimnames(m))

# using as.array directly on stanfit objects
m2 <- as.array(fit)

pairs(fit, pars = c("intercept","ln_mean_expression[1]","ln_expression_noise[1]","sd_sf","ln_sf[1]","beta_gene_variance","gene_variance[1]"))

sf_ <-  exp(apply(m,3,mean)[1000+seq_len(n)])


dt = data.table(melt(count_matrix))
setnames(dt, c("cell","gene","counts"))
dt[,total_counts:=sum(counts),by=.(cell)]

nf_ <- edgeR::calcNormFactors(t(count_matrix))

dt[,sf:=sf_[.GRP],by=.(cell)]
dt[,nf:=nf_[.GRP],by=.(cell)]

summary_genes = rbind(dt[,.(mean=mean(counts),var=var(counts),method="plain"),by=.(gene)],
                      dt[,.(mean=mean(counts/total_counts*mean(total_counts)),var=var(counts/total_counts*mean(total_counts)),method="total_counts normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts/sf),var=var(counts/sf),method="MCMCglmm normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts*nf),var=var(counts*nf),method="nf normalized"),by=.(gene)])

figure("mean-variance relationship of genes stan",
       ggplot(summary_genes,aes(x=mean,y=var,color=method))+geom_point(alpha=0.3)+
         scale_x_log10()+scale_y_log10()+geom_abline(color="black")+geom_smooth()+
         theme(legend.justification = c(1, 0), legend.position = c(1, 0))
)
