
# ---------------------------
library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix[1:100,1:100]#count_matrix#count_matrix[1:1000,1:100]
guide_matrix <- guide_matrix[1:100,]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)


library(rstan)
options(mc.cores = parallel::detectCores())






# ---------------------------
m1 <-  stan_model("scripts/model.stan")
# ======

#dat <- list(N_cells = nrow(count_matrix), N_genes = ncol(count_matrix), N_guides = ncol(guide_matrix), counts = as.vector(t(count_matrix)),
#            guide_detected =  lapply(1:ncol(guide_matrix), function(i) as.numeric(guide_matrix[,i])));
dat <- list(N_cells = nrow(count_matrix), N_genes = ncol(count_matrix), N_guides = ncol(guide_matrix), counts = as.vector(t(count_matrix)),
            guide_detected =  1*guide_matrix);

#fit <- stan(model_code = stanmodelcode, model_name = "example",
#            data = dat, iter = 12, chains = 3, sample_file = 'norm.csv',
#         #   control = list(adapt_delta = 0.99),
#            verbose = TRUE)
fit <- vb(m1, data = dat, sample_file = 'norm2.csv')
# ---------------------------------------
print(fit)
stan_trace(fit, pars = c("mean_ln_mean_expression","ln_mean_expression[1]","ln_expression_noise[1]","sd_sf","ln_sf[1]","beta_gene_variance","gene_variance[1]"))

# extract samples
e <- extract(fit, permuted = TRUE, inc_warmup = FALSE) # return a list of arrays
ggplot(melt(e$ln_mean_expression),aes(x=as.factor(Var2),y=value))+geom_boxplot()
ggplot(melt(e$ln_sf),aes(x=as.factor(Var2),y=value))+geom_boxplot()
ggplot(melt(e$llr_guide_present[,,1]),aes(x=as.factor(Var2),y=value))+geom_boxplot()

library(loo)

log_lik_1 <- extract_log_lik(fit)
loo_1 <- loo(log_lik_1)
print(loo_1)

sf_ <-  exp(colMeans(e$ln_sf))


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
