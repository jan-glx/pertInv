library(data.table)
library(ggplot2)
library(MCMCglmm)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

n_cells=500
n_genes=50
Y = count_matrix[seq_len(n_cells),seq_len(n_genes)]
X = guide_matrix[seq_len(n_cells),]


prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,)))



dat <- melt(cbind(data.table(Y, keep.rownames=TRUE), data.table(X)),
            measure.vars=seq_len(n_genes)+1,
            value.name="expression_count",
            variable.name = "gene")
setnames(dat,"rn","cell_id")

fit <- MCMCglmm(expression_count ~ gene+cell_id, family = "poisson", data = dat,
                pl = TRUE)
summary(fit)

plot(fit$VCV)
fit2 <- MCMCglmm(expression_count ~ 1, random= ~ gene+cell_id, family = "poisson", data = dat)
summary(fit2)

fit3 <- MCMCglmm(expression_count ~ 1, random= ~ gene+cell_id+m_Stat3_3:gene, family = "poisson", data = dat)
summary(fit3)


fit4 <- MCMCglmm(Y ~ 1, random= ~ trait+unit+m_Stat3_3:trait, family = "poisson", data = data.table(X))
summary(fit4)

fit5 <- MCMCglmm(expression_count ~ -1 +garbage, random= ~ m_Stat3_3:gene, rcov = ~ us(gene):cell_id, family = "poisson", data = copy(dat)[,garbage:=rnorm(.N)])
summary(fit5)

fit6 <- MCMCglmm(expression_count ~ 1, random= ~ us(gene)+cell_id+m_Stat3_3:gene, rcov = ~ us(gene):cell_id, family = "poisson", data = dat)
summary(fit6)
