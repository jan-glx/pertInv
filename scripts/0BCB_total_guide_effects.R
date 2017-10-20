covariates.dt <- cells_summary.dt[cellnames.dt, on="cell_id"]
Y <- stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt)

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)


guide_matrix <- cbc_gbc_dict.dt[cells_summary.dt[cellnames.dt, on="cell_id"], on="cell"
                                ][!is.na(guide_id), as.matrix(Matrix::sparseMatrix(cell_id2, guide_id, x=TRUE))]
colnames(guide_matrix) <- cbc_gbc_dict.dt[,guide[1],keyby=guide_id][[2]]

head(guide_matrix)

res <- lm(Y_adj~., as.data.table(guide_matrix))
aov(res)

res <- psych::corr.test(guide_matrix, Y, adjust = "BH")
minp <- arrayInd(which.min(res$p), dim(res$p))
minv <- min(res$p)
figure("most significant total guide effect",
ggplot(data.table(x=guide_matrix[,minp[1]], y=Y[,minp[2]]), aes(x=x,y=y))+
  geom_violin()+xlab(colnames(guide_matrix)[minp[1]])+ylab(colnames(Y)[minp[2]])
)

image(res$p)
plot_matrix(t(res$p))
plot(ecdf(res$p))
correlations.dt <- data.table(psych::corr.test(guide_matrix, Y, adjust = "none")$p,keep.rownames = TRUE)
correlations.dt <- melt(correlations.dt, id.vars = "rn")
ggplot(correlations.dt, aes(x=value))+stat_ecdf()+scale_x_log10()+scale_y_log10()


trans_neg_log10 <- scales::trans_new("neg_log10", function(x) -log10(x), function(y) 10^(-y), domain=c(0, Inf))

p_real_null_plot <- function(test, data) {
  p_real <- test(data)
  p_null <- test(data[sample(nrow(data)),])
  ggplot(data.table(resampled = -log10(sort(p_null)), real=-log10(sort(p_real))), aes(x=resampled, y=real))+geom_point()+geom_abline()
}

figure("qq-plot of pvalues of corr.test of guide effects against null",
p_real_null_plot(function(x) psych::corr.test(x, Y, adjust = "none", ci=F)$p, guide_matrix)+
  xlab("quantile null (resampled guide labels, -log10)")+
  ylab("quantile p value of cor.test (-log10)")
)


