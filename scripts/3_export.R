covariates.dt <- cells_summary.dt[cellnames.dt, on="cell_id"]

Y_tmp <- log10(count_matrix/rowSums(count_matrix)*1000+1)
fwrite(data.table(Y_tmp), "results/expression_logCp10k.csv")
fwrite(data.table(count_matrix), "results/expression_counts.csv")

Y <- stabilize_Anscombes(count_matrix) #count_matrix
fwrite(data.table(guide_matrix), "results/guide_matrix.csv")
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt)

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)
fwrite(data.table(Y_adj), "results/expression_adjusted.csv")

guide_matrix <- cbc_gbc_dict.dt[cells_summary.dt[cellnames.dt, on="cell_id"], on="cell"
                                ][!is.na(guide_id), as.matrix(Matrix::sparseMatrix(cell_id2, guide_id, x=TRUE))]
fwrite(data.table(guide_matrix), "results/guide_matrix.csv")
