# -----------
library(data.table)
library(ggplot2)
library(glmnet)
library(mvtnorm)
library(caret)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

Y = log2(1+count_matrix)
X = guide_matrix

n = nrow(Y)
p = ncol(Y)

n_folds = 5
folds_cells = createFolds(seq_len(n), k = n_folds, list = TRUE, returnTrain = FALSE)

folds_genes = createFolds(seq_len(p), k = n_folds, list = TRUE, returnTrain = FALSE)

pb = txtProgressBar(min = 0, max = n_folds, initial = 0,  style = 3)

dt = rbindlist(lapply( seq_len(n_folds), function (fold) {
  isTest = folds_cells[[fold]]
  fit = glmnet::glmnet(X[-isTest,], Y[-isTest,], alpha=0.5, lambda = 0.0005, family = "mgaussian") #, standardize = FALSE
  dt = data.table(
    fold=fold,
    res_ssq_intercept_only = sum((sweep(Y[isTest,], 2, colMeans(Y[-isTest,])))^2),
    res_ssq_model = sum((Y[isTest,]-predict(fit, X[isTest,])[,,1])^2)
  )
  setTxtProgressBar(pb,fold)
  dt
}))


dt[, .(R2 = 1-sum(res_ssq_model)/sum(res_ssq_intercept_only))
   ][, c("upper","lower"):=as.list(1-(1-R2)/qf(c(0.025,0.975),df1=n,df2=n))
     ][]
