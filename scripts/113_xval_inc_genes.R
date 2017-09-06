# -----------
library(data.table)
library(ggplot2)
library(glmnet)
library(mvtnorm)
library(caret)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

covariates_dt = fread("results/covariates.dt.csv")
batch = model.matrix(~batch-1,data=covariates_dt[, .(batch=unique(batch)), keyby=i.cell_id])


Y = log2(1+count_matrix)
X = cbind(guide_matrix,batch)

n = nrow(Y)
p = ncol(Y)

n_folds = 5
folds_cells = createFolds(seq_len(n), k = n_folds, list = TRUE, returnTrain = FALSE)

folds_genes = createFolds(seq_len(p), k = n_folds, list = TRUE, returnTrain = FALSE)

pb = txtProgressBar(min = 0, max = n_folds, initial = 0,  style = 3)

dt = rbindlist(lapply( seq_len(n_folds), function (fold) {
  isTest = folds_cells[[fold]]
  isTest_genes = folds_genes[[fold]]
  fit_with_guides = glmnet::glmnet(guide_matrix[-isTest,], Y[-isTest, isTest_genes], alpha=0.5, lambda = 0.0005, family = "mgaussian")
  fit_with_batch = glmnet::glmnet(X[-isTest,], Y[-isTest, isTest_genes], alpha=0.5, lambda = 0.0005, family = "mgaussian")
  fit_with_batch_lm = lm(Y[-isTest, isTest_genes]~., data.table(X[-isTest,]))

  cdr = rowMeans(Y[,-isTest_genes]>0)
  mean_count = rowMeans(Y[,-isTest_genes])
  cdr = as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))
  fit_with_batch_capture_efficiency = lm(Y[-isTest, isTest_genes]~.,data.table(cbind(X[-isTest,],cdr[-isTest,])))

  res_sum_of_squares_test = function(X) sum((Y[isTest,isTest_genes]-predict(lm()))^2)

  #glmnet::glmnet(madness, Y[-isTest, isTest_genes], alpha=0.5, lambda = 0.0005, family = "mgaussian")
  dt = data.table(
    fold=fold,
    res_ssq_intercept_only = res_sum_of_squares_test(colMeans(Y[-isTest,isTest_genes])),
    res_ssq_with_guides = res_sum_of_squares_test(predict(fit_with_guides, guide_matrix[isTest,])[,,1]),
    res_ssq_with_batch = res_sum_of_squares_test(predict(fit_with_batch, X[isTest,])[,,1]),
    res_ssq_with_batch_lm = res_sum_of_squares_test(predict(fit_with_batch_lm, data.table(X[isTest,]))),
    res_ssq_with_batch_capture_efficiency = res_sum_of_squares_test(predict(fit_with_batch_capture_efficiency, data.table(cbind(X[isTest,],cdr[isTest,]))))
  )
  setTxtProgressBar(pb,fold)
  dt
}))

melt(dt, id.vars=c("fold","res_ssq_intercept_only"))[, .(R2 = 1-sum(value)/sum(res_ssq_intercept_only)), by=.(variable)
   ][, c("upper","lower"):=as.list(1-(1-R2)/qf(c(0.025,0.975),df1=n,df2=n)),by=.(variable)
     ][]
