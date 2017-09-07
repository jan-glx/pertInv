# -----------
library(pertInv)
library(data.table)
library(cowplot)
library(glmnet)
library(mvtnorm)
library(caret)
library(boot)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

covariates_dt = fread("results/covariates.dt.csv")
batch = model.matrix(~batch-1,data=covariates_dt[, .(batch=unique(batch)), keyby=i.cell_id])


Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix

n = nrow(Y)
p = ncol(Y)

n_folds = 5
folds_cells = createFolds(seq_len(n), k = n_folds, list = TRUE, returnTrain = FALSE)

folds_genes = createFolds(seq_len(p), k = n_folds, list = TRUE, returnTrain = FALSE)

pb = txtProgressBar(min = 0, max = n_folds, initial = 0,  style = 3)

dt = rbindlist(lapply( seq_len(n_folds), function (fold) {
  isTest = folds_cells[[fold]]
  isTest_genes = folds_genes[[fold]]

  cdr = rowMeans(Y[,-isTest_genes]>0)
  mean_count = rowMeans(Y[,-isTest_genes])
  capture = as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))

  res_SSq = function(X) {
    fit =lm(Y[-isTest, isTest_genes]~.-1, as.data.table(X[-isTest,]))
    residuals = Y[isTest,isTest_genes]-predict(fit, as.data.table(X[isTest,]))
    rowSums(residuals^2)
  }

  #glmnet::glmnet(madness, Y[-isTest, isTest_genes], alpha=0.5, lambda = 0.0005, family = "mgaussian")
  dt = data.table(
    fold=fold,
    res_SSq = c(res_SSq(matrix(rep(1,nrow(guide_matrix),ncol=1))),
                res_SSq(cbind(1,guide_matrix)),
                res_SSq(cbind(1,capture)),
                res_SSq(cbind(batch)),
                res_SSq(cbind(capture,batch)),
                res_SSq(cbind(guide_matrix,capture,batch))
                ),
    method = rep(c("intercept_only","+guides","+capture","+batch","+capture+batch","+guides+capture+batch"), each=length(isTest)),
    cell = isTest
  )
  setTxtProgressBar(pb,fold)
  dt
}))

dt2=dt[method!="intercept_only"][dt[method=="intercept_only"], res_SSq_1 := i.res_SSq, on=.(cell)]

R2 = function(dt,idx) 1-dt[idx,sum(res_SSq)]/sum(dt[idx,sum(res_SSq_1)])
R2_dt = dt2[,  {boot.out= boot(.SD, R2, 1000); as.list(c(R2(.SD),boot.ci(boot.out,type="basic")$basic[4:5]))},by=method]
setnames(R2_dt,c("method","R^2","upper","lower"))

figure(
  "Variance explained through different adjustments - log2(1+count)",
  ggplot(R2_dt, aes(x=method,y=`R^2`)) + geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=lower, ymax=upper),width=2/3) +
    coord_flip() + ylab(expression(R[CV]^2)) + xlab("")
)

