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

n = nrow(count_matrix)
p = ncol(count_matrix)

guide_matrix = guide_matrix[sample(n),]

covariates_dt = fread("results/covariates.dt.csv")
batch = model.matrix(~batch-1,data=covariates_dt[, .(batch=unique(batch)), keyby=i.cell_id])


Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix

n_folds_cells = 1
n_folds_genes = 5

if (!(n_folds_genes==n_folds_cells || n_folds_genes==1 || n_folds_cells==1)) stop("n_fold_cells and n_fold_genes must be equal or 1")
n_folds = max(n_folds_genes, n_folds_cells)

folds_cells = createFolds(seq_len(n), k = n_folds_cells, list = TRUE, returnTrain = FALSE)
folds_genes = createFolds(seq_len(p), k = n_folds_genes, list = TRUE, returnTrain = FALSE)

pb = txtProgressBar(min = 0, max = n_folds, initial = 0,  style = 3)

dt = rbindlist(lapply( seq_len(n_folds), function (fold) {

  test_cells = folds_cells[[min(fold, n_folds_cells)]]
  train_cells = if (n_folds_cells>1) seq_len(n)[-test_cells] else test_cells
  test_genes = folds_genes[[min(fold, n_folds_genes)]]
  train_genes = if (n_folds_genes>1) seq_len(p)[-test_genes] else test_genes

  cdr = rowMeans(Y[,train_genes]>0)
  mean_count = rowMeans(Y[,train_genes])
  capture = as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))

  adj_guide_matrix = function(X, X_covariates) {
    n_guides = ncol(X)
    X_new = copy(X)
    X = cbind(X, X_covariates)
    sigma = sqrt(sum(matrixStats::colVars(Y[train_cells, train_genes])))
    LL_same <- function(res)    rowSums(dnorm(res, sd=sigma, log = TRUE))

    fit = lm(Y[train_cells, train_genes] ~., as.data.table(X[train_cells,]))
    residuals_correct <- predict(fit, as.data.table(X)) - Y[,train_genes]

    LL_correct = LL_same(residuals_correct)

    beta = coef(fit)

    for (guide in seq_len(n_guides)) {
      guide_detected = X_new[, guide]
      residuals_swapped = sweep(residuals_correct[guide_detected,,drop=F],2, beta[guide,,drop=F])
      LL_swapped = LL_same(residuals_swapped)
      LLR = LL_correct[guide_detected]-LL_swapped
      X_new[guide_detected,guide] <- LLR>0
    }
    colnames(X_new) <- sprintf("%s.adj", colnames(X_new))
    X_new
  }



  res_SSq = function(X) {
    fit = lm(Y[train_cells, test_genes]~.-1, as.data.table(X[train_cells,]))
    residuals = Y[test_cells,test_genes]-predict(fit, as.data.table(X[test_cells,]))
    rowSums(residuals^2)
  }
  1
  #glmnet::glmnet(madness, Y[train_cells, test_genes], alpha=0.5, lambda = 0.0005, family = "mgaussian")
  dt = data.table(
    fold=fold,
    res_SSq = c(res_SSq(matrix(rep(1,nrow(guide_matrix),ncol=1))),
                res_SSq(cbind(1,guide_matrix)),
                res_SSq(cbind(1,capture)),
                res_SSq(cbind(batch)),
                res_SSq(cbind(capture,batch)),
                res_SSq(cbind(guide_matrix,capture,batch)),
                res_SSq(cbind(guide_matrix,capture,batch, adj_guide_matrix(guide_matrix, cbind(capture,batch)))),
                res_SSq(cbind(1,adj_guide_matrix(guide_matrix, cbind(capture,batch))))
                ),
    method = rep(c("intercept_only","+guides","+capture","+batch","+capture+batch","+guides+capture+batch","+guides+capture\n+adj.guides+batch","+adj.guides"), each=length(test_cells)),
    cell = test_cells
  )
  setTxtProgressBar(pb,fold)
  dt
}))

dt2=dt[method!="intercept_only"][dt[method=="intercept_only"], res_SSq_1 := i.res_SSq, on=.(cell)]

R2 = function(dt,idx) 1-dt[idx,sum(res_SSq)]/sum(dt[idx,sum(res_SSq_1)])
R2_dt = dt2[,  {boot.out= boot(.SD, R2, 1000); as.list(c(R2(.SD),boot.ci(boot.out,type="basic")$basic[4:5]))},by=method]
setnames(R2_dt,c("method","R^2","upper","lower"))
setorder(R2_dt,"R^2")
R2_dt[,method:=as.character(method)]

cross_val_info =
  if(n_folds_cells>1){
    if (n_folds_genes>1) "gene&cell crossvalidated" else "cell crossvalidated"
  } else {
    if (n_folds_genes>1) "genes crossvalidated"     else "no crossvalidation"
  }

figure(
  paste0("Variance explained through different adjustments - ", cross_val_info,"\nresampled guide matrix2"),
  ggplot(R2_dt, aes(x=factor(method),y=`R^2`)) + geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=lower, ymax=upper),width=2/3) +
    coord_flip() + ylab(expression(R[CV]^2)) + xlab("")+
    scale_x_discrete(limits=rev(R2_dt[, unique(method)]))
)

