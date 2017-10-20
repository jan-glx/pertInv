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

#count_matrix = count_matrix[,1:100]
#guide_matrix = guide_matrix[,1:5]
guide_matrix = guide_matrix[, !colnames(guide_matrix) %in% c("m_Egr1_3", "m_MouseNTC_100_A_67005")] # this is to make the intercept meaningful, otherwise the guide matrix is almost colinear (schould remove rows without detected guides though)


n = nrow(count_matrix)
p = ncol(count_matrix)
n_guides = ncol(guide_matrix)

#guide_matrix = guide_matrix[sample(n),]

covariates_dt = fread("results/covariates.dt.csv")
batch = model.matrix(~batch-1,data=covariates_dt[, .(batch=unique(batch)), keyby=i.cell_id])


Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix

n_folds_cells = 3
n_folds_genes = 3

if (!(n_folds_genes==n_folds_cells || n_folds_genes==1 || n_folds_cells==1)) stop("n_fold_cells and n_fold_genes must be equal or 1")
n_folds = max(n_folds_genes, n_folds_cells)

folds_cells = createFolds(seq_len(n), k = n_folds_cells, list = TRUE, returnTrain = FALSE)
folds_genes = createFolds(seq_len(p), k = n_folds_genes, list = TRUE, returnTrain = FALSE)

pb = progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                total =  n_folds*(n_guides+1),
                                clear = FALSE, width= 60)
dt = rbindlist(lapply( seq_len(n_folds), function (fold) {

  test_cells = folds_cells[[min(fold, n_folds_cells)]]
  train_cells = if (n_folds_cells>1) seq_len(n)[-test_cells] else test_cells
  calib_cells = sample(train_cells, ceiling(length(train_cells)/5))
  train_cells = setdiff(train_cells,calib_cells)
  test_genes = folds_genes[[min(fold, n_folds_genes)]]
  train_genes = if (n_folds_genes>1) seq_len(p)[-test_genes] else test_genes

  cdr = rowMeans(Y[,train_genes]>0)
  mean_count = rowMeans(Y[,train_genes])
  capture = as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))

  LL <- function(res, sigma) apply(res, MARGIN=1, function(x) sum(dnorm(x, sd=sigma,log = TRUE)))

  adj_guide_matrix = function(X, X_covariates) {
    n_guides = ncol(X)
    X_adj = copy(X)
    X_ = as.data.table(cbind(X, X_covariates))
    sigma_dixit = apply(Y[train_cells, train_genes], MARGIN = 2, function(x) sqrt(sd(x)^2+1/4))

    fit = lm(Y[train_cells, train_genes] ~.-1, X_[train_cells,])

    residuals_correct <- predict(fit, X_) - Y[,train_genes]

    LL_correct = LL(residuals_correct, sigma_dixit)

    beta = coef(fit)


    theta0 = colMeans(X[train_cells,])
    LLR0 = log(theta0)-log(1-theta0)

    LLRs <- copy(X)
    for (guide in seq_len(n_guides)) {
      guide_detected = X[, guide]
      residuals_swapped = sweep(residuals_correct,2, beta[guide,,drop=F])
      LL_swapped = LL(residuals_swapped, sigma_dixit)
      LLRs[,guide] = LL_correct-LL_swapped
      LLRs[!guide_detected,guide] = -LLRs[!guide_detected,guide]
      pb$tick()
      # LLR = log(P(expression|guide detected) / P(expression|guide not detected))
    }


    cross_entropy_loss_ = function(scaling_factor) {
      LLR_tot = sweep(scaling_factor*LLRs[calib_cells,],2,LLR0,`+`)
      LLR_tot[!X[calib_cells,]] = -LLR_tot[!X[calib_cells,]]
      mean(log(1+exp(-LLR_tot)))
    }

    dt = data.table(melt(LLRs),melt(X)[,3], set=c("test","train","calib")[1+seq_len(n)%in%train_cells+2*seq_len(n)%in%calib_cells])
    setnames(dt,c("cell","guide","LLR","guide_detected","set"))
    if(fold==1) figure("theta density",
      ggplot(dt,aes(x=1/(1+exp(-LLR)),group=interaction(guide,guide_detected),color=guide_detected))+
        geom_density()+facet_grid(.~set)
    )
    fit_res <- optimize(f = cross_entropy_loss_, interval=c(-3,3))
    scaling_factor <- fit_res$minimum

    if(fold==1) figure("theta density calibrated",
      ggplot(dt,aes(x=1/(1+exp(-LLR*scaling_factor)), group=interaction(guide,guide_detected),color=guide_detected))+
        geom_density()+facet_grid(.~set)
    )
    X_adj =  1/(1+exp(-LLRs))
    X_adj[!X] = 0
    X_adj_calib = 1/(1+exp(-LLRs*scaling_factor))
    X_adj_calib[!X] = 0
    colnames(X_adj) <- sprintf("%s.adj", colnames(X_adj))
    colnames(X_adj_calib) <- sprintf("%s.adj.calib", colnames(X_adj_calib))
    list(X_adj, X_adj_calib)
  }



  res_SSq = function(X) {
    fit = lm(Y[train_cells, test_genes]~.-1, as.data.table(X[train_cells,]))
    residuals = Y[test_cells,test_genes]-predict(fit, as.data.table(X[test_cells,]))
    rowSums(residuals^2)
  }
  adj = adj_guide_matrix(X, cbind(capture,batch))
  X_adj = adj[[1]]
  X_adj_calib = adj[[2]]
  X_adj_thresh = X_adj>0.5
  X_adj2 = copy(X_adj)
  X_adj2[!X]=NA
  medians= apply(X_adj2, 2, median, na.rm=TRUE)

  X_adj_thresh2 = sweep(X_adj, 2, medians, `>`)


  dt = data.table(
    fold=fold,
    res_SSq = c(res_SSq(matrix(rep(1,n,ncol=1))),
                res_SSq(cbind(capture,batch)),
                res_SSq(cbind(X,capture,batch)),
                res_SSq(cbind(X_adj,capture,batch)),
                res_SSq(cbind(X_adj_calib,capture,batch)),
                res_SSq(cbind(X_adj_thresh,capture,batch)),
                res_SSq(cbind(X_adj_thresh2,capture,batch)),
                res_SSq(cbind(X,X_adj_thresh,capture,batch)),
                res_SSq(cbind(X,X_adj_thresh2,capture,batch))
                ),
    method = rep(c("intercept_only",
                   "+capture+batch",
                   "+guides+capture+batch",
                   "+adj.guides\n+capture+batch",
                   "+adj_calib.guides\n+capture+batch",
                   "+adj_thresh.guides\n+capture+batch",
                   "+adj_thresh2.guides\n+capture+batch",
                   "+guides+adj_thresh.guides\n+capture+batch",
                   "+guides+adj_thresh2.guides\n+capture+batch"
                   ), each=length(test_cells)),
    cell = test_cells
  )
  pb$tick()
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
    if (n_folds_genes>1) "gene crossvalidated"      else "no crossvalidation"
  }

figure(
  paste0("Variance explained through different adjustments - ", cross_val_info),
  ggplot(R2_dt, aes(x=factor(method),y=`R^2`)) + geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=lower, ymax=upper),width=2/3) +
    coord_flip() + ylab(expression(R[CV]^2)) + xlab("")+
    scale_x_discrete(limits=rev(R2_dt[, unique(method)])),
  height=6
)

