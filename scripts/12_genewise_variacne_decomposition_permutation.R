# load data -----------
# library(pertInv)
# library(data.table)
# library(cowplot)
# library(glmnet)
# library(mvtnorm)
# library(caret)
# library(boot)

library(pertInv)
data_set = "GSM2396858_k562_tfs_7"
# "GSM2396861_k562_ccycle"
# "GSM2396858_k562_tfs_7"
# "GSM2396859_k562_tfs_13"
# "GSM2396860_k562_tfs_highmoi"
# "GSM2396856_dc_3hr"
# "GSM2396857_dc_0hr"
data_folder <- paste0('data_processed/', data_set)

load(file = file.path(data_folder, "batch_matrix.RData"))
load(file = file.path(data_folder, "count_matrix.RData"))

load(file = file.path(data_folder, "guide_matrix.RData"))
covariates.dt <- fread(file.path(data_folder, "covariates.dt.csv"))

#count_matrix <- count_matrix[,1:100]#count_matrix#count_matrix[1:1000,1:100]
n_genes <- ncol(count_matrix)
n_cells <- nrow(count_matrix)
p <- n_genes
n <- n_cells
guide_matrix <- guide_matrix[1:n_cells,]
batch_matrix <- batch_matrix[1:n_cells,]

# transform --------------

Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix


wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))

# compute R2 cross validated ------------
n_folds_cells = 5
n_folds = n_folds_cells

folds_cells = caret::createFolds(seq_len(n), k = n_folds_cells, list = TRUE, returnTrain = FALSE)

pb = txtProgressBar(min = 0, max = n_folds, initial = 0,  style = 3)

dt = rbindlist(lapply( seq_len(n_folds), function (fold) {

  test_cells = folds_cells[[min(fold, n_folds_cells)]]
  train_cells = if (n_folds_cells>1) seq_len(n)[-test_cells] else test_cells

  cdr = rowMeans(count_matrix>0)
  mean_count = rowMeans(Y[,])
  capture = as.matrix(data.table(cdr,cdr^2,cdr^3,log(mean_count),log(mean_count)^2,log(mean_count)^3))


  res_SSq = function(method_name, X) {
    fit = lm(Y[train_cells, ]~.-1, as.data.table(X[train_cells,]))
    residuals = Y[test_cells,]-predict(fit, as.data.table(X[test_cells,]))
    data.table(method = method_name, res_SSq= colSums(residuals^2), gene = colnames(Y))
  }

  dt = rbind(res_SSq("intercept_only",matrix(rep(1,nrow(guide_matrix),ncol=1))),
             res_SSq("sgRNAs detected",cbind(1,guide_matrix)),
             res_SSq("sgRNAs (resampled)",cbind(1,guide_matrix[sample(nrow(guide_matrix)),]))
  )
  dt[, ':='(fold=fold)]
  setTxtProgressBar(pb,fold)
  dt
}))
# summarize results ----------
dt2=dt[, .(res_SSq=sum(res_SSq)), by=.(method, gene)]
dt2=dt2[method!="intercept_only"][dt2[method=="intercept_only"], res_SSq_0 := i.res_SSq, on=.(gene)]

dt2[, `R^2`:=1-res_SSq/res_SSq_0]
summary.dt <- dt2[, .(`mean R^2`=mean(`R^2`)), by=.(method)]
setorder(summary.dt, "mean R^2")[]


figure(
  paste0("Variance explained by guide (compared to resampled)"),
  ggplot(dt2, aes(x=factor(method),y=`R^2`)) +
    geom_hline(yintercept=0, linetype="dashed", color="gray")+
    geom_boxplot() +
    #stat_summary(fun.data = "mean_cl_boot", geom="errorbar",size=3, width=0.5, color="red") +
    #scale_shape_identity() +
    #geom_point(shape=124,size=5, alpha=0.5 ,color="black") +
    coord_flip() + ylab(expression(R[CV]^2)) + xlab("") +#+
    scale_x_discrete(limits=summary.dt[,method])
  , width=7, height=3
)

