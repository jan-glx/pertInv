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
n_genes <- ncol(count_matrix)
n_cells <- nrow(count_matrix)
load(file = file.path(data_folder, "guide_matrix.RData"))
covariates.dt <- fread(file.path(data_folder, "covariates.dt.csv"))





Y <- sweep(count_matrix, 2, colMeans(count_matrix), "/") # stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~0+I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch, data=covariates.dt) #+guide

lmfit <- lm.fit(x=cbind(X,guide_matrix), y=Y)
str(lmfit)

lmfit <- lm.fit(x=guide_matrix, y=Y)

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)



library(glmnet)
library(mvtnorm)
library(caret)



Y = log2(1+count_matrix)
X = cbind(guide_matrix,model.matrix(~0+I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch, data=covariates.dt)) #guide_matrix

n = nrow(Y)
p = ncol(Y)

n_folds = 5
folds_cells = createFolds(seq_len(n), k = n_folds, list = TRUE, returnTrain = FALSE)

folds_genes = createFolds(seq_len(p), k = n_folds, list = TRUE, returnTrain = FALSE)

pb = txtProgressBar(min = 0, max = n_folds, initial = 0,  style = 3)

dt = rbindlist(lapply( seq_len(n_folds), function (fold) {
  isTest = folds_cells[[fold]]
  fit = lm.fit(X[-isTest,], Y[-isTest,])#glmnet::glmnet, alpha=0.5, lambda = 0.0005, family = "mgaussian") #, standardize = FALSE
  dt = data.table(
    fold=fold,
    res_ssq_intercept_only = sum((sweep(Y[isTest,], 2, colMeans(Y[-isTest,])))^2),
    res_ssq_model = sum((Y[isTest,]-X[isTest,] %*%fit$coefficients)^2)
  )
  setTxtProgressBar(pb,fold)
  dt
}))


dt[, .(R2 = 1-sum(res_ssq_model)/sum(res_ssq_intercept_only))
   ][, c("upper","lower"):=as.list(1-(1-R2)/qf(c(0.025,0.975),df1=n,df2=n))][]











Y = log2(1+count_matrix)

plot(colMeans(Y), matrixStats::colVars(Y))
