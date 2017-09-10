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

#guide_matrix = guide_matrix[sample(n),]

covariates_dt = fread("results/covariates.dt.csv")
batch = model.matrix(~batch-1,data=covariates_dt[, .(batch=unique(batch)), keyby=i.cell_id])

capture = local({
  cdr = rowMeans(Y[,]>0)
  mean_count = rowMeans(Y[,])
  as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))
})

Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix
n_guides = ncol(X)


LLR_knockout = function(Y, X, X_covariates, train_cells, test_cells, pb) {
  n_guides = ncol(X)

  sigma = sqrt(sum(matrixStats::colVars(Y[train_cells, ])))
  LL_same <- function(res) rowSums(dnorm(res, sd=sigma, log = TRUE))

  X_ = as.data.table(cbind(X, X_covariates))
  fit = lm(Y[train_cells, ] ~.-1, X_[train_cells,])
  beta = coef(fit)

  residuals_correct <- predict(fit, X_[test_cells,]) - Y[test_cells,]
  LL_correct = LL_same(residuals_correct)
  pb$tick()

  rbindlist(lapply(seq_len(n_guides), function(guide) {
    theta0 = mean(X[train_cells, guide])
    guide_detected = X[test_cells, guide]
    residuals_swapped = residuals_correct - outer(2*guide_detected-1, beta[guide,])
    LL_swapped = LL_same(residuals_swapped)
    LLR = LL_correct-LL_swapped
    LLR[!guide_detected] = -LLR[!guide_detected]
    pb$tick()
    data.table(guide,LLR,theta0,guide_detected, cell=test_cells)
  }))
}

helper = function (n_folds_cells, pb) {
  folds_cells = createFolds(seq_len(n), k = n_folds_cells, list = TRUE, returnTrain = FALSE)

  rbindlist(lapply(seq_len(n_folds_cells), function (fold) {
    test_cells = folds_cells[[min(fold, n_folds_cells)]]
    train_cells = if (n_folds_cells>1) seq_len(n)[-test_cells] else test_cells
    LLR_knockout(Y, X, cbind(capture, batch), train_cells, test_cells, pb)[, fold:=fold]
  }))
}


pb = progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                total = 6*(n_guides+1),
                                clear = FALSE, width= 60)
dt = rbind(helper(5,pb)[,set:="test"], helper(1,pb)[,set:="train"])

# -----------
ggplot(dt, aes(x=LLR,  linetype=set, color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
  facet_wrap("guide",scales ="free")



figure(
  paste0("distribution of LLR=LL_ko-LL_no_ko"),#\n(resampled)
  ggplot(dt, aes(fill=guide_detected,y=LLR,x=set))+geom_split_violin()+
    facet_wrap("guide",scales ="free"),
  height=15,width=15
)

dt[, cross_entropy_loss(exp(LLR+log(theta0)-log(1-theta0)),guide_detected),by=.(set)]
dt[is.na(LLR), LLR:=0]
dt[,theta0_:=(theta0*.N+dt[,mean(guide_detected)])/(.N+1),by=.(fold,set,guide)]
dt[,cross_entropy_loss(exp(LLR+log(theta0_)-log(1-theta0_)),guide_detected),by=.(set)]

res = dt[set=="test", knnDemix::mixture.test(LLR[!guide_detected],LLR[guide_detected]),by=.(set,guide)]
res = res[,.(lower=min(conf.int),upper=max(conf.int)),by=setdiff(colnames(res),"conf.int")]
res[, mean(p.adjust(p.value,method="BH")<0.05, na.rm = TRUE)]
qqplot(-log10(seq(0,1,length.out = 60)[-c(1,60)]),res[, -log10(p.value)])
abline(0,1)



ggplot(dt[set=="test"][guide %in% res[estimate<0.9,guide]], aes(x=LLR,color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
  facet_wrap("guide",scales ="free")
