# -----------
library(pertInv)
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



Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix[, !colnames(guide_matrix) %in% c("m_Egr1_3", "m_MouseNTC_100_A_67005")] # this is to make the intercept meaningful, otherwise the guide matrix is almost colinear (schould remove rows without detected guides though)

capture = local({
  cdr = rowMeans(Y[,]>0)
  mean_count = rowMeans(Y[,])
  as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))
})


guides = colnames(X)
n_guides = length(guides)

LLR_knockout = function(Y, X_, train_cells, test_cells, pb) {
  sigma = sqrt(sum(matrixStats::colVars(Y[train_cells, ])))
  LL_same <- function(res) rowSums(dnorm(res, sd=sigma, log = TRUE))

  fit = lm(Y[train_cells, ] ~.-1, X_[train_cells,])
  beta = coef(fit)

  residuals_correct <- predict(fit, X_[test_cells,]) - Y[test_cells,]
  LL_correct = LL_same(residuals_correct)
  pb$tick()

  rbindlist(lapply(guides, function(guide) {
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

X_ = as.data.table(cbind(X, capture, batch))

helper = function (n_folds_cells, pb) {
  folds_cells = createFolds(seq_len(n), k = n_folds_cells, list = TRUE, returnTrain = FALSE)

  rbindlist(lapply(seq_len(n_folds_cells), function (fold) {
    test_cells = folds_cells[[min(fold, n_folds_cells)]]
    train_cells = if (n_folds_cells>1) seq_len(n)[-test_cells] else test_cells
    LLR_knockout(Y, X_, train_cells, test_cells, pb)[, fold:=fold]
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

knnDemix::mixture.test(runif(10000), res[,p.value]) # at least nguides*(1-upper(alpha)) guides work

setorder(res,p.value)

out =rbindlist(lapply(1:10, function(i) copy(res)[,rep:=i][,runif:=sort(runif(.N))]))

figure(
  "qq-plot of heterogenity-test",
  ggplot(out, aes(x=-log10(p.value), y= -log10(runif)))+stat_summary(fun.data=mean_sd)+
  geom_abline()+coord_flip()
  )

# figure(
#   "BH-plot of heterogenity-test",
#   ggplot(res, aes(y=-log10(p.value), x= rank(p.value))+
#     geom_abline(slope = 0.05/nrow(res))+coord_flip()
# )
# ggplot(res, aes(y=p.value, x= rank(-p.value)))+geom_point()+
#   geom_abline(slope = 0.05/nrow(res))+stat_function(fun = function(x) 1-(1-0.025)^(1/x))+stat_function(fun = function(x) 1-(1-0.975)^(1/x))+stat_function(fun = function(x)1-(1-0.5)^(1/x))+scale_y_sqrt()


figure(
  "LLR distribution of efficient guides",
  ggplot(dt[set=="test"][guide %in% res[p.adjust(p.value,method="BH")<0.05,guide]], aes(x=LLR,color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
    facet_wrap("guide",scales ="free")+geom_text(data=res[p.adjust(p.value,method="BH")<0.05,],aes(label=sprintf(" %s",signif(estimate,3)),x=-Inf,y=Inf),vjust=1,hjust=0,color="black"),
  width=12, height=12
)
