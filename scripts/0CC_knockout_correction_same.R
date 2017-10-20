# -----------
library(data.table)
library(ggplot2)
library(glmnet)
library(mvtnorm)
library(parallel)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

Y = log2(1+count_matrix)
X = guide_matrix

isTrain <- sample(nrow(X), ceiling(nrow(X)*2/3))
isNotTrain <- seq_len(nrow(X))[-isTrain]
isTestInNotTrain <- sample(rep(c(FALSE,TRUE), ceiling(length(isNotTrain)/2)), length(isNotTrain))

fit <- lm(Y[isTrain,] ~., as.data.table(X[isTrain,]))

residuals_train <- predict(fit, as.data.table(X[isTrain,])) - Y[isTrain,]

sigma = sqrt(sum(matrixStats::colVars(residuals_train)))


cl <- NULL
# Initiate cluster ---------
cl <- makeCluster(20)
clusterCall(cl, function() {
  library(glmnet);library(mvtnorm);library(parallel);library(data.table);library(stats)
})
# --------------

residuals_test_correct <- predict(fit, as.data.table(X[isNotTrain,])) - Y[isNotTrain,]
LL_same <- function(res)    rowSums(dnorm(res,    sd=sigma, log = TRUE))

LL_correct_same <- LL_same(residuals_test_correct)


if (!is.null(cl)) clusterExport(cl, ls())

ni <- seq_len(ncol(X)) # ni <- c(11,37,53,54)
pb = txtProgressBar(min = 0, max = length(ni), initial = 0,  style = 3)

workerfun <- function(guide) {
  X_  <-  copy(X[isNotTrain,])
  guide_detected <- X_[,guide]
  X_[,guide] <- !X_[,guide]
  residuals_test_swapped <- stats::predict(fit, as.data.table(X_)) - Y[isNotTrain,]

  LL_swapped_same  <- LL_same(residuals_test_swapped)

  dt1 <- data.table::data.table(
    guide = guide,
    guide_detected = guide_detected,
    LL_swapped = c(LL_swapped_same),
    LL_correct = c(LL_correct_same),
    LL_method = rep(c("same"), each=length(LL_swapped_same)),
    theta0 = mean(X[isTrain,guide]),
    isTest = isTestInNotTrain,
    guide_name=colnames(X)[guide]
  )
  setTxtProgressBar(pb,guide)
  dt1
}

dt <- rbindlist(if (!is.null(cl)) parLapply(cl,ni, workerfun) else lapply(ni, workerfun))


cross_entropy_lossf <-  function(p,label) -sum(log((!label)+(2*label-1)*p))/length(label)
cross_entropy_loss <-  function(ff,y, LLR) sum(log(1+exp(-(ff[1]*LLR+ff[2])*(2*y-1))))

dt[,LLR:=LL_correct-LL_swapped]
dt[(!guide_detected),LLR:=-LLR]

dt[, p_ko:=1/(1+exp(-(LLR+log(theta0)-log(1-theta0)))) ]

dt_bak <- dt
# dt_bak2 <- copy(dt_bak)
# dt_bak_overfit <- dt_bak
# --------------------------------------------------------------
dt <- copy(dt_bak)

dt <- dt[(isTest)]

res <- melt(dt,measure.vars = c("p_ko","p_ko_fitted"),value.name = "p_ko")[, method:=paste0(LL_method,variable)][(guide_detected) ,
   .(frac_confidently_perturbed = mean(p_ko>0.10)),by=.(guide,guide_name,method)]

ggplot(res[method=="samep_ko"], aes(y=frac_confidently_perturbed,x=method))+geom_violin()

dt[, .(x_entropy_loss=cross_entropy_lossf(p_ko,guide_detected)), by=.(LL_method)]


options(bitmapType='cairo')

ggplot(dt, aes(x= p_ko, linetype=guide_detected, color=LL_method, xintercept=theta0))+
  stat_ecdf() + facet_wrap("guide_name",scales="free_x") + geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")+
  ylab("ECDF")
ggsave(paste0("results/p_ecdf.png"), dpi = 400, width = 8, height=8)


ggplot(dt, aes(x= p_ko, linetype=guide_detected, color=LL_method, xintercept=theta0))+
  geom_density(aes(y=..scaled..),fill=NA, size=1,alpha=0.5) + facet_wrap("guide_name",scales="free_x") +
  geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")
ggsave(paste0("results/p_density.png"), dpi = 400, width = 8, height=8)



figure(
  paste0("split violin of LLR=LL_ko-LL_no_ko"),
  ggplot(dt, aes(color=guide_detected,y=LLR,x=set))+geom_split_violin()+
    facet_wrap("guide",scales ="free"),
  height=15,width=15
)


# ----------
stopCluster(cl)
