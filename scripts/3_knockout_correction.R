# -----------
library(data.table)
library(ggplot2)
library(glmnet)
library(mvtnorm)
library(parallel)
library(xgboost)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

Y = log2(1+count_matrix)
X = guide_matrix

isTrain <- sample(nrow(X), ceiling(nrow(X)*2/3))
isNotTrain <- seq_len(nrow(X))[-isTrain]
isTestInNotTrain <- sample(rep(c(FALSE,TRUE), ceiling(length(isNotTrain)/2)), length(isNotTrain))

fit <- glmnet::glmnet(X[isTrain,], Y[isTrain,], alpha=0.5, lambda = 0.0005, family = "mgaussian", standardize = FALSE)

residuals_train <- predict(fit, X[isTrain,])[,,1] - Y[isTrain,]

sigma <- cov(residuals_train)
sigma_res <- sqrt(diag(sigma))
sigma_genes <- apply(Y[isNotTrain,], MARGIN = 2, sd)

theta0 <- colSums(X[isNotTrain,])/length(isNotTrain)

cl <- NULL
# Initiate cluster ---------
cl <- makeCluster(20)
clusterCall(cl, function() {
  library(glmnet);library(mvtnorm);library(parallel);library(data.table);library(stats);library(xgboost)})
# --------------

residuals_test_correct <- predict(fit, X[isNotTrain,])[,,1] - Y[isNotTrain,]
LL_mnorm <- function(res) mvtnorm::dmvnorm(res, sigma=sigma,                    log=TRUE)
LL_norm  <- function(res)    rowSums(dnorm(res,    sd=sigma_res,                log = TRUE))
LL_dixit <- function(res)    rowSums(dnorm(res,    sd=sqrt(sigma_genes^2+0.25), log = TRUE))
LL_same <- function(res)    rowSums(dnorm(res,    sd=sqrt(sum(sigma_res^2)), log = TRUE))


LL_correct_mnorm <- LL_mnorm(residuals_test_correct)
LL_correct_norm  <- LL_norm(residuals_test_correct)
LL_correct_dixit <- LL_dixit(residuals_test_correct)
LL_correct_same <- LL_same(residuals_test_correct)


if (!is.null(cl)) clusterExport(cl, ls())

ni <- seq_len(ncol(X)) # ni <- c(11,37,53,54)

workerfun <- function(i) {
  tryCatch({
    X  <-  guide_matrix[isNotTrain,]
    guide_detected <- X[,i]
    X[,i] <- !X[,i]
    residuals_test_swapped <- stats::predict(fit, X)[,,1] - Y[isNotTrain,]

    LL_swapped_mnorm <- LL_mnorm(residuals_test_swapped)
    LL_swapped_norm  <- LL_norm( residuals_test_swapped)
    LL_swapped_dixit <- LL_dixit(residuals_test_swapped)
    LL_swapped_same  <- LL_same(residuals_test_swapped)

    dt1 <- data.table::data.table(
      guide = i,
      guide_detected = guide_detected,
      LL_swapped = c(LL_swapped_mnorm, LL_swapped_norm, LL_swapped_dixit, LL_swapped_same),
      LL_correct = c(LL_correct_mnorm, LL_correct_norm, LL_correct_dixit, LL_correct_same),
      LL_method = rep(c("mvnorm", "norm", "dixit", "same"), each=length(LL_swapped_mnorm)),
      theta0 = theta0[i],
      isTest = isTestInNotTrain
    )

    dt1[(guide_detected), LL_perturbed:=LL_correct]
    dt1[(guide_detected), LL_not_perturbed:=LL_swapped]
    dt1[(!guide_detected), LL_perturbed:=LL_swapped]
    dt1[(!guide_detected), LL_not_perturbed:=LL_correct]
    dt1[, p_ko:=1/(1+exp(-(LL_perturbed-LL_not_perturbed+log(theta0)-log(1-theta0))))]


    dtrain <- xgb.DMatrix(Y[isTrain,], label =  guide_matrix[isTrain, i])
    #xgb.cv(list(objective="binary:logistic", eta=0.02,eval_metric="logloss",lambda=1,alpha=1), dtrain, nfold=3, nrounds=100)
    #fitted_booster <- xgb.train(list(objective="binary:logistic", eta=0.02,eval_metric="logloss",lambda=1,alpha=1), dtrain, nrounds=100)
    fitted_booster <- xgb.train(list(booster="gbtree",objective="binary:logistic", eta=0.02,eval_metric="logloss",
                                     subsample=0.8,colsample_bytree=0.8,max_depth=1,eval_metric="logloss"), dtrain, nrounds=300)
    dt2 <- data.table(guide=i,
                      p_ko=predict(fitted_booster, Y[isNotTrain,]),
                      guide_detected=guide_matrix[isNotTrain, i],
                      LL_method="xgboost",
                      isTest = isTestInNotTrain)
    print(".")
    merge(dt1, dt2, all=TRUE)
  },
  error=function(e) e
  )
}

dt <- rbindlist(if (!is.null(cl)) parLapply(cl,ni, workerfun) else lapply(ni, workerfun))

dt[,guide_name:=colnames(guide_matrix)[guide]]


cross_entropy_lossf <-  function(p,label) -sum(log((!label)+(2*label-1)*p))/length(label)

cross_entropy_loss <-  function(ff,y, LLR) sum(log(1+exp(-(ff[1]*LLR+ff[2])*(2*y-1))))

dt[,LLR:=LL_perturbed-LL_not_perturbed]
dt[LL_method=="xgboost", LLR:=-log(1/p_ko-1)]
dt[LL_method!="xgboost", p_ko:=1/(1+exp(-(LLR+log(theta0)-log(1-theta0)))) ]



dt[, c("scaling_factor","offset"):=
     as.list(optim(par = c(1,0.001), fn = cross_entropy_loss, y=guide_detected[!isTest], LLR=LLR[!isTest])$par),
  by= .(LL_method,guide)]

dt[, p_ko_fitted := 1/(1+exp(-(scaling_factor*LLR+offset))), by= .(LL_method,guide)]

dt


dt_bak <- dt
# dt_bak2 <- copy(dt_bak)
# dt_bak_overfit <- dt_bak
# --------------------------------------------------------------
dt <- copy(dt_bak)

dt <- dt[(isTest)]

res <- melt(dt,measure.vars = c("p_ko","p_ko_fitted"),value.name = "p_ko")[, method:=paste0(LL_method,variable)][(guide_detected) ,
   .(frac_confidently_perturbed = mean(p_ko>0.10)),by=.(guide,guide_name,method)]

ggplot(res[method=="samep_ko_fitted"], aes(y=frac_confidently_perturbed,x=method))+geom_violin()


res <- dt[, .(x_entropy_loss=cross_entropy_lossf(p_ko,guide_detected),
              x_entropy_loss_fit=cross_entropy_lossf(p_ko_fitted,guide_detected)),
          by=.(LL_method,guide,guide_name)]
setorder(res,"x_entropy_loss_fit")
res
res[,r_fit:= rank(x_entropy_loss_fit),by=guide]
res[,r:= rank(x_entropy_loss),by=guide]
res[,sum(r_fit),by=LL_method]
res[,sum(r),by=LL_method]


dt <- dt[guide>50 | guide==36]

options(bitmapType='cairo')

ggplot(dt, aes(x= p_ko_fitted, linetype=guide_detected, color=LL_method, xintercept=theta0))+
  stat_ecdf() + facet_wrap("guide_name",scales="free_x") + geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")+
  ylab("ECDF")
ggsave(paste0("results/p_fitted_ecdf.png"), dpi = 400, width = 8, height=8)

ggplot(dt, aes(x= p_ko, linetype=guide_detected, color=LL_method, xintercept=theta0))+
  stat_ecdf() + facet_wrap("guide_name",scales="free_x") + geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")+
  ylab("ECDF")
ggsave(paste0("results/p_ecdf.png"), dpi = 400, width = 8, height=8)



ggplot(dt, aes(x= p_ko_fitted, linetype=guide_detected, color=LL_method, xintercept=theta0))+
  geom_density(aes(y=..scaled..),fill=NA, size=1,alpha=0.5) + facet_wrap("guide_name",scales="free_x") +
  geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")
ggsave(paste0("results/p_fitted_density.png"), dpi = 400, width = 8, height=8)

ggplot(dt, aes(x= p_ko, linetype=guide_detected, color=LL_method, xintercept=theta0))+
  geom_density(aes(y=..scaled..),fill=NA, size=1,alpha=0.5) + facet_wrap("guide_name",scales="free_x") +
  geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")
ggsave(paste0("results/p_density.png"), dpi = 400, width = 8, height=8)





for (method in c("mvnorm", "norm", "dixit","same","xgboost")) {
  ggplot(dt[LL_method == method], aes(x= p_ko_fitted, color=guide_detected, linetype=LL_method))+
    geom_density() + facet_wrap("guide_name", scales = "free")
  ggsave(paste0("results/p_fitted_density_",method,".png"), dpi = 400, width = 12, height=8)

  ggplot(dt[LL_method == method], aes(x= p_ko, color=guide_detected, linetype=LL_method))+
    geom_density() + facet_wrap("guide_name", scales = "free")
  ggsave(paste0("results/p_density_",method,".png"), dpi = 400, width = 12, height=8)

  ggplot(dt[LL_method == method],
         aes(x= LLR, color=guide_detected, linetype=LL_method))+
    geom_density() + facet_wrap("guide_name", scales = "free")
  ggsave(paste0("results/LLR_density_",method,".png"), dpi = 400, width = 12, height=8)

  ggplot(dt[LL_method == method][guide_name =="m_Stat1_3"][guide_detected==TRUE],
         aes(x=LLR, linetype=LL_method)) +
    geom_histogram(bins=20) + ggtitle("m_Stat1_3")
  ggsave(paste0("results/LLR_density_m_Stat1_3_",method,".png"), dpi = 400, width = 8, height=8)
}
# ----------
stopCluster(cl)
