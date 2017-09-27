# -----------
library(pertInv)
library(ggplot2)
library(glmnet)
library(mvtnorm)
library(parallel)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")
guide_matrix = guide_matrix[, !colnames(guide_matrix) %in% c("m_Egr1_3", "m_MouseNTC_100_A_67005")] # this is to make the intercept meaningful, otherwise the guide matrix is almost colinear (schould remove rows without detected guides though)


# count_matrix = count_matrix[1:1000,1:100]
# guide_matrix = guide_matrix[1:1000,1:10]
# -------------
n = nrow(count_matrix)
p = ncol(count_matrix)
n_guides = ncol(guide_matrix)

Y = log2(1+count_matrix)
X = guide_matrix

isTrain <- sample(nrow(X), ceiling(nrow(X)*2/3))
isNotTrain <- seq_len(nrow(X))[-isTrain]
isTestInNotTrain <- sample(rep(c(FALSE,TRUE), ceiling(length(isNotTrain)/2)), length(isNotTrain))

fit <- lm(Y[isTrain,] ~., as.data.table(X[isTrain,]))

residuals_train <- predict(fit, as.data.table(X[isTrain,])) - Y[isTrain,]

Sigma <- cov(residuals_train)
sigma_res <- sqrt(diag(Sigma))
sigma_genes <- apply(Y[isTrain,], MARGIN = 2, sd)
sigma_poisson <- sqrt(1/apply(2^Y[isTrain,]-1, MARGIN = 2, mean)+0.02138951)

# devtools::install_github("dgrtwo/ebbr")
theta0 <- ebbr::add_ebb_estimate(data.table(k = colSums(X[isNotTrain,]), n = length(isNotTrain)), k , n)[,.fitted]

cl <- NULL # cl <- makeCluster(20)
# Initiate cluster ---------

if (!is.null(cl)) clusterCall(cl, function() {
  library(glmnet);library(mvtnorm);library(parallel);library(data.table);library(stats)
})
# --------------

residuals_test_correct <- predict(fit, as.data.table(X[isNotTrain,])) - Y[isNotTrain,]


LL <- function(res, sigma) apply(res, MARGIN=1, function(x) sum(dnorm(x,  sd=sigma,log = TRUE)))

LL_norm  <- function(res, y) LL(res, sigma_res)
LL_dixit <- function(res, y) LL(res, sqrt(sigma_genes^2+0.25))
LL_same <- function(res, y)  LL(res, sum(sigma_res))
LL_poisson <- function(res, y) LL(res, sigma_poisson)

LL_mnorm <- function(res, y) mvtnorm::dmvnorm(res, sigma=Sigma, log=TRUE)


LL_correct_mnorm <- LL_mnorm(residuals_test_correct)
LL_correct_norm  <- LL_norm(residuals_test_correct)
LL_correct_dixit <- LL_dixit(residuals_test_correct)
LL_correct_same <- LL_same(residuals_test_correct)
LL_correct_poisson <- LL_poisson(residuals_test_correct)


if (!is.null(cl)) clusterExport(cl, ls())

ni <- seq_len(ncol(X))
pb = progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                total =  length(ni),
                                clear = FALSE, width= 60)

workerfun <- function(guide) {
  X_  <-  guide_matrix[isNotTrain,]
  guide_detected <- X_[,guide]
  X_[,guide] <- !X_[,guide]
  residuals_test_swapped <- stats::predict(fit, as.data.table(X_)) - Y[isNotTrain,]

  LL_swapped_mnorm <- LL_mnorm(residuals_test_swapped)
  LL_swapped_norm  <- LL_norm( residuals_test_swapped)
  LL_swapped_dixit <- LL_dixit(residuals_test_swapped)
  LL_swapped_same  <- LL_same(residuals_test_swapped)
  LL_swapped_poisson  <- LL_poisson(residuals_test_swapped)

  dt1 <- data.table::data.table(
    guide_id = guide,
    guide = colnames(guide_matrix)[guide],
    guide_detected =  rep(guide_detected, 5),
    LL_swapped = c(LL_swapped_mnorm, LL_swapped_norm, LL_swapped_dixit, LL_swapped_same, LL_swapped_poisson),
    LL_correct = c(LL_correct_mnorm, LL_correct_norm, LL_correct_dixit, LL_correct_same, LL_correct_poisson),
    LL_method = rep(c("mvnorm", "norm", "dixit", "same", "poisson"), each=length(LL_swapped_mnorm)),
    theta0 = theta0[guide],
    isTest = isTestInNotTrain
  )

  dt1[(guide_detected), LL_perturbed:=LL_correct]
  dt1[(guide_detected), LL_not_perturbed:=LL_swapped]
  dt1[(!guide_detected), LL_perturbed:=LL_swapped]
  dt1[(!guide_detected), LL_not_perturbed:=LL_correct]
  pb$tick()
  dt1
}

dt <- rbindlist(if (!is.null(cl)) parLapply(cl,ni, workerfun) else lapply(ni, workerfun))

dt[, LLR := LL_perturbed-LL_not_perturbed]
dt[, LLR0 := log(theta0)-log(1-theta0)]

dt[, p_ko:= 1/(1+exp(-(LLR+LLR0)))]



dt[,scaling_factor := {
  cross_entropy_loss_ = function(scaling_factor) {
    LLR_tot = -(scaling_factor*LLR+LLR0)
    LLR_tot[guide_detected] = -LLR_tot[guide_detected]
    mean(log(1+exp(-LLR_tot[!isTest])))
  }
  fit_res <- optimize(f = cross_entropy_loss_, interval=c(-3,3))
  fit_res$minimum
  },
  by=LL_method
  ]


dt[, p_ko_fitted:= 1/(1+exp(-(scaling_factor*LLR+LLR0)))]

dt


dt_bak <- dt
# dt_bak2 <- copy(dt_bak)
# dt_bak_overfit <- dt_bak
# --------------------------------------------------------------
dt <- copy(dt_bak)

dt <- dt[(isTest)]

res <- melt(dt,measure.vars = c("p_ko","p_ko_fitted"),value.name = "p_ko")[, method:=paste0(LL_method,variable)][(guide_detected) ,
   .(frac_confidently_perturbed = mean(p_ko>0.10)),by=.(guide,method)]

ggplot(res[method=="dixitp_ko_fitted"], aes(y=frac_confidently_perturbed,x=method))+geom_violin()


res <- dt[, .(x_entropy_loss=cross_entropy_loss(p_ko,guide_detected),
              x_entropy_loss_fit=cross_entropy_loss(p_ko_fitted,guide_detected)),
          by=.(LL_method,guide)]
setorder(res,"x_entropy_loss_fit")
res
res[,r_fit:= rank(x_entropy_loss_fit),by=guide]
res[,r:= rank(x_entropy_loss),by=guide]
res[,sum(r_fit),by=LL_method]
res[,sum(r),by=LL_method]

dt[, .(x_entropy_loss=cross_entropy_loss(p_ko,guide_detected),
       x_entropy_loss_fit=cross_entropy_loss(p_ko_fitted,guide_detected)),
   by=.(LL_method)]


options(bitmapType='cairo')

figure("p_fitted_ecdf",
       ggplot(dt, aes(x= p_ko_fitted, linetype=guide_detected, color=LL_method, xintercept=theta0))+
         stat_ecdf() + facet_wrap("guide",scales="free_x") + geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")+
         ylab("ECDF"),
       dpi = 400, width = 16, height=16
)

figure("p_ecdf",
       ggplot(dt, aes(x= p_ko, linetype=guide_detected, color=LL_method, xintercept=theta0))+
         stat_ecdf() + facet_wrap("guide",scales="free_x") + geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted")+
         ylab("ECDF"),
       dpi = 400, width = 16, height=16
)

figure("p_fitted_density",
       ggplot(dt, aes(x= p_ko_fitted, linetype=guide_detected, color=LL_method, xintercept=theta0))+
         geom_density(aes(y=..scaled..),fill=NA, size=1,alpha=0.5) + facet_wrap("guide",scales="free_x") +
         geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted"),
       dpi = 400, width = 16, height=16
)

figure("p_density",
       ggplot(dt, aes(x= p_ko, linetype=guide_detected, color=LL_method, xintercept=theta0))+
         geom_density(aes(y=..scaled..),fill=NA, size=1,alpha=0.5) + facet_wrap("guide",scales="free_x") +
         geom_vline(aes(xintercept=theta0), data=unique(dt, by=c("guide")), linetype="dotted"),
       dpi = 400, width = 16, height=16
)



for (method in c()){#"mvnorm", "norm", "dixit","same","xgboost")) {
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
