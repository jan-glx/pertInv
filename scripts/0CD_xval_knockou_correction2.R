
# import ------------
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

count_matrix <- count_matrix#[1:1000,1:100]#count_matrix#count_matrix[1:1000,1:100]
n_genes <- ncol(count_matrix)
n_cells <- nrow(count_matrix)
p <- n_genes
n <- n_cells


# adjustment ------------
wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))

Y <- count_matrix
Y <- sweep(Y,1,wMUC,"*")
Y <- log2(1+Y)#quantile_normalizen_cells(Y)

X = guide_matrix
XX <- data.table(guide_matrix, wMUC=log(wMUC))


isTrain <- sample(n_cells, ceiling(n_cells*2/3))
isNotTrain <- seq_len(n_cells)[-isTrain]
isTestInNotTrain <- sample(rep(c(FALSE,TRUE), ceiling(length(isNotTrain)/2)), length(isNotTrain))

fit <- lm(Y[isTrain,] ~., as.data.table(X[isTrain,]))

#glm.b_fits_sim <- apply(count_matrix[,1:4], 2, function(count_vector) MASS::glm.nb(count_vector ~., XX[,]))
#glm.b_fits_sim[[1]]$fitted.values <-  glm.b_fits_sim[[1]]$fitted.values+4*guide_matrix[,1]
#count_matrix <- sapply(glm.b_fits_sim, function(fit) simulate(fit)[[1]])

glm.b_fits <- apply(count_matrix[isTrain,1:4], 2, function(count_vector) MASS::glm.nb(count_vector ~., XX[isTrain,]))

residuals_train <- predict(fit, as.data.table(X[isTrain,])) - Y[isTrain,]

sigma_res <- apply(residuals_train, MARGIN = 2, sd)
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
LL_glm.nb <- function(fits, X, Y){
  logliks <- function(n, th, mu, y, w)
    w*(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
         y * log(mu + (y == 0)) - (th + y) * log(th + mu))
  LLS <- vapply(seq_along(fits), function(g) {
          logliks(n=length(Y), fits[[g]]$theta, mu=predict(fits[[g]], X, "response"), y=Y[,g], w=1)
    }, numeric(nrow(X)))
  rowMeans(LLS) # not sum as they are not independent
}

LL_correct_norm  <- LL_norm(residuals_test_correct)
LL_correct_dixit <- LL_dixit(residuals_test_correct)
LL_correct_same <- LL_same(residuals_test_correct)
LL_correct_poisson <- LL_poisson(residuals_test_correct)
LL_correct_glm.nb <- LL_glm.nb(glm.b_fits, XX[isNotTrain,], Y=count_matrix[isNotTrain, ])


if (!is.null(cl)) clusterExport(cl, ls())

ni <- 1:2#seq_len(ncol(X))
pb = progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                total =  length(ni),
                                clear = FALSE, width= 60)

workerfun <- function(guide) {
  X_  <-  guide_matrix[isNotTrain,]
  guide_detected <- X_[,guide]
  X_[,guide] <- !X_[,guide]
  guide_str <- colnames(guide_matrix)[guide]
  XX[,(guide_str):=!get(guide_str)]

  residuals_test_swapped <- stats::predict(fit, as.data.table(X_)) - Y[isNotTrain,]

  LL_swapped_norm  <- LL_norm( residuals_test_swapped)
  LL_swapped_dixit <- LL_dixit(residuals_test_swapped)
  LL_swapped_same  <- LL_same(residuals_test_swapped)
  LL_swapped_poisson  <- LL_poisson(residuals_test_swapped)
  LL_swapped_glm.nb <- LL_glm.nb(glm.b_fits, XX[isNotTrain,], Y=count_matrix[isNotTrain, ])
  XX[,(guide_str):=!get(guide_str)]

  dt1 <- data.table::data.table(
    guide_id = guide,
    guide = guide_str,
    guide_detected =  rep(guide_detected, 5),
    LL_swapped = c(LL_swapped_norm, LL_swapped_dixit, LL_swapped_same, LL_swapped_poisson, LL_swapped_glm.nb),
    LL_correct = c(LL_correct_norm, LL_correct_dixit, LL_correct_same, LL_correct_poisson, LL_correct_glm.nb),
    LL_method = rep(c("norm", "dixit", "same", "poisson", "glm.nb"), each=length(LL_swapped_norm)),
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
