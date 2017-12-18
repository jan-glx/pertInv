# -----------
library(pertInv)
library(glmnet)
library(mvtnorm)
library(caret)
library(boot)

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
batch_matrix <- batch_matrix[seq_len(n_cells),,drop=F]
batch_matrix <- batch_matrix[,colSums(batch_matrix)>0,drop=F]>0
guide_matrix <- guide_matrix[seq_len(n_cells),,drop=F]
guide_matrix <- guide_matrix[,colSums(guide_matrix)>0,drop=F]


target_genes <- stringr::str_match(colnames(guide_matrix),"^(?:c|m|p)_(?:sg)?((?:.*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]
is_intergenic <- rowSums(guide_matrix[,target_genes=="INTERGENIC"])

guide_matrix <- guide_matrix[,target_genes!="INTERGENIC",drop=F]
#guide_matrix <- guide_matrix[,1:3]

target_genes <- stringr::str_match(colnames(guide_matrix),"^(?:c|m|p)_(?:sg)?((?:.*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]

wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))

Y <- count_matrix
Y <- sweep(Y,1,wMUC,"*")
Y <- log2(1+Y)#quantile_normalizen_cells(Y)


capture = local({
  cdr = rowMeans(count_matrix[,]>0)
  mean_count = log2(1+rowMeans(count_matrix[,]))
  as.matrix(data.table(cdr,cdr^2,cdr^3,mean_count,mean_count^2,mean_count^3))
})

X_ <- if (ncol(batch_matrix)>1) data.table(guide_matrix, wMUC=log(wMUC),capture, batch_matrix[,-1]) else data.table(guide_matrix, wMUC=log(wMUC),capture)




guides = colnames(guide_matrix)
n_guides = length(guides)

# -------
LL <- function(res, sigma) apply(res, MARGIN=1, function(x) sum(dnorm(x,  sd=sigma,log = TRUE)))


LLR_knockout = function(Y, X_, train_cells, test_cells, pb) {
  calib_cells_idx <- sample(length(train_cells), round(length(train_cells)/5))
  calib_cells <- train_cells[calib_cells_idx]
  train_cells <- train_cells[-calib_cells_idx]
  set = character(nrow(X_))
  set[train_cells] = "train"
  set[test_cells] = "test"
  set[calib_cells] = "calib"

  fit = lm(Y[train_cells, ] ~.+1, X_[train_cells,])
  #fits = glmnet::cv.glmnet(X_, Y, family="mgaussian")


  beta = coef(fit) # coef(fits, s="lambda.min")
  sigma_res <- matrixStats::colSds(residuals(fit))
  LL_same <- function(res, y)  LL(res, mean(sigma_res))

  residuals_correct <- predict(fit, X_[c(test_cells,calib_cells),]) - Y[c(test_cells,calib_cells),]
  #residuals_correct <- predict(fits, X_[test_cells,], s="lambda.min") - Y[test_cells,]
  LL_correct = LL_same(residuals_correct)
  pb$tick()

  dt <- rbindlist(lapply(guides, function(guide) {
    guide_detected = guide_matrix[c(test_cells,calib_cells), guide]
    theta0 = (sum(guide_matrix[train_cells, guide])+0.5)/length(train_cells)
    residuals_swapped = residuals_correct - outer(2*guide_detected-1, beta[paste0(guide,"TRUE"),])
    LL_swapped = LL_same(residuals_swapped)
    LLR = LL_correct-LL_swapped
    LLR[!guide_detected] = -LLR[!guide_detected]
    LLR0 = log(theta0)-log(1-theta0)
    pb$tick()
    data.table(guide,LLR,theta0,LLR0,guide_detected, set=set[c(test_cells,calib_cells)], cell=rownames(guide_matrix)[c(test_cells,calib_cells)])
  }))

  dt[,scaling_factor := {
    fit_res <- optimize(f = function(sf) cross_entropy_loss_llr(sf*LLR[set=="calib"]+LLR0[set=="calib"], guide_detected[set=="calib"]), interval=c(-0.0001,5))
    fit_res$minimum
  }
  ]
  dt[set=="test",]
}


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
dt = rbind(helper(5,pb), helper(1,pb)[,set:="train"])
dt_bak <-  dt
# -----------
dt <-  copy(dt_bak)

ggplot(dt, aes(x=LLR,  linetype=set, color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
  facet_wrap("guide",scales ="free")

figure(
  paste0("distribution of LLR=LL_ko-LL_no_ko"),#\n(resampled)
  ggplot(dt, aes(fill=guide_detected,y=LLR,x=set))+geom_split_violin()+
    facet_wrap("guide",scales ="free"),
  height=15,width=15
)

dt[, p_ko:=1/(1+exp(-LLR))]
dt[, `p_ko (calibrated)`:=1/(1+exp(-(scaling_factor*LLR+LLR0)))]
ggplot(dt[set=="train", .(fraction_confidently_perturbed=mean(p_ko>0.9)), by=.(guide,guide_detected)], aes(x=guide_detected, y=fraction_confidently_perturbed)) + geom_violin()
ggplot(dt[set=="test", .(fraction_confidently_perturbed=mean(p_ko>0.60)), by=.(guide,guide_detected)], aes(x=guide_detected, y=fraction_confidently_perturbed)) + geom_violin()
ggplot(dt[set=="test", .(fraction_confidently_perturbed=mean(`p_ko (calibrated)`>0.90)), by=.(guide,guide_detected)], aes(x=guide_detected, y=fraction_confidently_perturbed)) + geom_violin()

dt <-   dt[set=="test",]

dt[, cross_entropy_loss_llr(scaling_factor*LLR+LLR0,guide_detected),by=.(set)]



dt[, cross_entropy_loss_llr(scaling_factor*LLR+LLR0,guide_detected),by=.(set)]

res = dt[, knnDemix::mixture.test(LLR[!guide_detected],LLR[guide_detected]),by=.(set,guide)]
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
  ggplot(dt[set=="test"][guide %in% res[p.adjust(p.value,method="BH")<0.08,guide]], aes(x=scaling_factor*LLR,color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
    facet_wrap("guide",scales ="free")+
    geom_text(data=res[p.adjust(p.value,method="BH")<0.08,],aes(
      label=sprintf("%s%%",signif(1-estimate,3)*100),#latex2exp::TeX(sprintf("%s\\%%",signif(1-estimate,3)*100), output="character"), #P(\\mathit{K_{r}}|\\mathit{D_{r}}=1) \n \\geq{}
      x=Inf,y=Inf), vjust=1.01, hjust=1.01, color="black", parse=FALSE)+
    xlab(latex2exp::TeX("log  \\frac{P(\\mathit{Y_{c}}|\\mathit{D_{cr}}=1)}{P(\\mathit{Y_{c}}|\\mathit{D_{cr}}=0)}"))+
    ylab("scaled densitiy") +
    scale_color_manual(name="", values=cbbPalette[2:3], breaks=c(TRUE, FALSE), labels=c("sgRNA detected", "not detected"))+
    scale_y_continuous(limits=c(0, 1), expand = expand_scale(c(0, 0),c(0, 0.02))) +# expand_scale(c(0, 0),c(0, 0.2))) +
    theme(legend.position=c(1,-0.25), legend.justification=c(1,0)),
  width=9, height=5
)

dt <- dt[res, on="guide"]
dt[,theta0_hat:=1-estimate]

dt[,guide_target:=stringr::str_match(guide,"^(?:c|m|p)_(?:sg)?((?:.*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]]
dt[,guide_variant:=as.integer(factor(guide)), by=guide_target]
dt[, `p_ko (lower bound)`:=1/(1+exp(-(scaling_factor*LLR+log(theta0_hat)-log(1-theta0_hat))))]

ggplot(dt[(guide_detected)], aes(x=`p_ko (lower bound)`,color=factor(guide_variant)))+facet_wrap("guide_target",scales="free_x")+geom_density(aes(y=..scaled..))

#dt[, `p_ko (lower bound)`:=1/(1+exp(-(log(pmax(1E-100,1+(exp(scaling_factor*LLR)-1)/theta0_hat))+log(theta0_hat)-log(1-theta0_hat))))]
#dt[, `p_ko (lower bound)`:=1-exp(log(estimate)+scaling_factor*LLR)]
figure("Estimated knockout probability",
ggplot(dt[(guide_detected)], aes(x=`p_ko (lower bound)`, group=factor(guide_variant)))+facet_wrap("guide_target",scales="free_x",ncol=5)+
  stat_density(aes(y=..scaled..), position="identity", geom="line", color= cbbPalette[3])+ylab("scaled density")+
  scale_y_continuous(limits=c(0, 1), expand = expand_scale(c(0, 0),c(0, 0.02))) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab(latex2exp::TeX("\\mathrm{estimated knockout probability in cell with sgRNA detected,}   \\hat{P}(\\mathit{K_{cr}}|\\mathit{Y_{c:},D_{cr}})")) +
  theme(axis.text.x=element_text(size=9))
, width=10,height=4)


dt
ggplot(dt[guide_target=="YY1"], aes(x=LLR,  linetype=set, color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
  facet_wrap("guide",scales ="free")
ggplot(dt[], aes(x=LLR,  linetype=set, color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
  facet_wrap("guide",scales ="free")


tmp <- dcast(dt[guide_target %in% c("YY1","GABPA")], cell~guide,value.var=c("LLR","guide_detected"))
tmp
ggplot(tmp, aes(x=LLR_p_sgYY1_10,y=LLR_p_sgYY1_3,  color= interaction(guide_detected_p_sgYY1_10, guide_detected_p_sgYY1_3)))+geom_density_2d()
ggplot(tmp, aes(x=LLR_p_sgYY1_10,y=LLR_p_sgYY1_3)) + facet_grid(guide_detected_p_sgYY1_10 ~ guide_detected_p_sgYY1_3)+geom_point()+geom_abline()



ggplot(tmp, aes(x=LLR_p_sgYY1_10,y=LLR_p_sgGABPA_9,  color= interaction(guide_detected_p_sgYY1_10, guide_detected_p_sgGABPA_9)))+geom_density_2d()+ facet_grid(guide_detected_p_sgYY1_10 ~ guide_detected_p_sgGABPA_9)+geom_abline()+geom_abline(slope=-1)




setorder(dt,-LLR)
roc <- dt[,.(true_positive_rate=cumsum(guide_detected)/sum(guide_detected),false_positive_rate=cumsum(!guide_detected)/sum(!guide_detected)),
          by=.(set)]
figure("ROC guide inference same",
       ggplot(roc, aes(x=false_positive_rate,y=true_positive_rate))+
         geom_line(color=cbbPalette[3])+
         geom_abline(linetype="dashed", color=cbbPalette[1])+
         xlab("false positive rate (1-specificity)")+
         ylab("true positive rate (sensitivity)")+
         scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) +
         scale_x_continuous(limits=c(0, 1), expand = c(0, 0))
       , width=4,height=4)

figure("ROC guide inference same (zoomed)",
       ggplot(roc, aes(x=false_positive_rate,y=true_positive_rate))+
         geom_line(color=cbbPalette[3])+
         geom_abline(linetype="dashed", color=cbbPalette[1])+
         xlab("false positive rate (1-specificity)")+
         ylab("true positive rate (sensitivity)")+
         coord_cartesian(xlim=c(0,0.1), ylim=c(0,0.2))
       , width=4,height=4)

#AUC
roc[,sum(true_positive_rate*c(0,diff(false_positive_rate)))]


setorder(dt,-LLR)
dt[,p.val:=(1+cumsum(!guide_detected))/(1+sum(!guide_detected)),by=guide]
dt
dtfdr <-  dt[(guide_detected), {res <- fdrtool::fdrtool(p.val, statistic="pvalue", plot=FALSE,verbose=FALSE)
  data.table(lfdr=res$lfdr, eta0=res$param[1,"eta0"], eta0.SE=res$param[1,"eta0.SE"], LLR=LLR)},by=c("guide","guide_target","guide_variant")]
res2 <- unique(dtfdr, by="guide")

figure(
  "LLR distribution of efficient guides fdr",
  ggplot(dt[set=="test"][guide %in% res2[eta0<0.70,guide]], aes(x=scaling_factor*LLR,color= guide_detected))+geom_density(aes(y=..scaled..),fill=NA)+
    facet_wrap("guide",scales ="free")+
    geom_text(data=res2[eta0<0.70,],aes(
      label=sprintf("%s%%",signif(1-eta0,3)*100),#latex2exp::TeX(sprintf("%s\\%%",signif(1-estimate,3)*100), output="character"), #P(\\mathit{K_{r}}|\\mathit{D_{r}}=1) \n \\geq{}
      x=Inf,y=Inf), vjust=1.01, hjust=1.01, color="black", parse=FALSE)+
    xlab(latex2exp::TeX("log  \\frac{P(\\mathit{Y_{c}}|\\mathit{D_{cr}}=1)}{P(\\mathit{Y_{c}}|\\mathit{D_{cr}}=0)}"))+
    ylab("scaled densitiy") +
    scale_color_manual(name="", values=cbbPalette[2:3], breaks=c(TRUE, FALSE), labels=c("sgRNA detected", "not detected"))+
    scale_y_continuous(limits=c(0, 1), expand = expand_scale(c(0, 0),c(0, 0.02))) +# expand_scale(c(0, 0),c(0, 0.2))) +
    theme(legend.position=c(1,-0.25), legend.justification=c(1,0)),
  width=9, height=5
)

figure1("fdrtool p_sgYY1_10",
        fdrtool::fdrtool(dt[guide=="p_sgYY1_10"][(guide_detected),p.val], statistic="pvalue", plot=TRUE),
        width=8, height=8
)

figure("Estimated knockout probability - lfdr",
       ggplot(dtfdr, aes(x=1-lfdr, group=factor(guide_variant)))+
         facet_wrap("guide_target",scales="free_x",ncol=5)+
         stat_density(aes(y=..scaled..), position="identity", geom="line", color= cbbPalette[3])+
         ylab("scaled density")+
         scale_y_continuous(limits=c(0, 1), expand = expand_scale(c(0, 0),c(0, 0.02))) +
         scale_x_continuous(expand = c(0, 0)) +
         xlab(latex2exp::TeX("\\mathrm{estimated knockout probability in cell with sgRNA detected,}   \\hat{P}(\\mathit{K_{cr}}|\\mathit{Y_{c:},D_{cr}})")) +
         theme(axis.text.x=element_text(size=9))
       , width=10,height=4)

figure("Estimated knockout probability - lfdr ECDF",
       ggplot(dtfdr, aes(x=1-lfdr, group=factor(guide_variant)))+
         facet_wrap("guide_target",scales="fixed",ncol=5)+
         stat_ecdf(color= cbbPalette[3],pad=FALSE)+
         ylab("ECDF")+
         coord_cartesian(xlim=c(0, 1),ylim=c(0,1) ) +
         #scale_x_continuous(limits=c(0, 1)) +
         xlab(latex2exp::TeX("\\mathrm{estimated knockout probability in cell with sgRNA detected,}   \\hat{P}(\\mathit{K_{cr}}|\\mathit{Y_{c:},D_{cr}})")) +
         theme(axis.text.x=element_text(size=9))
       , width=10,height=4)


dtfdr[,mean(1-lfdr>0.9),by=guide]
dtfdr[,mean(1/(1+exp(-LLR))>0.9),by=guide]

dt2 <- copy(dt_bak)
dt2[, set:=factor(set,c("train","test"))]

figure("dixit confidently perturbed",
ggplot(dt2[,mean(1/(1+exp(-LLR))>0.9),by=c("guide","set","guide_detected")], aes(y=V1,x=set,color=(guide_detected),fill=(guide_detected)))+
 # facet_grid(.~set) +
  geom_violin(position="identity",trim=F) +
  geom_violin(position="identity",fill=NA,trim=F) +
  geom_jitter(shape=21,color="black") +
  coord_cartesian(ylim=c(0, 1), expand=FALSE) +
  ylab("fraction \"confidently\" perturbed") +
  scale_color_manual(name="", values=cbbPalette[2:3], breaks=c(TRUE, FALSE), labels=c("sgRNA detected", "not detected")) +
  scale_fill_manual(name="", values=cbbPalette[2:3], breaks=c(TRUE, FALSE), labels=c("sgRNA detected", "not detected")) +
  theme(legend.position=c(1,1), legend.justification=c(1,1))
,width=4,height=4)

