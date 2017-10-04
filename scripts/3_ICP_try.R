library("InvariantCausalPrediction")
library("pertInv")

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")
covariates.dt = fread("results/covariates.dt.csv")

Y <- stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=unique(covariates.dt,by="cell")) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

ExpInd <- unique(covariates.dt,by="cell")[, guide] #interaction(batch, guide)]#target_gene]#
ExpInd[is.na(ExpInd)] <- "none"

#res <- InvariantCausalPrediction::ICP(X, Y, ExpInd,"exact", maxNoObs = 1000, maxNoVariables=12)
#res
ii<-0

pb = txtProgressBar(min = ii*10+1, max = min(ncol(Y_adj),(ii+1)*10), initial = 0, style=3)
resp= list()
res = list()
for (i in seq(ii*10+1, min(ncol(Y_adj),(ii+1)*10))) {
  YY <- Y_adj[,i]
  XX <- Y_adj[,-i]
  res[[colnames(Y_adj)[i]]] <-ICP(XX, YY, ExpInd,
                                  alpha=1,
                                  gof=0.000000001,
                 test = "exact",#function(x,y) t.test(x,y)$p.value,
                 showAcceptedSets = FALSE,
                 showCompletion = FALSE,
                 stopIfEmpty=FALSE) #hiddenICP(XX, YY, ExpInd, intercept = TRUE)
  tmp <- res[[colnames(Y_adj)[i]]]
  resp[[colnames(Y_adj)[i]]] <- data.table(
    "cause"=tmp$colnames[tmp$usedvariables],
    "effect"= colnames(Y_adj)[i],
    "pvalue"=tmp$pvalues[tmp$usedvariables])
  setTxtProgressBar(pb,i)
}
close(pb)
#save(res, file=paste0("results/ICP_run_",ii,".Rdata"))

# start_time <- proc.time()
# res <- hiddenICP(XX, YY, ExpInd, intercept = TRUE)
# proc.time()-start_time
#
# start_time <- proc.time()
# res <- ICP(XX, YY, ExpInd)
# proc.time()-start_time
#
#
# start_time <- proc.time()
# res <- ICP(XX, YY, ExpInd, tes="exact")
# proc.time()-start_time
#
#
# library(sisVIVE)
# cv.sisVIVE(YY,XX[,1], model.matrix(~ExpInd))


