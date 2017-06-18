library("InvariantCausalPrediction")


Y <- stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

i=3
YY <- Y_adj[,i]
XX <-Y_adj[,-i]
ExpInd <- covariates.dt[, target_gene] #interaction(batch, guide)]#target_gene]#
ExpInd[is.na(ExpInd)] <- "none"

#res <- InvariantCausalPrediction::ICP(X, Y, ExpInd,"exact", maxNoObs = 1000, maxNoVariables=12)
#res

pb = txtProgressBar(min = 0, max = ncol(Y_adj), initial = 0, style=3)
res = list()
for (i in seq_len(ncol(Y_adj))) {
  YY <- Y_adj[,i]
  XX <- Y_adj[,-i]
  res[[i]] <-ICP(XX, YY, ExpInd,
                 showAcceptedSets = FALSE,
                 showCompletion = FALSE,
                 stopIfEmpty=TRUE) #hiddenICP(XX, YY, ExpInd, intercept = TRUE)
  setTxtProgressBar(pb,i)
}
close(pb)
save(res, file="ICP_run.Rdata")

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


