library(pertInv)
# --------------

#
#               X_5
#               +
#               v
#   X_2        X_4
#   +          +
#   +---->Y<---+
#   |     +
#   v     v
#   X_6   X_3



n <- 3000
dt <- local({
  set.seed(1)
  dts <- list()

  e <- 1
  null_model <- list(
    X_5 = quote(X_5 <- rnorm(n)),
    X_4 =  quote(X_4 <- (2*X_5 + rnorm(n)) / sqrt(5)),
    X_2 = quote(X_2 <- rnorm(n)),
    Y = quote(Y <- (-2*X_2 + 2*X_4 + rnorm(n)) /sqrt(9)),
    X_3 = quote(X_3 <- (2*Y + rnorm(n)) / sqrt(5)),
    X_6 =  quote(X_6 <- (2*X_2 + rnorm(n)) / sqrt(5)),
    add = quote(dts[[e]] <- data.table(e, Y, X_2, X_3, X_4, X_5, X_6))
  )

  e <- 1
  model <- copy(null_model)
  for (inst in model) eval(inst)

  e=2
  model[["X_2"]] <- quote(X_2 <- rnorm(n) + 1)
  model[["X_3"]] <- quote(X_3 <- rnorm(n) + 1)
  for (inst in model) eval(inst)

  model <- copy(null_model)
  e=2
  model[["X_4"]] <- quote(X_4 <- rnorm(n) + 1)
  for (inst in model) eval(inst)

  rbindlist(dts)
})

dt[,e:=as.factor(e)]
dt[,i:=.I]

XX <-  dt[, cbind(X_2, X_3, X_4, X_5, X_6)]
YY <-  dt[, Y]
ExpInd <- dt[, e]

res <- InvariantCausalPrediction::hiddenICP(XX, YY, ExpInd,
    alpha=1,
    intercept=TRUE)
res

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






