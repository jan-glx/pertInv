library(pertInv)
# --------------
# DAG:
#             X_5
#              +
#              v
#  X_2        X_4
#   +          +
#   +---->Y<---+
#   |     +
#   v     v
#  X_1   X_3
set.seed(2)
N=1000
pb =  progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                 total =  N,
                                 clear = FALSE, width= 60)
res.dt <- rbindlist( replicate(N,{
# inner loop --------------
n <- 50

dt <- local({
  dts <- list()

  e <- 1
  null_model <- list(
    X_5 = quote(X_5 <- rnorm(n)),
    X_4 =  quote(X_4 <- (3*X_5 + rnorm(n)) / sqrt(10)),
    X_2 = quote(X_2 <- rnorm(n)),
    Y = quote(Y <- (-3*X_2 + 3*X_4 + rnorm(n)) /sqrt(19)),
    X_3 = quote(X_3 <- (3*Y + rnorm(n)) / sqrt(10)),
    X_1 =  quote(X_1 <- (3*X_2 + rnorm(n)) / sqrt(10)),
    add = quote(dts[[e]] <- data.table(e, Y, X_1, X_2, X_3, X_4, X_5))
  )

  e <- 1
  model <- copy(null_model)
  for (inst in model) eval(inst)

  e=2
  model <- copy(null_model)
  model[["X_2"]] <- quote(X_2 <- rnorm(n) + 1)
  model[["X_3"]] <- quote(X_3 <- rnorm(n) + 1)
  for (inst in model) eval(inst)

  e=3
  model <- copy(null_model)
  model[["X_4"]] <- quote(X_4 <- rnorm(n) + 1)
  for (inst in model) eval(inst)

  e=4
  model <- copy(null_model)
  model[["X_1"]] <- quote(X_1 <- rnorm(n) + 1)
  #model[["X_5"]] <- quote(X_5 <- rnorm(n) + 1)
  for (inst in model) eval(inst)
  rbindlist(dts)
})

dt[,e:=as.factor(e)]
dt[,i:=.I]

f <- sample(c(-1,-0.4,0.3,1.1,2,3,4), 1)
XX0 <- dt[, cbind(X_1, X_2, X_3, X_4, X_5)]
XX <-  matrix(log(1+rpois(n=prod(dim(XX0)),lambda=exp(XX0+f))), nrow = nrow(XX0))

R2_poisson <-  pmax(0,1-mean(matrixStats::colVars(log(1+exp(XX0+f)))/
  sapply(seq_len(ncol(XX)), function(i) var(log(1+rpois(n=nrow(XX0)*30,lambda=exp(XX0[,i]+f)))))))



YY0 <-  dt[, Y]
YY <-  log(1+rpois(length(YY0),lambda=exp(YY0+f)))
ExpInd <- dt[, e]
res <- tryCatch({
res <- InvariantCausalPrediction::ICP(XX, YY, ExpInd, selection="all", showAcceptedSets = F, showCompletion = F, test="ranks" )
data.table(id=sample(.Machine$integer.max,1),  f=f,R2_poisson, p.val=res$pvalues, node=c("X_1", "X_2", "X_3", "X_4", "X_5"))
}, error = function(e) return(data.table(id=integer(0),  f=numeric(0), R2_poisson=numeric(0), p.val=numeric(0), node=character(0))) )
pb$tick()
res
# end   ----------
}, simplify=FALSE))

res.dt[,mean_R2_poisson:=mean(R2_poisson),by=f]


res.dt[, hypotheses:="true null"]
res.dt[node %in% c("X_2", "X_4"), hypotheses:="true alternative"]



alpha= 0.1
figure("ICP size/power over counting noise",
ggplot(res.dt[hypotheses == "true null"], aes(x = mean_R2_poisson, y = as.numeric(p.val<alpha), color = hypotheses)) +
  stat_summary(data = res.dt[hypotheses != "true null"], fun.data = binom_mean_ci_jeffreys) +
  stat_summary(fun.data = binom_mean_ci_jeffreys) +
  geom_hline(linetype="dashed",yintercept=alpha, color="dimgray")+
  annotate("text", x = res.dt[,max(mean_R2_poisson)], y = alpha, hjust=1, vjust=1, color="dimgray",label = paste0("nominal size: ",alpha)) +
  expand_limits(y = 0) +
  ylab("rejection rate") +
  xlab("counting noise as fraction of total variance") +
  scale_color_manual(values = cool_warm(2),
                     breaks = c("true null","true alternative"),
                     labels = c("true null (size)","true alternative (power)")) +
  theme(legend.justification = c(0,1))
, width=7, height=4)
#  -------------

# label switichig ---------------------------------------------------------------------------------
set.seed(2)
N=1000
pb =  progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                 total =  N,
                                 clear = FALSE, width= 60)
res.dt <- rbindlist( replicate(N,{
  # inner loop --------------
  n <- 25

  dt <- local({
    dts <- list()


    null_model <- list(
      X_5 = quote(X_5 <- rnorm(n)),
      X_4 =  quote(X_4 <- (3*X_5 + rnorm(n)) / sqrt(10)),
      X_2 = quote(X_2 <- rnorm(n)),
      Y = quote(Y <- (-3*X_2 + 3*X_4 + rnorm(n)) /sqrt(19)),
      X_3 = quote(X_3 <- (3*Y + rnorm(n)) / sqrt(10)),
      X_1 =  quote(X_1 <- (3*X_2 + rnorm(n)) / sqrt(10)),
      add = quote(dts[[e]] <- data.table(e, Y, X_1, X_2, X_3, X_4, X_5))
    )

    e <- 1L
    model <- copy(null_model)
    for (inst in model) eval(inst)

    e=2L
    model <- copy(null_model)
    model[["X_2"]] <- quote(X_2 <- rnorm(n) + 1)
    model[["X_3"]] <- quote(X_3 <- rnorm(n) + 1)
    for (inst in model) eval(inst)

    e=3L
    model <- copy(null_model)
    model[["X_4"]] <- quote(X_4 <- rnorm(n) + 1)
    for (inst in model) eval(inst)

    e=4L
    model <- copy(null_model)
    model[["X_1"]] <- quote(X_1 <- rnorm(n) + 1)
    #model[["X_5"]] <- quote(X_5 <- rnorm(n) + 1)
    for (inst in model) eval(inst)
    rbindlist(dts)
  })

  dt[,i:=.I]

  f <- sample(c(0,0.1,0.2,0.3,0.4,0.5), 1)
  sel <- runif(nrow(dt))<f
  dt[sel, e:=1L]#sample(levels(e),sum(sel),replace=TRUE)
  XX0 <- dt[, cbind(X_1, X_2, X_3, X_4, X_5)]
  XX <-  XX0



  YY0 <-  dt[, Y]
  YY <-  YY0
  ExpInd <- dt[, as.factor(e)]
  res <- tryCatch({
    res <- InvariantCausalPrediction::ICP(XX, YY, ExpInd, selection="all", showAcceptedSets = F, showCompletion = F, test="ranks" )
    data.table(id=sample(.Machine$integer.max,1),  f=f, p.val=res$pvalues, node=c("X_1", "X_2", "X_3", "X_4", "X_5"))
  }, error = function(e) return(data.table(id=integer(0),  f=numeric(0), R2_poisson=numeric(0), p.val=numeric(0), node=character(0))) )
  pb$tick()
  res
  # end   ----------
}, simplify=FALSE))



res.dt[, hypotheses:="true null"]
res.dt[node %in% c("X_2", "X_4"), hypotheses:="true alternative"]



alpha= 0.1
figure("ICP size/power over label switching",
       ggplot(res.dt[hypotheses == "true null"], aes(x = 1-f, y = as.numeric(p.val<alpha), color = hypotheses)) +
         stat_summary(data = res.dt[hypotheses != "true null"], fun.data = binom_mean_ci_jeffreys) +
         stat_summary(fun.data = binom_mean_ci_jeffreys) +
         geom_hline(linetype="dashed",yintercept=alpha, color="dimgray")+
         annotate("text", x = res.dt[,max(f)], y = alpha*0.95, hjust=1, vjust=1, color="dimgray",label = paste0("nominal size: ",alpha)) +
         expand_limits(y = 0) +
         scale_y_continuous(expand=c(0,0)) +
         ylab("rejection rate") +
         xlab("knockout efficiency") +
         scale_x_reverse(labels = scales::percent) +
         scale_color_manual(values = cool_warm(2),
                            breaks = c("true null","true alternative"),
                            labels = c("true null (size)","true alternative (power)")) +
       theme(legend.justification = c(0,1))
       , width=7, height=4)






