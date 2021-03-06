
# ---------------------------
library(pertInv)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# helper functions --
lp__vb <- function(fit_vb, e_vb=extract(fit_vb)) {
  sapply(seq_len(dim(e_vb[[1]])[1]), function(i) {
    pars <- index_sample(e_vb, i)
    log_prob(fit_vb, unconstrain_pars(fit_vb, pars))
  })
}
lp_R_term <- function(fit) {
  e <- extract(fit, pars="lp_R_terms")[[1]]
  e <- aperm(e, c(1,3,2))
  #dim(e) <- c(dim(e)[1], prod(dim(e)[2:3]))
  e
}
summarize_sample <- function(xx) {
  xx <- xx[!names(xx)%in% c("D","lp__","R","Y")]
  xx <- lapply(xx,function(x){
    if(!is.null(dim(x))&&length(dim(x))>1) {
      dim(x) <- c(dim(x)[1],prod(dim(x)[-1]))
    }
    x
  })
  melt(as.data.table(xx),id.vars=integer(0))[,.(mean=mean(value),median=median(value),lower=quantile(value,0.025),upper=quantile(value,0.975)),by=variable]
}

# clear/init results lists
mm <- list()
fit_mc <- list()
fit_vb <- list()
ee <- list()
# ---------------------------
set.seed(1)
with_mc <-  TRUE

dat <- list(n_c = 64, n_g = 16, n_r = 4);

ii <- 1
mm[[ii]] <-  stan_model_builder(system.file("stan_files","1C_simulate.stan"))
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, chains=1, iter = 1000, warmup=0, algorithm="Fixed_param")
ee[[ii]] <- extract(fit_mc[[ii]])

dat  <-  c(dat, index_sample(ee[[ii]], i=1))
R_true <- dat$R

ii <- 2
mm[[ii]] <-  stan_model_builder(system.file("stan_files","1C_fit_vec_nb_disc.stan"))
if (with_mc) fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
fit_vb[[ii]] <- vb(mm[[ii]], data = dat)



ee[[ii]] <- extract(fit_vb[[ii]])
if (with_mc) ee[[ii]] <- extract(fit_mc[[ii]])
datt <- lapply(index_sample(ee[[ii]], i=1),function(x){
  x <- as.array(x)
  if(length(x)>1) {
    dim(x)<-c(1,dim(x))
  }
  x})






dt <- rbind(
  summarize_sample(extract(fit_mc[[1]]))[,method:="prior"],
  summarize_sample(extract(fit_vb[[ii]]))[,method:="VB"],
  if (with_mc)   summarize_sample(extract(fit_mc[[ii]]))[,method:="MC"] else data.table(),
  summarize_sample(datt)[,method:="ground truth"]
  )

dt[,c("variable_class","number"):=tstrsplit(variable,'.V',fixed=TRUE,keep=1:2)]
dt[,number:=as.integer(number)]
dt[is.na(number),number:=0]

figure("parameter recovery VBMC",
ggplot(dt[number<11], aes(x=variable,y=mean,ymin=lower,ymax=upper,color=method))+
  geom_errorbar(width=0,position=position_dodge(0.2))+
  geom_point(position=position_dodge(0.2))+
  facet_wrap("variable_class",scales="free")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(1, 0),
        legend.title=element_blank())+
  ylab("value")+coord_flip(),
width=9,height=5)

# plot 1D paramters ---
dt.single <- dt[!(variable_class %in% dt[number==3, unique(variable_class)])][!(variable %in% c("logit_p_D_given_R.V1", "logit_p_D_given_R.V2"))]

figure("1D parameter recovery VB/MC",
ggplot(dt.single, aes(x=variable,y=mean,ymin=lower,ymax=upper,color=method))+
  geom_rect(data=dt.single[method=="prior"], xmin=-Inf,xmax=Inf,alpha=0.3,fill="gray",color=NA)+
  geom_hline(data=dt.single[method=="prior"], aes(yintercept=mean),color="gray")+
  geom_crossbar(data=dt.single[method %in% "ground truth"],aes(y=mean,color=method))+
  geom_errorbar(data=dt.single[method %in% c("VB","MC")],width=0,position=position_dodge(0.3))+
  geom_point(data=dt.single[method %in% c("VB","MC")],position=position_dodge(0.3))+
  facet_wrap("variable",scales="free", nrow=2,strip.position="right")+
 # coord_flip()+
  theme(axis.title.x=element_blank(),
        strip.text.y = element_text(angle = 90,margin = margin(0,0.1,0,0.1, "cm")),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        legend.position="none",
        )+
  ylab("value"),
width=9,height=5)



# plot nD paramters ---
dt.multi <- dt[(variable_class %in% dt[number==3, unique(variable_class)])][number<11][!(variable_class %in% c("logit_p_R_r", "lp_R_terms"))]
varsorted <- setorder(unique(dt.multi[method=="ground truth"], key="variable"),"mean")[,variable]
dt.multi[,variable_factor:=factor(variable,levels=varsorted),]
dt.multi[dt.multi[method=="ground truth",.(variable_coded=rank(mean),variable),by=variable_class],variable_coded:=variable_coded,on="variable"]

dt.prior <- dt.multi[method=="prior"][,lapply(.SD,mean),by=.(method,variable_class), .SDcols=c("mean","median","lower","upper")]
dt.prior[,variable_coded:=1L]

figure("n-D parameter recovery VB/MC",
       ggplot(dt.multi, aes(x=variable_coded,y=mean,ymin=lower,ymax=upper,color=method))+
         geom_point(data=dt.multi[method %in% c("VB","MC")],position=position_dodge(0.3))+
      #   geom_rect(data=dt.prior, xmin=-Inf, xmax=Inf, alpha=0.3,fill="gray",color=NA)+
         geom_hline(data=dt.prior, aes(yintercept=mean),color="gray")+
           geom_rect(data=dt.prior, xmin=-Inf, xmax=Inf, alpha=0.3,fill="gray",color=NA)+
         geom_crossbar(data=dt.multi[method %in% "ground truth"],aes(y=mean,color=method))+
         geom_errorbar(data=dt.multi[method %in% c("VB","MC")],width=0,position=position_dodge(0.3))+
         geom_point(data=dt.multi[method %in% c("VB","MC")],position=position_dodge(0.3))+
         facet_wrap("variable_class",scales="free", nrow=2)+
        scale_x_discrete(limits=varsorted)+
         theme(axis.title.x=element_blank(),
               strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               legend.title=element_blank(),
               legend.justification=c(1,0),
               legend.position=c(0.9,0.1))+
         ylab("value"),
       width=9,height=5)


# initialization from prior ---------------------------------------------------
inits <- index_samples(ee[[1]], sample(dim((ee[[1]])[[1]])[1]-1,4)+1)

ii <- 3
mm[[ii]] <-  stan_model_builder(system.file("stan_files","1C_fit_vec_nb_disc.stan"))
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, init=inits)
fit_vb[[ii]] <- vb(mm[[ii]], data = dat, output_samples=100, tol_rel_obj=0.001)
ee[[ii]] <- extract(fit_vb[[ii]])


# initialize knockout matrix with guide matrix
ii <- 4
mm[[ii]] <-  stan_model_builder(system.file("stan_files","1C_fit_vec_nb_disc.stan"))

# ----------------------------------
dt <- data.table()
dat$R <- R_true
N <- 1000
p_Z_given_D <- numeric(4)
n <- 500
pb = progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                total =  (n),
                                clear = FALSE, width= 60)
for (i in seq_len(n)) {
  point <- sample(interesting_points,1)[[1]]
  c_ <- point[1]# sample(dat$n_c, 1)
  r <-  point[2]
  dat$R <- R_true
  dat$R[c_,r] <- 0
  capture.output(fit_vb_t <- vb(mm[[ii]], data = dat, output_samples=N, tol_rel_obj=0.0005))
  dat$R[c_,r] <- 1
  capture.output(fit_vb_star <- vb(mm[[ii]], data = dat, output_samples=N, tol_rel_obj=0.0005))

  lp_D_X_Zt_at_X_given_D_Zstar <- lp__vb(fit_vb_t, extract(fit_vb_star))
  lp_D_X_Zstar_at_X_given_D_Zt <- lp__vb(fit_vb_star, extract(fit_vb_t))
  lp_D_X_Zt_at_X_given_D_Zt <- lp__vb(fit_vb_t, extract(fit_vb_t))
  lp_D_X_Zstar_at_X_given_D_Zstar <- lp__vb(fit_vb_star, extract(fit_vb_star))


  #log_p_Zstar_given_D_over_p_Zt_given_D <- log_mean_exp(lp_D_X_Zt_at_X_given_D_Zstar)-log_mean_exp(lp_D_X_Zstar_at_X_given_D_Zt)
  p_Z_given_D[1] <- 1/(1+exp((log_sum_exp(c(0,log_mean_exp(lp_D_X_Zstar_at_X_given_D_Zt - lp_D_X_Zt_at_X_given_D_Zt)))-log_sum_exp(c(0,log_mean_exp(lp_D_X_Zt_at_X_given_D_Zstar - lp_D_X_Zstar_at_X_given_D_Zstar))))))
  p_Z_given_D[2] <- 1/(1+exp(-(log_mean_exp(lp_D_X_Zt_at_X_given_D_Zstar)-log_mean_exp(lp_D_X_Zstar_at_X_given_D_Zt))))
  p_Z_given_D[3] <- mean(1/(1+exp(-(lp_D_X_Zstar_at_X_given_D_Zstar-lp_D_X_Zt_at_X_given_D_Zt)))<runif(N))
  p_Z_given_D[4] <- 1/(1+exp((log_mean_exp(-lp_D_X_Zt_at_X_given_D_Zt)-log_mean_exp(-lp_D_X_Zstar_at_X_given_D_Zstar))))

  dt <- rbind(dt,data.table(p_Z_given_D=p_Z_given_D, method=1:4, R=R_true[c_,r], D=dat$D[c_,r]))
  print(dt[,as.data.table(mean_sd(p_Z_given_D))[,total_acceptance:=mean(p_Z_given_D>runif(.N))],by=.(method,R,D)])
  pb$tick()
}

ggplot(dt,aes(y=1-p_Z_given_D,x=as.factor(method)))+
  geom_violin(color=NA,fill=cbbPalette[2],scale = "width")+
  geom_jitter_normal(shape=3, width=0.1,size=3,alpha=0.5,stroke=0,color="black")+
  #stat_summary(fun.ymin="mean",fun.ymax="mean",fun.y="mean",geom="errorbar",color="black",width=0.5,size=2)+
  stat_summary(fun.data=mean_se,color="red",size=2,geom="pointrange")+
  facet_grid(R~D,labeller = "label_both")

figure("computing posterior of discrete paramter",
ggplot(dt[method!=3],aes(y=1-p_Z_given_D,x=factor(method,1:4,c("stabilized","inverse","averaged direct","direct"))))+
  geom_violin(color=NA,fill=cbbPalette[2],scale = "width")+
  geom_jitter_normal(shape=3, width=0.1,size=3,alpha=0.5,color="black")+
  #stat_summary(fun.ymin="mean",fun.ymax="mean",fun.y="mean",geom="errorbar",color="black",width=0.5,size=2)+
  stat_summary(fun.data=pertInv::quantile_ci, color="red",size=1,geom="pointrange")+
  facet_grid(R~D,labeller = "label_both")+xlab("method")+ylab(latex2exp::TeX("P(R_{cr}|R_{-cr},D,Y)")),
width=6.5,height=4.5)
#+  scale_shape_identity()
#+  scale_shape_identity()
# ----------------------------------

# -----------
dat$R <- dat$D
dat$R <-  R_true
N <- 99
fit_vb_tm1 <- vb(mm[[ii]], data = dat, output_samples=N, tol_rel_obj=0.001)
#lp_tm1 <-  lp__vb(fit_vb_tm1)
#lp_R_term_tm1 <- lp_R_term(fit_vb_tm1)

Rt <- list()
n <- 30
p_acceptance <- numeric(n)
p_acceptance2 <- numeric(n)
p_acceptance3 <- numeric(n)
p_acceptance4 <- numeric(n)
positions <- list()
output <- list()
t <- 1
pb = progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                total =  (n),
                                clear = FALSE, width= 60)
interesting_points <- data.table(D=as.vector(dat$D),R=melt(R_true))[,.(list(c(R.Var1[1],R.Var2[1]))), by= .(R.value,D)][,V1]
for (t in seq_len(n)){
  point <- sample(interesting_points,1)[[1]]
  c_ <- point[1]# sample(dat$n_c, 1)
  r <-  point[2]#sample(dat$n_r, 1)
  cat("\tD:", dat$D[c_,r],"\tRtm1:",dat$R[c_,r],"\tRtrue:",R_true[c_,r])
  positions[[t]] <- c("c"=c_, "r"=r)
  output[[t]] <- capture.output(
    fit_vb_tm1 <- vb(mm[[ii]], data = dat, output_samples=N, tol_rel_obj=0.001,eta=1)
  )

  dat$R[c_,r] <- !dat$R[c_,r]
  output[[t]] <- c(output[[t]],capture.output(
    fit_vb_t <- vb(mm[[ii]], data = dat, output_samples=N, tol_rel_obj=0.001,eta=1)
  ))

  #lp_t <- lp__vb(fit_vb_t)
  #lp_R_term_t <- lp_R_term(fit_vb_t)
  #p_acceptance[t] <- 1/(1+exp(lp_R_term_tm1[c_,r]-lp_R_term_t[c_,r]))
  #p_acceptance[t] <- exp(lp_R_term_t[c_,r])
  #lp_D_X1_Z1 <- lp__vb(fit_vb_t, extract(fit_vb_t))
  #lp_D_X0_Z0 <- lp__vb(fit_vb_tm1, extract(fit_vb_tm1))
  lp_D_X_Ztm1_at_X_given_D_Zt <- lp__vb(fit_vb_tm1, extract(fit_vb_t))
  lp_D_X_Zt_at_X_given_D_Ztm1 <- lp__vb(fit_vb_t, extract(fit_vb_tm1))
  lp_D_X_Ztm1_at_X_given_D_Ztm1 <- lp__vb(fit_vb_tm1, extract(fit_vb_tm1))
  lp_D_X_Zt_at_X_given_D_Zt <- lp__vb(fit_vb_t, extract(fit_vb_t))
  log_p_Zt_given_D_over_p_Ztm1_given_D <- log_mean_exp(lp_D_X_Ztm1_at_X_given_D_Zt)-log_mean_exp(lp_D_X_Zt_at_X_given_D_Ztm1)
  p_acceptance2[t] <- exp(-log_p_Zt_given_D_over_p_Ztm1_given_D)

  lp_D_X_Ztm1_at_X_given_D_Zm1 <- lp__vb(fit_vb_tm1, extract(fit_vb_tm1))
  lp_D_X_Zt_at_X_given_D_Z <- lp__vb(fit_vb_t, extract(fit_vb_t))
  p_acceptance3[t] <- mean(-(lp_D_X_Zt_at_X_given_D_Z-lp_D_X_Ztm1_at_X_given_D_Zm1)<log(runif(N)))
  p_acceptance4[t] <- exp(log_mean_exp(-lp_D_X_Ztm1_at_X_given_D_Zm1)-log_mean_exp(-lp_D_X_Zt_at_X_given_D_Z))



  p_acceptance[t] <- 1/(1+exp(-(log_sum_exp(c(0,log_mean_exp(lp_D_X_Zt_at_X_given_D_Ztm1 - lp_D_X_Ztm1_at_X_given_D_Ztm1)))-log_sum_exp(c(0,log_mean_exp(lp_D_X_Ztm1_at_X_given_D_Zt - lp_D_X_Zt_at_X_given_D_Zt))))))

  #P_Zt_given_D <- 1/mean(1/(exp(lp_R_term_t[,c_,r])))
  #P_Zt_given_D <- 1/(2 * mean(1/c(exp(lp_R_term_t[,c_,r]),1-exp(lp_R_term_tm1[,c_,r]))))
  #P_Ztm1_given_D <- 1/(2 * mean(1/c(exp(lp_R_term_tm1[,c_,r]),1-exp(lp_R_term_t[,c_,r]))))
  #P_Ztm1_given_D <- 1/mean(1/(exp(lp_R_term_tm1[,c_,r])))
   # should be the same 1/(2*exp(log_mean_exp(log(1)-c(lp_R_term_t[,c_,r],log(1-exp(lp_R_term_tm1[,c_,r]))))))
  #1/(2*exp(log_mean_exp(log(1)-c(lp_R_term_tm1[,c_,r],log(1-exp(lp_R_term_t[,c_,r]))))))
  cat("\tp_acceptance:",p_acceptance[t], "\tp_acceptance2:", p_acceptance2[t], "\tp_acceptance3:",p_acceptance3[t], "\tp_acceptance4:",p_acceptance4[t])
  accpeted <- runif(1) < p_acceptance[t]
  if (accpeted) {
    cat("\t+\n")
    #lp_tm1 <- lp_t
    fit_vb_tm1 <- fit_vb_t
    #lp_R_term_tm1 <- lp_R_term_t
  } else {
    cat("\t-\n")
    dat$R[c_,r] <- !dat$R[c_,r]
  }
  Rt[[t]] <- dat$R>0
  pb$tick()
}
dat$R <- R_true
mR <- Reduce(`+`,Rt)/n
res <- data.table(D=melt(dat$D),Rtrue=as.vector(R_true),mR=as.vector(mR))
res[,mean(mR),by=.(D.value,Rtrue)]
res[data.table(do.call(rbind,interesting_points)),on=c("D.Var1"="V1","D.Var2"="V2")]

#-----------------

#check runtimes
lapply(fit_mc, function(fit_mc) rowMeans(sapply(fit_mc@sim$samples, function(sample) attr(sample, "elapsed_time"))))




ii <- 3
mm[[ii]] <-  stan_model_builder(system.file("stan_files","1_fit_vec.stan"))
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat)
ee[[ii]] <- extract(fit_mc[[ii]])

gRNA_effects_s <-  ee[[ii]]$gRNA_effects
dim(gRNA_effects_s) <-  c(dim(gRNA_effects_s)[1], prod(dim(gRNA_effects_s)[-1]) )
gRNA_effects_s <- apply(gRNA_effects_s, MARGIN=2,quantile, probs=c(0.05,0.95))

mean(as.vector(dat$gRNA_effects)>gRNA_effects_s[1,] & as.vector(dat$gRNA_effects)<=gRNA_effects_s[2,])


sqrt(sum((dat$gRNA_effects - colMeans(ee[[ii]]$gRNA_effects))^2))
cor(as.vector(dat$D),as.vector(dat$R))


ii <- 4
mm[[ii]] <-  mm[[ii-1]]
dat_dist <- dat
dat_dist$R <- dat_dist$D
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat_dist)
ee[[ii]] <- extract(fit_mc[[ii]])

sqrt(sum((dat$gRNA_effects - colMeans(ee[[ii]]$gRNA_effects))^2))

ggplot(
  melt(data.table(true_value=as.vector(dat$gRNA_effects),
                  guides_known=as.vector(colMeans(ee[[3]]$gRNA_effects)),
                  guides_not_known=as.vector(colMeans(ee[[4]]$gRNA_effects)))[,i:=.I],
       id.vars=c("i","true_value"),value.name="estimated_value",variable.name="method"),
  aes(x=true_value,y=estimated_value,color=method,group=as.factor(i)))+
  geom_line(color="gray")+geom_point()+geom_smooth(aes(group=method),method="lm")+geom_abline()


ii <- 5
mm[[ii]] <-  mm[[ii-1]]
inits = index_samples(ee[[1]], is=2:5)
fit_mc[[ii]] <- sampling(mm[[ii]], data = dat, init = inits)
ee[[ii]] <- extract(fit_mc[[ii]])

dt <- rbindlist(lapply(fit_mc[[ii]]@sim$samples,function(x) as.data.table(c(x,attr(x, "sampler_params")))[,chain__:=stringi::stri_rand_strings(1, 3)]))
dt[,i__:=seq_len(.N),by=chain__]
scalar_pars <- stringr::str_subset(colnames(dt), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
dt <- dt[,scalar_pars, with=FALSE]
figure("low_information_reparametrized2",
       ggplot(melt(dt[i__>1000],id.vars=c("energy__","chain__","i__")),aes(y=energy__,x=value))+geom_point()+facet_wrap("variable", scales="free_x")+
         geom_smooth(method="lm"))

scalar_pars <- stringr::str_subset(names(fit_mc[[ii]]@sim$samples[[1]]), "(?:[^\\]]|(?:\\[1\\])|(?:\\[1,1\\]))$")
pairs(fit_mc[[ii]],pars=scalar_pars[1:ceiling(length(scalar_pars)/2)])
pairs(fit_mc[[ii]],pars=scalar_pars[-(1:ceiling(length(scalar_pars)/2))])
pairs(fit_mc[[ii]],pars=c("E_c[1]","sd_E","lp__"))
pairs(fit_mc[[ii]],pars=c("sd_mu_X","mu_X_g[1]","mu_X_g[2]","mu_X","lp__"))

sstan = shinystan::as.shinystan(fit_mc[[ii]], pars=fit_mc[[ii]]@sim$pars_oi[!(fit_mc[[ii]]@sim$pars_oi %in% c("X", "X_train", "X_test"))])

shinystan::launch_shinystan(sstan)






dat <- index_sample(e[[ii]], 1)
dat$n_c  <-  ncol(dat$X) # number of cells
dat$n_g  <-  nrow(dat$X) # number of genes
dat$n_r  <-  ncol(dat$D) # number of sgRNAs
n <- dat$n_c * dat$n_g
dat$ii_test <-  sample(n, ceiling(n/5))
dat$n_test <- length(dat$ii_test)





# ----------------- snippets

# show merged source
cat(fit_mc[[ii]]@stanmodel@model_code, file="results/stan_model_code.stan")

# attach data to globalenv for shinystan to find
list2env(dat,globalenv())

# check predictions with histogramms # deprecated
for (i in names(ee[[ii]])) {
  tmp = ee[[ii]][[i]]
  dim(tmp) = c(dim(tmp)[1],prod(dim(tmp)[-1]))
  hist(tmp[,1], main=i)
  abline(v=dat[[i]][1],col="red")
}

#shinystan subset to not n_c x n_g parameters
sstan = shinystan::as.shinystan(fit_mc[[ii]], pars=fit_mc[[ii]]@sim$pars_oi[!(fit_mc[[ii]]@sim$pars_oi %in% c("X", "X_train", "X_test"))])
shinystan::launch_shinystan(sstan)


lapply(ee[[ii]], mean)

