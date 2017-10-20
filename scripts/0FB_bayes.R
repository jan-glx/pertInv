library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix[1:1000,1:100]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)

dt = data.table(melt(count_matrix))
setnames(dt, c("cell","gene","counts"))
dt[,total_counts:=sum(counts),by=.(cell)]
nf_ <- edgeR::calcNormFactors(t(count_matrix))

#
# y = rowSums(count_matrix)[1:10]
# n=10#nrow(count_matrix)
# xx=tcrossprod(count_matrix[1:10,])
# diag(xx) <- 0
# xx=rbind(xx,rep(1,n))
#
# MASS::ginv(xx) %*%c(rep(0,n),n)
# MASS::ginv(xx)




library(MCMCglmm)

#
# data <-  data.table(log2(1+count_matrix))
# responses = parse(text=paste0("cbind(",paste(colnames(count_matrix), collapse = ","),")"))[[1]]
#
# fit <- eval(substitute(
#   MCMCglmm(
#     responses ~ trait-1,
#     random = ~ units ,
#     data = data,
#     rcov = ~ idh(trait):units,
#     family=rep("gaussian", ncol(count_matrix)),
#     verbose = TRUE,
#     pr=TRUE,
#     nitt=13000, thin=10, burnin=3000,
#   ),
#   list(responses=responses)))
# sf_ <- 2^colMeans(fit$Sol)[-(1:n)]



data <-  data.table(count_matrix)
responses = parse(text=paste0("cbind(",paste(colnames(count_matrix), collapse = ","),")"))[[1]]
prior <- list(R = list(V = diag(p), nu = 0.002),
              G = list(G1= list(V = 1, nu = 0.002),
                       G2= list(V = 1, nu = 0.002)))#

fit <- eval(substitute(
  MCMCglmm(responses ~ 1,
           random =  ~ trait+units, # ~ units
           data = data,
           rcov = ~ idh(trait):units,
           family=rep("poisson", ncol(count_matrix)),
           prior = prior,
           pr=TRUE,
           verbose = TRUE),
  list(responses=responses)))
sf_ <- exp(colMeans(fit$Sol)[-(1:(p+1))])

dt[,sf:=sf_[.GRP],by=.(cell)]
dt[,nf:=nf_[.GRP],by=.(cell)]

summary_genes = rbind(dt[,.(mean=mean(counts),var=var(counts),method="plain"),by=.(gene)],
                      dt[,.(mean=mean(counts/total_counts*mean(total_counts)),var=var(counts/total_counts*mean(total_counts)),method="total_counts normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts/sf),var=var(counts/sf),method="MCMCglmm normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts*nf),var=var(counts*nf),method="nf normalized"),by=.(gene)])

figure("mean-variance relationship of genes MCMCglmm",
       ggplot(summary_genes,aes(x=mean,y=var,color=method))+geom_point(alpha=0.3)+
         scale_x_log10()+scale_y_log10()+geom_abline(color="black")+geom_smooth()+
         theme(legend.justification = c(1, 0), legend.position = c(1, 0))
)
