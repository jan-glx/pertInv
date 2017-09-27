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


library(glmmTMB)

ddt = as.data.table(melt(count_matrix))
setnames(ddt,c("cell","gene","counts"))

fit1 = glmmTMB(counts~1+(1|gene)+(1|cell), data=dt)
fit1
fit2 = update(fit1, family="poisson")
fit2
fit3 = update(fit2, family="nbinom2")
fit3
fit4 = update(fit3, dispformula = ~gene)
fit4
sf_ <- exp(ranef(fit4)$cond$cell[[1]])

library(glmnet)
fit <- cv.glmnet(cbind(sparse.model.matrix(~ gene-1,ddt),sparse.model.matrix(~ cell-1,ddt)),ddt$counts,family="poisson",alpha=0,lambda=1/3^(1:10),nfolds=5)
sf2_ <- exp(coef(fit)[-(1:(p+1))])


dt[,sf:=sf_[.GRP],by=.(cell)]
dt[,sf2:=sf2_[.GRP],by=.(cell)]
dt[,nf:=nf_[.GRP],by=.(cell)]

summary_genes = rbind(dt[,.(mean=mean(counts),var=var(counts),method="plain"),by=.(gene)],
                      dt[,.(mean=mean(counts/total_counts*mean(total_counts)),var=var(counts/total_counts*mean(total_counts)),method="total_counts normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts/sf),var=var(counts/sf),method="glmmTMB normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts/sf_),var=var(counts/sf_),method="glmnet normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts*nf),var=var(counts*nf),method="nf normalized"),by=.(gene)])

figure("Mean-variance relationship of genes glmmTMB _",
       ggplot(summary_genes,aes(x=mean,y=var,color=method))+geom_point(alpha=0.3)+
         scale_x_log10()+scale_y_log10()+geom_abline(color="black")+geom_smooth()+
         theme(legend.justification = c(1, 0), legend.position = c(1, 0))

)
