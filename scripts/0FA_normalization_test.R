

# ---------------------
mmv <- copy(log(1+count_matrix))
mmv <- mmv[,order(colMeans(mmv))]
mmv <- mmv[order(rowMeans(mmv)),]
mm<- sapply((1:(dim(mmv)[1]/1000))-1,function(i) sapply((1:(dim(mmv)[2]/10))-1,function(j)
  mean((mmv[1+((i*1000):((i+1)*1000-1)),1+((j*10):((j+1)*10-1))]))))
ddm <- melt(as.data.table(mm)[,i:=.I],id.vars = "i")
ggplot(ddm,aes(x=i,y=variable,fill=log(value)))+geom_raster()+
  viridis::scale_fill_viridis(name="log mean")+
  xlab("low to high mean counts of cell")+
  ylab("low to high mean counts of gene")

#mmv <- mmv/colMeans(mmv)
#mmv <- mmv/rowMeans(mmv)
#mmv <- log(1+mmv*1E4)
ms<- sapply((1:(dim(mmv)[1]/1000))-1,function(i) sapply((1:(dim(mmv)[2]/10))-1,function(j)
  mean(sqrt(matrixStats:::colVars(mmv[1+((i*1000):((i+1)*1000-1)),1+((j*10):((j+1)*10-1))])))))

dd <- melt(as.data.table(ms/mm)[,i:=.I],id.vars = "i")

ggplot(dd,aes(x=i,y=variable,fill=(value)))+geom_raster()+
  viridis::scale_fill_viridis(name="rel sd ")+
  xlab("low to high mean counts of cell")+
  ylab("low to high mean counts of gene")

mv<- sapply((1:(dim(mmv)[1]/1000))-1,function(i) sapply((1:(dim(mmv)[2]/10))-1,function(j)
  mean((matrixStats:::colVars(mmv[1+((i*1000):((i+1)*1000-1)),1+((j*10):((j+1)*10-1))])))))

dd <- melt(as.data.table(mv/mm)[,i:=.I],id.vars = "i")

ggplot(dd,aes(x=i,y=variable,fill=(value)))+geom_raster()+
  viridis::scale_fill_viridis(name="rel var ")+
  xlab("low to high mean counts of cell")+
  ylab("low to high mean counts of gene")





#----------------------
counts.dt[,y:=log2(1+count)]
dt <- merge(covariates.dt[,.(tmp=1L,cell_id2,guide, target_gene , batch, CDR, total_counts_scaled)],
            genenames.dt[,.(tmp=1L,gene_id2)],on="tmp",all=TRUE,allow.cartesian = TRUE)
dt <-counts.dt[dt,on=c("gene_id2","cell_id2")]
dt[is.na(count),count:=0]
dt[,gene:=factor(gene_id2)]
dt[CDR>0.166,CDR_group:="high"]
dt[CDR<0.133,CDR_group:="low"]
dt[CDR>=0.133 & CDR<=0.166,CDR_group:="medium"]


figure("Mean-Variance relationship across CDR",
       ggplot(dt[,.(mean=mean(y), std=sd(y)),by=.(gene_id,CDR_group)],aes(x=mean,y=std,color=CDR_group))+
         geom_point()+geom_smooth(method = 'loess')+
         expand_limits(y=0)
)


# expression matrix -----------
transformation_method = "Ancombes" #"Log(1+counts)" # "Log(1+counts)" #"Ancombes" #"Log(1+counts)"

Y <- switch(transformation_method,
            Ancombes=stabilize_Anscombes(count_matrix),
            `Log(1+counts)`=log2(1+count_matrix)
)
Y[1:3,1:4]

X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
anova(fit)
Y_adj <- residuals(fit)

gene_gene <- cor(Y)
gene_gene_adj <- cor(Y_adj)
if (!exists("new_order")){
  new_order <- (hclust(dist(gene_gene)))$order
}

diag(gene_gene) <- NA
figure1(paste(transformation_method,"raw"),
        plot_corr_matrix(gene_gene, new_order), title=FALSE)


diag(gene_gene_adj) <- NA
figure1(paste(transformation_method,"adjusted"),
        plot_corr_matrix(gene_gene_adj,new_order), title=FALSE)

figure1(paste(transformation_method, " adjusted p-val distribution"),
        hist(rcorr(Y_adj, type="pearson")$P), TRUE)

figure1(paste(transformation_method, " adjusted null p-val distribution"),
        hist(rcorr(apply(Y_adj, 2, sample) , type="pearson")$P), TRUE)

figure1(paste(transformation_method, " raw p-val distribution"),
        hist(rcorr(Y, type="pearson")$P), TRUE)

figure1(paste(transformation_method, " raw null p-val distribution"),
        hist(rcorr(apply(Y, 2, sample) , type="pearson")$P), TRUE)

figure1(paste(transformation_method, "QQ plot"),
        qqplot(-log_rcorr(Y_adj, type="pearson")$log_P, -log_rcorr(Y, type="pearson")$log_P))

figure1(paste(transformation_method, "QQ plot null"),
        qqplot(-log_rcorr(apply(Y_adj, 2, sample), type="pearson")$log_P, -log_rcorr(apply(Y, 2, sample), type="pearson")$log_P))

figure1(paste(transformation_method, "QQ plot null"),
        qqplot(-log_rcorr(apply(Y_adj, 2, sample), type="pearson")$log_P, -log_rcorr(apply(Y, 2, sample), type="pearson")$log_P))

figure1(paste(transformation_method, "QQ plot data vs null, adjusted"),
        qqplot(-log_rcorr(Y_adj, type="pearson")$log_P, -log_rcorr(apply(Y_adj, 2, sample), type="pearson")$log_P))

figure1(paste(transformation_method, "QQ plot data vs null, not adjusted"),
        qqplot(-log_rcorr(Y, type="pearson")$log_P, -log_rcorr(apply(Y, 2, sample), type="pearson")$log_P))

# ----------------------------------
guides <- covariates.dt[,guide]
guides[is.na(guides)] <-"none"
guides_f <- factor(guides)

X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide
fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)
pvals_kruskal_adj <- apply(Y_adj, MARGIN=2L, function(x) kruskal.test(x, guides_f)$p.val)
pvals_kruskal <- apply(Y, MARGIN=2L, function(x) kruskal.test(x, guides_f)$p.val)
pvals_kruskal_adj_0 <- apply(apply(Y_adj, 2, sample), MARGIN=2L, function(x) kruskal.test(x, guides_f)$p.val)

figure1(paste(transformation_method, "QQ plot: effect significance -"),
        qqplot(-log(pvals_kruskal_adj), -log(pvals_kruskal)))


figure1(paste(transformation_method, "QQ plot: effect significance -0"),
        qqplot(-log(pvals_kruskal_adj), -log(pvals_kruskal_adj_0)))


X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide
fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)
pvals_kruskal_adj <- apply(Y_adj, MARGIN=2L, function(x) kruskal.test(x, guides_f)$p.val)
pvals_kruskal <- apply(Y, MARGIN=2L, function(x) kruskal.test(x, guides_f)$p.val)

kruskal.test(covariates.dt$CDR, guides_f)

covariates.dt[is.na(target_gene),target_gene:="none"]
covariates.dt[,target_gene_f:=as.factor(target_gene)]
figure("CDR-guide effects",
       ggplot(covariates.dt,aes(x=guide,y=CDR,fill=target_gene))+
         # geom_segment(aes(x=CDR,xend=CDR,y=target_gene_f,
         #                 yend=as.integer(target_gene_f)-1,color=target_gene),alpha=0.3)
         #+
         geom_violin()+
         stat_summary(fun.data = "mean_cl_boot")+
         coord_flip())

figure1(paste(transformation_method, "QQ plot: effect significance -"),
        qqplot(-log(pvals_kruskal_adj), -log(pvals_kruskal)))

X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled,data=covariates.dt) #+guide
fit <- lm(Y ~ .,as.data.table(X))
Y_adj_no <- residuals(fit)
pvals_kruskal_no_batch <- apply(Y_adj_no, MARGIN=2L, function(x) kruskal.test(x, guides_f)$p.val)

figure1(paste(transformation_method, "QQ plot: effect significance - adj. batch/no batch"),
        {qqplot(-log(pvals_kruskal_adj), -log(pvals_kruskal_no_batch)); abline(0,1)})

#-------------------------------
library(Rtsne)

subset_ <- sample(nrow(Y), 10000)

res<- Rtsne(Y[subset_,])
res <- data.table(res$Y)
setnames(res,c("t-SNE1","t-SNE2"))
res <- cbind(res,covariates.dt[subset_])
figure("tSNE-unadjusted",
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=batch))+geom_point())


X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled,data=covariates.dt) #+guide
fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

res<- Rtsne(Y_adj[subset_,])
res <- data.table(res$Y)
setnames(res,c("t-SNE1","t-SNE2"))
res <- cbind(res,covariates.dt[subset_])
figure("tSNE-adjusted (not_batch)",
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=batch))+geom_point()
)



X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide
fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

res<- Rtsne(Y_adj[subset_,])
res <- data.table(res$Y)
setnames(res,c("t-SNE1","t-SNE2"))
res <- cbind(res,covariates.dt[subset_])
figure("tSNE-adjusted (batch too)",
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=batch))+geom_point())


figure("tSNE-adjusted by targeted",
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=target_gene))+geom_point())
figure("tSNE-adjusted by is targeted",
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=is.na(target_gene)))+geom_point())
figure("tSNE-adjusted by CDR",
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=CDR))+geom_point()+
         scale_color_viridis())

