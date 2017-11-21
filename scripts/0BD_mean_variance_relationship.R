# load data ------------------
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
# -----------

dt = data.table(melt(count_matrix))
setnames(dt, c("cell","gene","counts"))
dt[,total_counts:=sum(counts), by=.(cell)]


figure("variance explained by count noise (over genes)",
  ggplot(dt[,.(pois_var=mean(counts), var=var(counts)), by=.(gene)][, R2_pois:=pois_var/var], aes(x=R2_pois)) +
    geom_density() +
    xlab(latex2exp::TeX("R^2_{Pois}")) + xlim(0,1),
  height=3, width=8
)


sf.dt <-  data.table(
  cell = rownames(count_matrix),
  raw = 1,
  TMM = edgeR::calcNormFactors(t(count_matrix)),
  total_UMIs = mean(count_matrix)/rowMeans(count_matrix)
)
sf.dt[,mine:=count_matrix %*% (1/matrixStats::colVars(count_matrix))]
sf.dt[,mine:=mean(mine)/mine]


sf.dtm <- melt(sf.dt,id.vars="cell", value.name="sf", variable.name="method")

summary_genes <- sf.dtm[,dt[copy(.SD),.(mean=mean(counts*sf), var=var(counts*sf)),by="gene",on="cell"], by=method]

#summary_genes[,{f <- isoreg(x=log(mean),log(var)/log(mean)); gam(f$yf)},by=method]

figure("mean-variance relationship of genes",
       ggplot(summary_genes, aes(x=mean, y=var, group=method, color=method)) + geom_point(alpha=0.3,stroke=0,size=2) +
         scale_x_log10() + scale_y_log10() + geom_abline(color="black", linetype="dashed") +
         geom_abline(intercept=log10(0.1), slope=2, color="black", linetype="dotted") +
         geom_abline(intercept=log10(1), slope=2, color="black", linetype="dotted") +
         geom_abline(intercept=log10(10), slope=2, color="black", linetype="dotted") +
         geom_smooth(formula=y ~ poly(x, 3), method="gam", size=1.5, color="black",fill="NA")+
         geom_smooth(formula=y ~ poly(x, 3), method="gam") +
         theme(legend.justification = c(1, 0), legend.position = c(1, 0)),
       width=5, height=4
)

figure("mean-variance relationship of genes - zoom",
       ggplot(summary_genes, aes(x=log(mean), y=log(log(var)/log(mean)), group=method, color=method)) + geom_point(alpha=0.3) +
         scale_x_log10() +scale_y_log10() +
         geom_smooth(formula=y ~ poly(x, 3), method="gam", size=1.5, color="black",fill="NA")+
         geom_smooth(formula=y ~ poly(x, 3), method="gam") +
         theme(legend.justification = c(1, 0), legend.position = c(1, 0)),
       width=5, height=4
)
# #setorder(summary_genes,mean)
# summary_genes[,var_fitted := {f <- isoreg(x=log(mean),log(var)-log(mean)); t=f$yf;t[f$ord]=t;exp(t+log(mean))}, by=method]
# #ggplot(summary_genes,aes(x=mean))+geom_point(aes(y=var))+geom_point(aes(y=var_fitted),color="red")
# ggplot(summary_genes[method=="mine"],aes(x=log(mean)))+geom_point(aes(y=(var)/mean))+geom_point(aes(y=(var_fitted)/(mean)),color="red")
#
# summary_genes[,var_fitted := {f <- lm(log(var/mean^2)~log(mean)); exp(fitted.values(f))*mean^2}, by=method]
# ggplot(summary_genes[method=="mine"],aes(x=log(mean)))+geom_point(aes(y=log(var)))+geom_point(aes(y=log(var_fitted)),color="red")
# ggplot(summary_genes[method=="mine"],aes(x=log(mean)))+geom_point(aes(y=log(var/mean^2)))+geom_point(aes(y=log(var_fitted/mean^2)),color="red")
# ggplot(summary_genes[method=="mine"],aes(x=log(mean)))+geom_point(aes(y=log(var)/log(mean)))
#
# summary_genes[,var_fitted := {f <- mgcv::gam(log(log(var/mean)/mean)~log(mean)); exp(exp(fitted.values(f))*mean)*mean}, by=method]
# ggplot(summary_genes[method=="mine"],aes(x=(mean)))+geom_point(aes(y=log(var)))+geom_point(aes(y=log(var_fitted)),color="red")
# ggplot(summary_genes[method=="mine"],aes(x=log(mean)))+geom_point(aes(y=log(log(var/mean)/mean)))+geom_point(aes(y=log(log(var_fitted/mean)/mean)),color="red")
#
# figure("mean-variance relationship of genes",
#        ggplot(summary_genes, aes(x=mean, y=var, group=method, color=method)) + geom_point(alpha=0.3) +
#          scale_x_log10() + scale_y_log10() + geom_abline(color="black", linetype="dashed") +
#          geom_abline(intercept=log10(0.1), slope=2, color="black", linetype="dotted") +
#          geom_abline(intercept=log10(1), slope=2, color="black", linetype="dotted") +
#          geom_abline(intercept=log10(10), slope=2, color="black", linetype="dotted") +
#          geom_line(aes(y=var_fitted), size=1.5, color="black")+
#          geom_line(aes(y=var_fitted))+
#          #geom_smooth(formula=y ~ poly(x, 3),aes(y=var_fitted), method="gam", size=1.5, color="black",fill="NA")+
#          #geom_smooth(formula=y ~ poly(x, 3),aes(y=var_fitted), method="gam") +
#          theme(legend.justification = c(1, 0), legend.position = c(1, 0)),
#        width=5, height=4
# )
figure("mean-rel. variance relationship of genes",
       ggplot(summary_genes, aes(x=mean, y=var/mean, group=method, color=method)) + geom_point(alpha=0.3) +
         scale_x_log10() + scale_y_log10() + geom_hline(yintercept=1,color="black") +
         geom_smooth(formula=y ~ poly(x, 6), method="gam", size=1.5, color="black",fill="NA")+
         geom_smooth(formula=y ~ poly(x, 6), method="gam")+
         theme(legend.justification = c(1, 0), legend.position = c(1, 0)),
       width=5, height=4
)


dt[dt[,.(total_counts=sum(counts)), by=.(cell)][, `# UMIs`:= c("low","medium","high")[ceiling(rank(total_counts)/.N*3)]], `# UMIs`:=`# UMIs`, on=.(cell)]
summary_genes <- dt[,.(mean=mean(count), var=var(counts)),by=.(gene,`# UMIs`)]

figure("mean-variance relationship of genes by groups of cells",
       ggplot(summary_genes, aes(x=mean, y=var, color=`# UMIs`)) + geom_point(alpha=0.3) +
         scale_x_log10() + scale_y_log10() + geom_abline(color="black") + geom_smooth()+
         theme(legend.justification = c(1, 0), legend.position = c(1, 0))
)





dt[, log2c := log2(0.5+counts)]
summary_genes =  rbind(dt[,.(mean=mean(log2c),var=var(log2c),method="plain"),by=.(gene)],
                      dt[,.(mean=mean(log2c+log2(mean(total_counts)/total_counts)),var=var(log2c+log2(mean(total_counts)/total_counts)),method="total_counts normalized"),by=.(gene)],
                      dt[,.(mean=mean(log2c+log2(nf)),var=var(log2c+log2(nf)),method="nf normalized"),by=.(gene)])

figure("mean-variance relationship of log2 genes + 0.5 ",
       ggplot(summary_genes,aes(x=mean,y=var^0.25,color=method))+geom_point(alpha=0.3)+
         scale_x_log10()+scale_y_log10()+geom_abline(color="black")+geom_smooth()+
         theme(legend.justification = c(1, 0), legend.position = c(1, 0))
)
