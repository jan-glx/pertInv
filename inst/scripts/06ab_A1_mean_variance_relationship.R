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

raw_library_size = rowMeans(count_matrix)/mean(count_matrix)

sf.dt <-  data.table(
  cell = rownames(count_matrix),
  raw = 1,
  TMM = 1/(edgeR::calcNormFactors(t(count_matrix), method="TMM")*raw_library_size),
  RLE = 1/(edgeR::calcNormFactors(t(count_matrix), method="RLE")*raw_library_size),
  UQ = 1/(edgeR::calcNormFactors(t(count_matrix), p=0.97, method="upperquartile")*raw_library_size),
  MUC = mean(count_matrix)/rowMeans(count_matrix)
)
sf.dt[,wMUC:=count_matrix %*% (1/matrixStats::colVars(count_matrix))]
sf.dt[,wMUC:=mean(wMUC)/wMUC]
sf.dt[,wMUC:=exp(log(wMUC)-mean(log(wMUC)))]
sf.dt[,MUC:=exp(log(MUC)-mean(log(MUC)))]

cor(sf.dt[,-(1:2)])


sf.dtm <- melt(sf.dt,id.vars="cell", value.name="sf", variable.name="method")

summary_genes <- sf.dtm[,dt[copy(.SD),.(mean=mean(counts*sf), var=var(counts*sf)),by="gene",on="cell"], by=method]

#summary_genes[,{f <- isoreg(x=log(mean),log(var)/log(mean)); gam(f$yf)},by=method]

figure("mean-variance relationship of genes",
       ggplot(summary_genes, aes(x=mean, y=var, group=method, color=method)) +
         geom_abline(color="gray", linetype="dashed") + #Poisson
         geom_abline(intercept=log10(1), slope=2, color="gray", linetype="dotted") + #lognormal CV=1
         geom_abline(intercept=log10(1/2^2), slope=2, color="gray", linetype="dotted") +
         geom_abline(intercept=log10((1/4)^2), slope=2, color="gray", linetype="dotted") +
         geom_point(alpha=0.5, stroke=0, size=1) +
         scale_x_log10() + scale_y_log10() +
         ylab("variance") +
         geom_smooth(formula=y ~ poly(x, 3), method="gam", size=1.1, color="black",fill="NA")+
         geom_smooth(formula=y ~ poly(x, 3), method="gam") +
         theme(legend.justification = c(0, 1), legend.position = c(0.05, 1))+
         annotation_logticks(),
       width=4.5, height=3.5
)

figure("mean-variance relationship of genes - zoom",
       ggplot(summary_genes, aes(x=log(mean), y=log(log(var)/log(mean)), group=method, color=method)) +
         geom_point(alpha=0.5,stroke =0,size=1) +
         scale_x_log10(breaks = log10_minor_break()) +scale_y_log10(breaks = log10_minor_break()) +
         geom_smooth(formula=y ~ poly(x, 3), method="gam", size=1.1, color="black",fill="NA")+
         geom_smooth(formula=y ~ poly(x, 3), method="gam") +
       #  geom_hline(yintercept=log10(log(1)), color="black", linetype="dashed") +
         stat_function(fun=function(x) log10(log((log(10)+2*x)/x)), color="black", linetype="dotted") +
         stat_function(fun=function(x) log10(log((log(1)+2*x)/x)), color="black", linetype="dotted") +
         stat_function(fun=function(x) log10(log((log(0.1)+2*x)/x)), color="black", linetype="dotted") +
         coord_cartesian(ylim=c(0.1,1.2)),#+
       #  theme(legend.justification = c(1, 1), legend.position = c(1, 1)),
       width=6, height=4
)

figure("mean-CV relationship of genes",
       ggplot(summary_genes, aes(x=mean, y=sqrt(var)/mean, group=method, color=method)) +
         geom_hline(yintercept=1, color="gray", linetype="dotted") +
         geom_hline(yintercept=1/2, color="gray", linetype="dotted") +
         geom_hline(yintercept=1/4, color="gray", linetype="dotted") +
         geom_abline(intercept=0, slope=-1/2, color="gray", linetype="dashed") + #Poisson
         geom_point(stroke =0,size=1,alpha=0.5) +
         scale_x_log10() +scale_y_log10(breaks = log10_minor_break())+
         ylab("coefficient of variation (CV)") +
         geom_smooth( formula = y ~ s(x, bs = "cs"), method="gam", size=1.1, color="black")+
         geom_line(stat="smooth", formula = y ~ s(x, bs = "cs"), method="gam",size=1,alpha=1) +
         theme(legend.justification = c(0.0, 0), legend.position = "none")+
         annotation_logticks(),
       width=3.5, height=3.5
)
