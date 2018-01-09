# load data -----------

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

#count_matrix <- count_matrix[,1:100]#count_matrix#count_matrix[1:1000,1:100]
n_genes <- ncol(count_matrix)
n_cells <- nrow(count_matrix)
p <- n_genes
n <- n_cells
guide_matrix <- guide_matrix[1:n_cells,]
batch_matrix <- batch_matrix[1:n_cells,]


dt <- data.table(melt(count_matrix))
setnames(dt, c("cell","gene","count"))

# transform --------------

Y = log2(1+count_matrix)
X = guide_matrix


wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))


# deconvolve -------------

ii <- c(975L, sample(ncol(count_matrix), 4), 1991L)
ii <- ii[c(1,2,4,6)]
dt <- rbindlist(lapply(ii, function(i) {
  y = count_matrix[,i]
  as.data.table(deconv(exp(seq(log(mean(y)/20), log(mean(y)+6*sd(y)), length.out = 200)), X = 1, y = y, #/wMUC
                       family ="PoissonExp", ignoreZero = FALSE, pDegree=10,c0=20)$stats)[,gene:=colnames(count_matrix)[i]]
}))
# plot results -------
dt2 <-data.table(melt(count_matrix[,ii]))
setnames(dt2, c("cell","gene", "count"))

dt[,panel:="1"]
dt[,density:="lambda (estimated)"]
dt[,`:=`(lambda=exp(theta),f=g/exp(theta),f.min=(g-SE.g)/exp(theta),f.max=(g+SE.g)/exp(theta))]
dt[,`:=`(f_scaled=f/max(f),f.min_scaled=f.min/max(f),f.max_scaled=f.max/max(f)),by=gene]
dt[,`:=`(g_scaled=g/max(g),g.min_scaled=(g-SE.g)/max(g),g.max_scaled=(g+SE.g)/max(g)),by=gene]

dt2[,panel:="2"]
dt2[,density:="Y (empirical)"]
dt[, gene_short :=tstrsplit(as.character(gene),"_")[[2]]]
dt2[, gene_short :=tstrsplit(as.character(gene),"_")[[2]]]

dt3 <- dt2[,.(..count..=.N), by=c("gene","gene_short","count","density")]
dt3[,ncount:=..count../max(..count..),by="gene"]


figure("Bayesian deconvolution",
       ggplot(data = dt3, aes(color=density,fill=density)) +
         geom_bar(aes(x=count, y=ncount), stat="identity")+
         geom_ribbon(data=dt, aes(x=theta,ymin=g.min_scaled,ymax=g.max_scaled),color=NA,alpha=0.5)+
         geom_line(data=dt, aes(x=theta,y=g_scaled))+
         labs(x = expression(lambda), y = "scaled density") +
         facet_wrap("gene_short", scales="free_x",nrow=2) +
         scale_color_manual(values=cbPalette[-1]) +
         scale_fill_manual(values=cbPalette[-1]) +
        theme(legend.justification = c(1, 1), legend.position = c(1, 1)),
       width=8,height=4
)

