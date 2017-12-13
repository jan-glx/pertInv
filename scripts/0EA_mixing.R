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

Y = log2(1+count_matrix) #stabilize_Anscombes(count_matrix) #log2(1+count_matrix)
X = guide_matrix


wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))






# deconvolve -------------

ii <- c(975L, sample(ncol(count_matrix), 4), 1991L)#c(867L, 1461L, 871L, 1120L, 47L, 1991L)#c(1135L, 1134L, 1461L, 1133L, 1187L, 1876L)#sample(ncol(count_matrix), 6)
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


1







## LOG scale bullshit -----------------------


# deconvolve -------------

ii <- sample(ncol(count_matrix), 6)

dt <- rbindlist(lapply(ii, function(i) {
  y = count_matrix[,i]
  as.data.table(deconv(seq(log(mean(y)/30), log(mean(y)+4*sd(y)), length.out = 100), X = 1/wMUC, y = y,
                        family ="PoissonExpLog", ignoreZero = FALSE, pDegree=10,c0=20)$stats)[,gene:=colnames(count_matrix)[i]]
}))
# plot results -------
dt2 <-data.table(melt(count_matrix[,ii]))
setnames(dt2, c("cell","gene", "count"))

dt[,panel:="1"]
dt[,density:="lambda (estimated)"]
dt[,`:=`(lambda=exp(theta),f=g/exp(theta),f.min=(g-SE.g)/exp(theta),f.max=(g+SE.g)/exp(theta))]
dt[,`:=`(f_scaled=f/max(f),f.min_scaled=f.min/max(f),f.max_scaled=f.max/max(f)),by=gene]

dt2[,panel:="2"]
dt2[,density:="Y (empirical)"]
dt[, gene_short :=tstrsplit(as.character(gene),"_")[[2]]]
dt2[, gene_short :=tstrsplit(as.character(gene),"_")[[2]]]

dt3 <- dt2[,.(..count..=.N), by=c("gene","gene_short","count","density")]
dt3[,ncount:=..count../max(..count..),by="gene"]


figure("Bayes deconv: hist & g",
       ggplot(data = dt3, aes(color=density,fill=density)) +
         geom_bar(aes(x=count, y=ncount), stat="identity")+
         geom_ribbon(data=dt, aes(x=lambda,ymin=f.min_scaled,ymax=f.max_scaled),color=NA,alpha=0.5)+
         geom_line(data=dt, aes(x=lambda,y=f_scaled))+
         labs(x = expression(lambda), y = "scaled density") +
         facet_wrap("gene_short", scales="free_x",nrow=3) +
         scale_color_manual(values=cbPalette[-1]) +
         scale_fill_manual(values=cbPalette[-1])
)














ddt = rbindlist(list(dt,dt2),fill=TRUE)


ddt[,max_g:=max(g, na.rm = TRUE), by=gene]

#expression(hat(g)(lambda))
figure("Bayes deconv: hist & g",
       ggplot(data = ddt, aes(color=panel,fill=panel)) +
         stat_summary(data=ddt[panel=="2"], geom="bar",fun.data = function(data){data.table(x=1:5)[,.(..count..=.N),by=x][,y:=..count../max(..count..)]}, mapping=aes(x=count, y=count,color=panel,fill=panel))+
         geom_smooth(data= ddt[panel=="1"],stat="identity", mapping = aes(x = exp(theta), y = g/exp(theta), ymin = (g-SE.g)/exp(theta), ymax = (g+ SE.g)/exp(theta), color=panel)) +
         labs(x = expression(lambda), y = "scaled density") +
         facet_wrap("gene_short", scales="free_x",shrink=FALSE) +
         scale_color_manual(values=cbPalette[-1]) +
         scale_fill_manual(values=cbPalette[-1])
)

figure("Bayes deconv: hist & g (log scale)",
       ggplot(data = ddt) +
         geom_smooth(data= ddt[panel=="1"],stat="identity", mapping = aes(x = theta, y = g, ymin = g-SE.g, ymax =g+ SE.g, color=gene)) +
         labs(x = expression(log(theta)), y = expression(hat(g)(theta))) +
         geom_histogram(data=ddt[panel=="2"], aes(x=log(pmax(count,1/100)), y=..count.., color=gene), position=position_dodge(0.1),fill=NA,breaks = log(c(1/100,1:100)))+ #(log(c(1/100,1:10))+log(1:11))/2
         coord_cartesian(xlim=c(-5,5))+facet_grid(panel~.,scales="free_y")
)
figure("Bayes deconv: ECDF & G",
       ggplot(data = ddt) +
         geom_smooth(data= ddt[panel=="1"],stat="identity", mapping = aes(x = exp(theta), y = G, ymin = G-SE.G, ymax =G+ SE.G, color=gene)) +
         labs(x = expression(theta), y = expression(hat(G)(theta))) +
         stat_ecdf(data=ddt[panel=="2"], aes(x=count, color=gene))+ #(log(c(1/100,1:10))+log(1:11))/2   log(pmax(count,1/10000))
         coord_cartesian(xlim=c(0,10))
)
figure("Bayes deconv: ECDF & G (log scale)",
       ggplot(data = ddt) +
         geom_smooth(data= ddt[panel=="1"],stat="identity", mapping = aes(x = exp(theta), y = G, ymin = G-SE.G, ymax =G+ SE.G, color=panel)) +
         labs(x = expression(theta), y = expression(hat(G)(theta))) +
         stat_ecdf(data=ddt[panel=="2"], aes(x=pmax(count,0.0001), color=panel))+ #(log(c(1/100,1:10))+log(1:11))/2   log(pmax(count,1/10000))
         coord_cartesian(xlim=exp(range(dt$theta)))+scale_x_log10()+facet_grid(gene_short~.)
)

# ---------------

y = as.vector(count_matrix)
dt = as.data.table(deconv2(seq(log(1/100), log(10*mean(y)), length.out = 100), X= 1, y = y,
                           family ="PoissonExpLog", ignoreZero = FALSE, pDegree=10,c0=20)$stats)
# -------------

i <- sample(ncol(count_matrix), 1)
i_guide <-  sample(ncol(guide_matrix), 1)
ss <- guide_matrix[,i_guide]

X <- rowSums(count_matrix)
#X <- rowSums(count_matrix>0)
#X <- rep(1, nrow(count_matrix))

y <- count_matrix[,i]
# hist(y, breaks = 30)
f <- deconv2(seq(0, mean(y)*5, length.out = 100)[-1], X= (X/mean(X))[!ss], y = y[!ss],
             family ="PoissonExp", ignoreZero = FALSE, pDegree=10,c0=20)
f <- deconv2(seq(0, mean(y)*5, length.out = 100)[-1], X= (X/mean(X))[!ss], y = y[!ss], aStart = f$mle,
             family ="PoissonExp", ignoreZero = FALSE, pDegree=10,c0=5)
f2 <- deconv2(seq(0, mean(y)*5, length.out = 100)[-1], X= (X/mean(X))[ss], y = y[ss], aStart = f$mle,
              family ="PoissonExp", ignoreZero = FALSE, pDegree=10,c0=5)

#figure("Poisson-P-Spline-mixture (Bayesian deconvolution) exposure=1",
ggplot(data = rbind(as.data.table(f$stats)[,guide:=paste0(colnames(guide_matrix)[i_guide],"-")],as.data.table(f2$stats)[,guide:=paste0(colnames(guide_matrix)[i_guide],"+")])) +
  geom_smooth(stat="identity", mapping = aes(x = theta, y = g, ymin = g-SE.g, ymax =g+ SE.g,color=guide)) +
  labs(x = expression(theta), y = expression(hat(g)(theta)))+ggtitle(paste0(i,": ", colnames(count_matrix)[i]))#, paste0("Counts of gene ", i))








