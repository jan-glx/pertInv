# -----------
library(pertInv)
library(SCnorm)

load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

dt = data.table(melt(count_matrix))
setnames(dt, c("cell","gene","counts"))
dt[,total_counts:=sum(counts),by=.(cell)]

ggplot(dt[gene %in% sample(unique(dt,by="gene")[,gene], 5)],aes(x=total_counts,y=counts,color=gene))+geom_jitter()+geom_smooth(method="lm")
ggplot(dt[gene %in% sample(unique(dt,by="gene")[,gene], 5)],aes(x=total_counts,y=counts))+geom_bin2d()+facet_wrap("gene")

Y = count_matrix
X = guide_matrix

plotCountDepth(Data = t(count_matrix))

DataNorm <- SCnorm(Data = t(count_matrix), Conditions= rep(1,nrow(count_matrix)),
                   PrintProgressPlots = TRUE,
                   FilterCellNum = 10,
                   PropToUse = .1,
                   Thresh = .1,
                   ditherCounts = TRUE)


plotCountDepth(Data =DataNorm)
