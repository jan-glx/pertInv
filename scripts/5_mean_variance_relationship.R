library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix#[1:1000,1:100]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)

dt = data.table(melt(count_matrix))
setnames(dt, c("cell","gene","counts"))
dt[,total_counts:=sum(counts),by=.(cell)]
nf_ <- edgeR::calcNormFactors(t(count_matrix))


dt[,nf:=nf_[.GRP],by=.(cell)]

summary_genes = rbind(dt[,.(mean=mean(counts),var=var(counts),method="plain"),by=.(gene)],
                      dt[,.(mean=mean(counts/total_counts*mean(total_counts)),var=var(counts/total_counts*mean(total_counts)),method="total_counts normalized"),by=.(gene)],
                      dt[,.(mean=mean(counts*nf),var=var(counts*nf),method="nf normalized"),by=.(gene)])

figure("mean-variance relationship of genes",
       ggplot(summary_genes,aes(x=mean,y=var,color=method))+geom_point(alpha=0.3)+
         scale_x_log10()+scale_y_log10()+geom_abline(color="black")+geom_smooth()+
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
