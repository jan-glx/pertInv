library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix#[1:1000,1:100]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)

top10_expressed <- order(colSums(count_matrix))[seq_len(10)]

dt <-  data.table(n_guides = rowSums(guide_matrix),
                  total_counts = rowSums(count_matrix),
                  top10_counts = rowSums(count_matrix[,top10_expressed])
                  )



#ggplot(dt, aes(x=n_guides, y=total_counts))+
#  layer(geom="point", stat="identity", position=position_jitter_normal(width=0.15), data=function(dt) sample.dt(dt, 1000)) +
#  stat_summary(fun.data = quantile_ci, color="red", geom = "smooth")+
#  stat_summary(fun.data = quantile_ci, color="blue", geom = "smooth", fun.args=list(p=0.75)) +
#  stat_summary(fun.data = quantile_ci, color="green", geom = "smooth", fun.args=list(p=0.25))+
#  geom_hline(yintercept = quantile(dt[n_guides==1,total_counts], probs=c(0.25,0.5,0.75)), linetype="dashed", alpha=0.5)
figure("total_count-guide_detection relationship",
ggplot(dt, aes(x=n_guides, y=total_counts))+
  layer(geom="point", stat="identity", position=position_jitter_normal(width=0.15), data=function(dt) sample.dt(dt, 1000)) +
  layer( stat="identity", position="identity",data = dt[,quantile_ci(total_counts,p=c(0.25,0.5,0.75)),keyby=n_guides], aes(y=y,ymin=ymin,ymax=ymax,x=n_guides,color=as.factor(p)), geom="smooth") +
  geom_hline(yintercept = quantile(dt[n_guides==1,total_counts], probs=c(0.25,0.5,0.75)), linetype="dashed", alpha=0.5) + labs(colour = "quantile")
)
