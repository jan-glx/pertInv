library(pertInv)
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

count_matrix <- count_matrix#[1:1000,1:100]#count_matrix#count_matrix[1:1000,1:100]
p <- ncol(count_matrix)
n <- nrow(count_matrix)

top_expressed <- order(colSums(count_matrix), decreasing = TRUE)

dt <-  data.table(n_guides = rowSums(guide_matrix),
                  total_counts = rowSums(count_matrix),
                  top10_counts = rowSums(count_matrix[,top_expressed[seq_len(10)]]),
                  low100_counts = rowSums(count_matrix[,rev(top_expressed)[seq_len(100)]]),
                  top100_counts = rowSums(count_matrix[,top_expressed[seq_len(100)]])
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


ggplot(dt, aes(y=n_guides, x=total_counts))+
  geom_smooth()+
  geom_smooth(method="lm")+
 # geom_smooth(method="loess")+
  layer(geom="point", stat="identity", position=position_jitter_normal(height=0.030,width=0), data=function(dt) sample.dt(dt, 1000)) +scale_x_log10()+scale_y_log10()


ggplot(dt, aes(y=total_counts, x=top10_counts))+geom_point(data=sample.dt(dt,13000))+geom_smooth()+geom_abline(slope=dt[,sum(total_counts)/sum(top10_counts)])
ggplot(dt, aes(y=low100_counts, x=top10_counts))+geom_jitter_normal(data=sample.dt(dt,13000))+geom_smooth(method="lm")+scale_x_log10()+scale_y_log10()
+geom_abline(slope=dt[,sum(total_counts)/sum(top10_counts)])



total_counts = rowSums(count_matrix)

hist(total_counts)
library(fitdistrplus)
fitdistcens
fitdist(total_counts,"pois")
fitdist(total_counts,"nbinom")
descdist(total_counts, discrete = TRUE)
descdist(total_counts, discrete = FALSE)
descdist(log(total_counts), discrete = FALSE)
fit.total_counts <- fitdist(total_counts, "nbinom")
plot(fit.total_counts)
fit.total_counts <- fitdist(total_counts, "gamma")
plot(fit.total_counts)
descdist(total_counts)

descdist(log(count_matrix[,4]+0.1), discrete = FALSE)

descdist(count_matrix[,4], discrete = FALSE)
