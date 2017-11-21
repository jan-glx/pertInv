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
n_genes <- ncol(count_matrix)
n_cells <- nrow(count_matrix)
load(file = file.path(data_folder, "guide_matrix.RData"))
covariates.dt <- fread(file.path(data_folder, "covariates.dt.csv"))

Y <- sweep(count_matrix, 2, colMeans(count_matrix), "/") # stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~0+I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch, data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

lmfit <- lm.fit(x=guide_matrix, y=Y_adj)

dt <- melt(data.table(Y_adj, keep.rownames = TRUE), id.vars="rn")


dt <- covariates.dt[dt,on=c("cell"="rn")]

dt_means <- dt[,.(mean=mean(value)),by=.(variable,target_gene,guide)]
t_tests <- dt_means[!is.na(target_gene),{null<-mean[target_gene=="INTERGENIC"];.SD[target_gene!="INTERGENIC",t.test(mean,null)$p.value,by=target_gene]},by=variable]
ggplot(t_tests,aes(x=variable,y=target_gene,fill=(p.adjust(V1,"BH"))))+geom_raster()+viridis::scale_fill_viridis()


ggplot(dt_means,aes(x=variable,y=mean,color=target_gene,group=guide))+stat_summary()


dt[,rep:=as.integer(as.factor(guide)),by=target_gene]

ss_genes <- sample(dt[,unique(variable)],10)
figure("guide_effects_subset",
ggplot(dt[variable%in%ss_genes],aes(x=variable,y=value,color=as.factor(rep),group=guide))+
  stat_summary()+stat_summary(geom="line")+facet_grid(target_gene~.,scales="free_y")+
  geom_hline(yintercept = 0, linetype="dotted")+scale_color_discrete(name="guide variant")
)


dt[is.na(guide),guide:="unknown"]
dt[is.na(target_gene),target_gene:="unknown"]

intra_target_tests <- dt[target_gene!="unknown",.(p.value=kruskal.test(value,as.factor(guide))$p.value),by=.(variable,target_gene)]
intra_target_tests[,p.adj := p.adjust(p.value,method="BH")]
intra_target_tests[,min(p.adjust(p.value)),by=target_gene][,p.adj:=p.adjust(V1,method="BH")][p.adj<0.1]

#intra_target_tests[,p.adjust(p.adj_ind,method="BH")]
#intra_target_tests[p.adj<0.5]
#dt[,1,by=.(target_gene,guide)][,.N,by=target_gene]


# significant guide effects on TF level
wilcox_tests <-dt[!is.na(target_gene),{null<-value[target_gene=="INTERGENIC"];.SD[target_gene!="INTERGENIC",wilcox.test(value,null)$p.value,by=target_gene]}, by=variable]

setorder(wilcox_tests[,min(p.adjust(V1)),by=target_gene][,p.adj:=p.adjust(V1,"BH")],by="V1")[]

wilcox_tests[p.adjust(V1,"BH")<0.05,]

# significant guide effect on guide level

wilcox_tests <-dt[!is.na(target_gene),{null<-value[target_gene=="INTERGENIC"];.SD[target_gene!="INTERGENIC",wilcox.test(value,null)$p.value,by=guide]}, by=variable]

setorder(wilcox_tests[,min(p.adjust(V1)),by=guide][,p.adj:=p.adjust(V1,"BH")],by="V1")[]

wilcox_tests[p.adjust(V1,"BH")<0.05,]

