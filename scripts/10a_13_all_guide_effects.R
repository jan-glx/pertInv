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
n_g <- ncol(count_matrix)
n_c <- nrow(count_matrix)
load(file = file.path(data_folder, "guide_matrix.RData"))
n_r <- ncol(guide_matrix)
covariates.dt <- fread(file.path(data_folder, "covariates.dt.csv"))

#guide_matrix <- guide_matrix[sample(n_c),]

Y <- log2(1+count_matrix) # stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~0+I(CDR^2)+I(CDR^3)+CDR+log(total_counts_scaled)+log(total_counts_scaled)^2+log(total_counts_scaled)^3+batch, data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

lmfit <- lm.fit(x=guide_matrix, y=Y_adj)

dt0 <- melt(data.table(Y_adj, keep.rownames = TRUE), id.vars="rn")


dt0 <- covariates.dt[dt0,on=c("cell"="rn")]

target_genes <- stringr::str_match(colnames(guide_matrix),"^(?:c|m|p)_(?:sg)?((?:.*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]
is_intergenic <- rowSums(guide_matrix[,target_genes=="INTERGENIC"])

dt <- rbindlist(lapply(11, function(r) {
  if (target_genes[r]=="INTERGENIC")
    return(data.table())
  detected <-  guide_matrix[,r]
  ret <- data.table(gene=colnames(count_matrix),
                    p.value_adj=sapply(seq_len(n_g), function(g) t.test(Y_adj[detected&!is_intergenic,g],Y_adj[is_intergenic&!detected,g])$p.value),
                    eff_adj=sapply(seq_len(n_g), function(g) diff(t.test(Y_adj[detected&!is_intergenic,g],Y_adj[is_intergenic&!detected,g])$estimate)),
                    ci_adj=sapply(seq_len(n_g), function(g) diff(t.test(Y_adj[detected&!is_intergenic,g],Y_adj[is_intergenic&!detected,g])$conf.int)),
                    p.value_raw=sapply(seq_len(n_g), function(g) t.test(    Y[detected&!is_intergenic,g],    Y[is_intergenic&!detected,g])$p.value)
                    )
  ret[,guide:=colnames(guide_matrix)[r]]
  ret[,target:=target_genes[r]]
  ret
}))

dt <- data.table(gene=colnames(count_matrix),overdispersion= matrixStats::colVars(count_matrix)/colMeans(count_matrix))[dt,on="gene"]
dt <- data.table(gene=colnames(count_matrix),powerinfo= ((matrixStats::colVars(count_matrix)-colMeans(count_matrix))/colMeans(count_matrix)))[dt,on="gene"]

GABPAdt <- dt[guide=="p_sgGABPA_1",]
GABPAdt[, p.ihw := IHW::adj_pvalues(IHW::ihw(p.value_adj, log(powerinfo), nbins=2, alpha=0.05))]

ggplot(GABPAdt, aes(x=-log10(p.value_adj),y=log(powerinfo),color=-log10(p.ihw)>1))+geom_point()

ggplot(dt[], aes(gene, guide)) +
  geom_raster(aes(fill = eff_adj,alpha=min(ci_adj)/ci_adj)) +
  scale_fill_gradientn(colors=cool_warm(500))


dt[guide=="p_sgGABPA_1",hist(p.value_adj)]
dt[guide=="p_sgGABPA_1",hist(p.value_adj)]

res <-  dt[, .(
  "BH with tec. cov." = sum(p.adjust(p.value_adj,method="BH")<0.05),
  "BH unadjusted" = sum(p.adjust(p.value_raw,method="BH")<0.05),
  "IHW with tec. cov." = sum(IHW::adj_pvalues(IHW::ihw(p.value_adj, powerinfo, alpha=0.05))<0.05),
  "IHW unadjusted" = sum(IHW::adj_pvalues(IHW::ihw(p.value_raw, powerinfo, alpha=0.05))<0.05),
  "IHW with tec. cov. per guide" = sum(IHW::adj_pvalues(IHW::ihw(p.value_raw, factor(guide), covariate_type="nominal", alpha=0.05))<0.05),
)]
res <- melt(res,id=integer(0),variable.name = "method",value.name="# discoveries (FDR<5%)")
res
figure("increasing power",
  ggplot(res, aes(x=method,y=`# discoveries (FDR<5%)`))+
    geom_col()+coord_flip(),
  width=4.5,height=2.25
)

dt[,p.adjj := IHW::adj_pvalues(IHW::ihw(p.value_adj, factor(guide), alpha=0.05))]
filtered <- unique(dt[IHW::adj_pvalues(IHW::ihw(p.value_adj, powerinfo, alpha=0.05))<0.05],by=c("gene","guide"))
filtered <-dt[gene %in% filtered[,gene]]
ggplot(filtered, aes(gene, guide)) +
  geom_raster(aes(fill = eff_adj/max(abs(eff_adj)),alpha=1-p.adjj)) +
  scale_fill_gradientn(limits=c(-1,1), colors=cool_warm(101))

dt1 <- rbindlist(lapply(seq_len(n_r), function(r) {
  detected <-  guide_matrix[,r]
  ret <- data.table(gene=colnames(count_matrix),
                    mean=colMeans(Y_adj[detected,]),
                    vars=matrixStats::colVars(Y_adj[detected,])
  )
  ret[,guide:=colnames(guide_matrix)[r]]
  ret[,target:=target_genes[r]]
  ret
}))


ggplot(dt1, aes(gene, guide)) +
  geom_raster(aes(fill = mean, alpha=min(vars)/vars)) +
  scale_fill_gradientn(colors=cool_warm(500))

em <- as.matrix(dcast(filtered, gene~guide, value.var="eff_adj")[,-1])

dim(em)
em[1:3,1:4]

pheatmap::pheatmap(t((em)),color =cool_warm(501), breaks=seq(-0.4,0.4,length.out=501),clustering_distance_rows="correlation", clustering_distance_columns="correlation")



hist(dt$p.value_adj)
hist(dt$p.value_raw)

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


# significant guide effects on TF level
wilcox_tests <-dt[!is.na(target_gene),{null<-value[target_gene=="INTERGENIC"];.SD[target_gene!="INTERGENIC",wilcox.test(value,null)$p.value,by=target_gene]}, by=variable]

setorder(wilcox_tests[,min(p.adjust(V1)),by=target_gene][,p.adj:=p.adjust(V1,"BH")],by="V1")[]

wilcox_tests[p.adjust(V1,"BH")<0.05,]

# significant guide effect on guide level

wilcox_tests <-dt[!is.na(target_gene),{null<-value[target_gene=="INTERGENIC"];.SD[target_gene!="INTERGENIC",wilcox.test(value,null)$p.value,by=guide]}, by=variable]

setorder(wilcox_tests[,min(p.adjust(V1)),by=guide][,p.adj:=p.adjust(V1,"BH")],by="V1")[]

wilcox_tests[p.adjust(V1,"BH")<0.05,]



