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

sf.dtm <- melt(sf.dt,id.vars="cell", value.name="sf", variable.name="method")

#guide_matrix <- guide_matrix[sample(n_c),]
guide_matrix2 <- guide_matrix[sample(n_c),]

target_genes <- stringr::str_match(colnames(guide_matrix),"^(?:c|m|p)_(?:sg)?((?:.*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]
is_intergenic <- rowSums(guide_matrix[,target_genes=="INTERGENIC"])
is_intergenic2 <- rowSums(guide_matrix2[,target_genes=="INTERGENIC"])

dt <- rbindlist(lapply(1:11, function(r) {
  if (target_genes[r]=="INTERGENIC")
    return(data.table())
  detected <-  guide_matrix[,r]
  detected2 <-  guide_matrix2[,r]
  ret <- rbindlist(lapply(seq_len(n_g), function(g) {
                      sf.dtm[, {
                        Y <- log2(1+count_matrix[,g] *sf)
                      .(p.value_adj=t.test(Y[detected&!is_intergenic],Y[is_intergenic&!detected])$p.value,
                        p.value_null=t.test(Y[detected2&!is_intergenic2],Y[is_intergenic2&!detected2])$p.value
                      )
                      },by=method]}))
  ret[,gene:=colnames(count_matrix)]
  ret[,guide:=colnames(guide_matrix)[r]]
  ret[,target:=target_genes[r]]
  ret
}))

dt <- data.table(gene=colnames(count_matrix),overdispersion= matrixStats::colVars(count_matrix)/colMeans(count_matrix))[dt,on="gene"]
dt <- data.table(gene=colnames(count_matrix),powerinfo= ((matrixStats::colVars(count_matrix)-colMeans(count_matrix))/colMeans(count_matrix)))[dt,on="gene"]

GABPAdt <- dt[guide=="p_sgGABPA_1",]
GABPAdt[, p.ihw := IHW::adj_pvalues(IHW::ihw(p.value_adj, log(powerinfo), nbins=2, alpha=0.05))]



#BiocInstaller::biocLite("IHW")
dt2 <- copy(dt)
dt2[,`:=`(mlp_1=-log10(sort(p.value_adj)), mlp_0= -log10(sort(p.value_null)),ii=.I),by=method]
ggplot(dt2,aes(color=method,x=mlp_1,y=mlp_0))+geom_line()+geom_abline()
ggplot(melt(dt2,measure.vars = c("mlp_0","mlp_1")),aes(color=method,x=1+ii,y=value,linetype=variable))+geom_line()+scale_x_log10()


res <-  melt(dt, measure.vars = c("p.value_adj", "p.value_null"))[, .(
  "BH" = sum(p.adjust(value,method="BH")<0.05),
  "IHW" = sum(IHW::adj_pvalues(IHW::ihw(value, factor(guide), alpha=0.05))<0.05)
  ),
  by=.(method,variable)
]
res <- melt(res,id=c("method","variable"),variable.name = "method2",value.name="# discoveries (FDR<5%)")
res

figure("power of normalization",
       ggplot(res, aes(x=interaction(variable,method2,method),y=`# discoveries (FDR<5%)`,fill=method2))+
         geom_col()+coord_flip(),
       width=4.5,height=2.25
)
