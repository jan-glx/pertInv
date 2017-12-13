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

dt = data.table(melt(count_matrix))
setnames(dt, c("cell","gene","count"))

#ggplot(dt[,.(CDR = mean(count>0), total_counts = sum(count)), by=.(cell)], aes(x=CDR, y=total_counts)) + geom_point()

wMUC <- count_matrix %*% (1/matrixStats::colVars(count_matrix))
wMUC <- mean(wMUC)/wMUC
wMUC <- exp(log(wMUC)-mean(log(wMUC)))




library(Rtsne)
transformation_method = "Log(1+Y)" #"wLog(1+Y)" #"Log(1+counts)" # "Log(1+counts)" #"Ancombes" #"Log(1+counts)"
size_normalized = "size_normalized" # "unormalized" # "size_normalized"

Y <- count_matrix
if(size_normalized=="size_normalized")
  Y <- sweep(Y,1,wMUC,"*")

Y <- switch(transformation_method,
            Ancombes=stabilize_Anscombes(Y),
            `Log(1+Y)`=log2(1+Y),
            `wLog(1+Y)`=sweep(log2(1+Y),2,
                                 log(matrixStats::colVars(Y)/colMeans(Y))/sqrt(matrixStats::colVars(log2(1+Y))),
                                 FUN="*")
)


Y[1:3,1:4]
subset_ <- sample(nrow(Y), 10000)
subset_ <- seq_len(nrow(Y))

res<- Rtsne(Y[subset_,])
res <- data.table(res$Y)
setnames(res,c("t-SNE1","t-SNE2"))
res <- cbind(res,covariates.dt[subset_])
figure(paste0("tSNE - ", size_normalized),
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=batch))+
         geom_point(size=1,shape=3)+
         scale_color_manual(labels = 1:6, values = cbbPalette[-1])+
         theme(legend.justification = c(0, 1), legend.position = c(0, 1)),
       sub_title=paste0(data_set,"-",transformation_method),
       width=4.5,height=4.5)
figure(paste0("tSNE - ", size_normalized, " - by CDR"),
       ggplot(res, aes(x=`t-SNE1`,y=`t-SNE2`,color=CDR))+geom_point()+
         viridis::scale_color_viridis(),
       sub_title=paste0(data_set,"-",transformation_method))



res.dt <- covariates.dt[,{
  res <- Rtsne(Y[.I,])
  res <- data.table(res$Y)
  setnames(res,c("t-SNE1","t-SNE2"))
  res[,CDR:=CDR]
  res[,guides:=guides]
  res[,n_guides:=n_guides]
}, by=batch]

figure(paste0("per batch tSNE - ", size_normalized, " - by CDR"),
       ggplot(res.dt, aes(x=`t-SNE1`,y=`t-SNE2`,color=CDR))+
         facet_wrap("batch") +
         geom_point()+
         viridis::scale_color_viridis(),
       sub_title=paste0(data_set, "-", transformation_method))
mcgs <- names(rev(sort(table(res.dt[,guides])))[c(2,5,9,7)])
res.dt[,common_guides:=factor(guides,levels=mcgs,labels=mcgs)]

figure(paste0("per batch tSNE - ", size_normalized, " - by common guide"),
  ggplot(res.dt, aes(x=`t-SNE1`,y=`t-SNE2`,color=common_guides))+
    facet_wrap("batch") +
    geom_point(stroke=0,aes(alpha=0.5-0.4*is.na(common_guides)))+
    scale_alpha(guide = 'none') +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0)),
  sub_title=paste0(data_set, "-", transformation_method)
)


fit <- lm(Y ~ .,as.data.table(model.matrix(~0+batch, data=covariates.dt)))
Y_adj <- residuals(fit)

res2 <- Rtsne(Y_adj[subset_,])
res2 <- data.table(res2$Y)
setnames(res2,c("t-SNE1","t-SNE2"))
res2 <- cbind(res2,covariates.dt[subset_])
figure(paste0("tSNE - batch corrected, ", size_normalized),
       ggplot(res2, aes(x=`t-SNE1`,y=`t-SNE2`,color=batch))+
         geom_point(size=1,shape=3)+
         scale_color_manual(labels = 1:6, values = cbbPalette[-1])+
         theme(legend.justification = c(0, 1), legend.position = c(0, 1)),
       sub_title=paste0(data_set,"-",transformation_method),
       width=4.5,height=4.5)

figure(paste0("tSNE - batch corrected, ", size_normalized, " - by CDR"),
       ggplot(res2, aes(x=`t-SNE1`,y=`t-SNE2`,color=CDR))+geom_point()+
         viridis::scale_color_viridis(),
       sub_title=paste0(data_set,"-",transformation_method))



res2
mcgs <- names(rev(sort(table(res2[,guides])))[c(2,5,9,7)])
res2[,common_guides:=factor(guides,levels=mcgs,labels=mcgs)]
res2[,p_sgGABPA_1:=guide_matrix[,colnames(guide_matrix)=="p_sgGABPA_1"]]

figure(paste0("tSNE - batch corrected", size_normalized, " - p_sgGABPA_1"),
       ggplot(res2, aes(x=`t-SNE1`,y=`t-SNE2`))+
         geom_point( size=1,shape=3,aes(color=p_sgGABPA_1))+
         scale_alpha(guide = 'none') +
         theme(legend.justification = c(0, 1), legend.position = c(0, 1))+
         scale_color_manual(name = "sgGABPA",labels = c("not detected","present"), values = cbPalette[]),
       sub_title=paste0(data_set, "-", transformation_method),
       width=4.5,height=4.5)
)
