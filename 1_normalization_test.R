# includes -------------
library(data.table)
library(cowplot)
library(stringr)
library(Matrix)
library(scran)

source("0_common.R")

# Load count data -----------------------------------
data_set = "GSM2396859_k562_tfs_13"
# "GSM2396861_k562_ccycle"
# "GSM2396858_k562_tfs_7"
# "GSM2396859_k562_tfs_13"
# "GSM2396860_k562_tfs_highmoi"
# "GSM2396856_dc_3hr"
# "GSM2396859_k562_tfs_13"

counts.dt <- fread(paste0('gzip -dc data/',data_set,'.mtx.txt.gz'))
setnames(counts.dt,c("gene_id","cell_id","count"))
n_genes <- counts.dt[1,gene_id]
n_cells <- counts.dt[1,cell_id]
n_not_no_reads <- counts.dt[1,count]
counts.dt <- counts.dt[-1]

# Load ID-cell and ID-gene mapping
cellnames.dt <- fread(paste0('gzip -dc data/',data_set,'_cellnames.csv.gz'), header=TRUE)
setnames(cellnames.dt, c("cell_id", "cell"))
cellnames.dt[,cell_id:=cell_id+1] # changing to one_based

cellnames.dt[,batch:=str_match(cell,".*_(.*_.*)")[,2]]
cellnames.dt[,cell_bc:=str_match(cell,"(.*)_.*_.*")[,2]]

genenames.dt <- fread(paste0('gzip -dc data/',data_set,'_genenames.csv.gz'), header=TRUE)
setnames(genenames.dt, c("gene_id", "gene"))
genenames.dt[,gene_id:=gene_id+1] # changing to one_based




# Load cell-guide mapping
cbc_gbc_dict.dt <-
  if(file.exists(paste0('data/',data_set,'_cbc_gbc_dict.csv.gz'))){
    fread(paste0('gzip -dc data/',data_set,'_cbc_gbc_dict.csv.gz'), header= FALSE)
  }else{
    fread(paste0('gzip -dc data/',data_set,'_cbc_gbc_dict_strict.csv.gz'), header= FALSE) # _strict _lenient
  }
setnames(cbc_gbc_dict.dt, c("guide", "cells"))
cbc_gbc_dict.dt[, target_gene:=str_match(guide,"^p_(?:sg)?((?:(?<=sg).*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]]

cbc_gbc_dict.dt <- cbc_gbc_dict.dt[,.(cell=str_split(cells, ", ")[[1]]), by=.(guide, target_gene)]
cbc_gbc_dict.dt[,.N,by=cell][,.N,keyby=N]
## add batch info and subset observed cells
#cbc_gbc_dict.dt <- cbc_gbc_dict.dt[cellnames.dt,on="cell"]

# find targt genes of guides
genenames.dt[, targted_by:=str_match(
  genenames.dt[,gene],
  paste0("(", paste0(cbc_gbc_dict.dt[,unique(target_gene)],collapse="|"), ")")
)[,2]]

# inspection --------------------

guides_per_cell_dt <- unique(cbc_gbc_dict.dt,by=c("cell","guide"))[,.(unique_guides_per_cell=sum(!is.na(guide))),by=cell]

figure("guides per cell",
       ggplot(guides_per_cell_dt,aes(x=as.factor(unique_guides_per_cell)))+geom_bar(aes(y=..count../sum(..count..)))+
         xlab("unique guides per cell (detected)")+ylab("relative frequency")
)

figure("distribution of number of cells with guide",
       ggplot(cbc_gbc_dict.dt[,.N,by=guide], aes(x=N))+
         geom_density(color='red')+
         geom_segment(aes(x=N,xend=N,y=0,yend=0.0002))+
         xlab('# cells with respective guide')+ expand_limits(x=0)
)

NUM_GENES_TO_KEEP = 2000


genes_summary.dt <-
  counts.dt[, .(total_counts=sum(count), cells_with_counts=.N),by=.(gene_id)
            ][, keep := rank(-cells_with_counts)<=NUM_GENES_TO_KEEP]

figure("Distribution counts over gene",
       ggplot(genes_summary.dt, aes(x=total_counts))+
         geom_density()+
         geom_segment(aes(x=total_counts,xend=total_counts,y=0,yend=-0.02,color=keep))+
         scale_x_log10()
)
genes_summary.dt[,.N,by=keep]


cells_summary.dt <-
  counts.dt[,.(total_counts=sum(count),genes_with_counts=.N), by=cell_id
            ][, keep := !isOutlier(total_counts, nmads=2, type="lower", log=TRUE)]
cells_summary.dt[,CDR:=genes_with_counts/n_genes]
cells_summary.dt[,total_counts_scaled:=total_counts/max(total_counts)]


figure("Distribution of counts over cells",
       ggplot(cells_summary.dt, aes(x=total_counts))+
         geom_density()+scale_x_log10(minor_breaks= (1:1000)*5E3)+
         geom_segment(aes(x=total_counts, xend=total_counts, y=0, yend=-0.2, color=keep))
)
cells_summary.dt[,.N,by=keep]


# filtering
counts.dt <- counts.dt[gene_id %in% genes_summary.dt[keep==TRUE,gene_id]]
counts.dt <- counts.dt[cell_id %in% cells_summary.dt[keep==TRUE,cell_id]]


counts.dt[, gene_id2:=.GRP, keyby=gene_id]
counts.dt[, cell_id2:=.GRP, keyby=cell_id]


count_matrix <- counts.dt[,as.matrix(sparseMatrix(cell_id2, gene_id2, x=count))]

genenames.dt <- genenames.dt[unique(counts.dt, by=c("gene_id","gene_id2")), on="gene_id"]
cellnames.dt <- cellnames.dt[unique(counts.dt, by=c("cell_id","cell_id2")), on="cell_id"]

colnames(count_matrix) <- setorder(genenames.dt,gene_id2)[, gene]
rownames(count_matrix) <- setorder(cellnames.dt,cell_id2)[, cell]

# covariate matrix -----------------

covariates.dt <- cbc_gbc_dict.dt[cells_summary.dt[cellnames.dt, on="cell_id"], on="cell"]
covariates.dt[is.na(guide),guide:="none"]
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide
X[c(1:3,(nrow(X)-3):nrow(X)),]


# expression matrix
transformation_method = "Log(1+counts)" # "Log(1+counts)" #"Ancombes" #"Log(1+counts)"

Y <- switch(transformation_method,
            Ancombes=stabilize_Anscombes(count_matrix),
            `Log(1+counts)`=log2(1+count_matrix)
)
Y[1:3,1:4]

adjustment_terms <- colnames(X)


fit <- lm(Y ~ .,as.data.table(X))# lm(x=X,y=Y)
anova(fit)
Y_adj <- residuals(fit)

fit <- lm(Y ~ .,as.data.table(X))# lm(x=X,y=Y)
anova(fit)

library(viridis)
gene_gene <- cor(Y)
gene_gene_adj <- cor(Y_adj)
if (!exists("new_order")){
  new_order <- (hclust(dist(gene_gene)))$order
}

figure1(paste(transformation_method,"raw"),
image(gene_gene[new_order,][,new_order],zlim=c(-1,1),col=viridis(100), axes=F, xlab="gene", ylab="gene"))

figure1(paste(transformation_method, "adjusted"),
image(gene_gene_adj[new_order,][,new_order],zlim=c(-1,1),col=viridis(100), axes=F, xlab="gene", ylab="gene"))



library(Hmisc)
figure1(paste(transformation_method, " adjusted p-val distribution"),
hist(rcorr(Y_adj, type="pearson")$P))

figure1(paste(transformation_method, " adjusted null p-val distribution"),
hist(rcorr(apply(Y_adj, 2, sample) , type="pearson")$P))

figure1(paste(transformation_method, " raw p-val distribution"),
hist(rcorr(Y, type="pearson")$P))

figure1(paste(transformation_method, " raw null p-val distribution"),
hist(rcorr(apply(Y, 2, sample) , type="pearson")$P))

figure1(paste(transformation_method, "QQ plot"),
qqplot(rcorr(Y_adj, type="pearson")$P, rcorr(Y, type="pearson")$P))

figure1(paste(transformation_method, "QQ plot null"),
qqplot(rcorr(apply(Y_adj, 2, sample), type="pearson")$P, rcorr(apply(Y, 2, sample), type="pearson")$P))
