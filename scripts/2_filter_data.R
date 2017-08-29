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
            ][, `:=`(keep = rank(-cells_with_counts)<=NUM_GENES_TO_KEEP,
                     frac_cells_with_counts=total_counts/.N,
                     mean_counts_per_cell_with_counts=total_counts/cells_with_counts)]

figure("Distribution fraction of cells with counts over genes",
       ggplot(genes_summary.dt, aes(x=frac_cells_with_counts))+
         geom_density()+
         scale_x_log10()
)
median( genes_summary.dt[,frac_cells_with_counts])

figure("Distribution mean counts per cells with counts over genes",
       ggplot(genes_summary.dt, aes(x=mean_counts_per_cell_with_counts-1))+
         geom_density()+
         scale_x_log10()
)
median( genes_summary.dt[,mean_counts_per_cell_with_counts])


figure("Distribution counts over gene",
       ggplot(genes_summary.dt, aes(x=total_counts))+
         geom_density()+
         geom_segment(aes(x=total_counts,xend=total_counts,y=0,yend=-0.02,color=keep))+
         scale_x_log10()
)
genes_summary.dt[,.N,by=keep]


cells_summary.dt <-
  counts.dt[,.(total_counts=sum(count),genes_with_counts=.N), by=cell_id
            ][, keep := total_counts>quantile(total_counts, 50/100)]
cells_summary.dt[,n_guides:=0L]
cells_summary.dt[cbc_gbc_dict.dt[,.(n_guides=.N),by=cell_id], n_guides:=i.n_guides, on="cell_id"]

cells_summary.dt[n_guides>1, keep := FALSE]

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


count_matrix <- counts.dt[,as.matrix(Matrix::sparseMatrix(cell_id2, gene_id2, x=count))]

genenames.dt <- genenames.dt[unique(counts.dt, by=c("gene_id","gene_id2")), on="gene_id"]
cellnames.dt <- cellnames.dt[unique(counts.dt, by=c("cell_id","cell_id2")), on="cell_id"]

colnames(count_matrix) <- setorder(genenames.dt,gene_id2)[, gene]
rownames(count_matrix) <- setorder(cellnames.dt,cell_id2)[, cell]
# -------------
# covariate matrix -----------------
cbc_gbc_dict.dt[,guide_id := as.integer(factor(guide))]
covariates.dt <- cbc_gbc_dict.dt[cells_summary.dt[cellnames.dt, on="cell_id"], on="cell"]

guide_matrix <- covariates.dt[!is.na(guide_id), as.matrix(Matrix::sparseMatrix(cell_id2, guide_id, x=TRUE))]

covariates.dt[is.na(guide),guide:="none"]

fwrite(covariates.dt,"results/covariates.dt.csv")
save(count_matrix,file="results/count_matrix.RData")
