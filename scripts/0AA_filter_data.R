# inspection --------------------
logfile= file.path("results",paste0("filter_out_",data_set,".log"))
cat("",file=logfile)
logprint <-  function(x) {
  xx <-  substitute(x)
  res <- eval(xx)
  cat(paste0("\"",deparse(xx),"\""), res,"\n", file=logfile, append=TRUE)
  res
}


true_false_scale  <- {
  dd <- c(TRUE, FALSE)
  dd.col <- cool_warm(2)
  names(dd.col)  <- dd
  dd.lab <- c("yes","no")
  names(dd.lab)  <- dd
  scale_color_manual(name="kept",values = dd.col, labels=dd.lab)
}


guides_per_cell_dt <- unique(cbc_gbc_dict.dt,by=c("cell","guide"))[,.(unique_guides_per_cell=sum(!is.na(guide))),by=cell]

figure(paste0("guides per cell – ", data_set),
       ggplot(guides_per_cell_dt,aes(x=as.factor(unique_guides_per_cell)))+geom_bar(aes(y=..count../sum(..count..)))+
         xlab("# unique guides per cell (detected)")+ylab("relative frequency"), height=3,width=5
)

guides_summary <- cbc_gbc_dict.dt[,.N,by=guide]
guides_summary[,keep:=N>30]
cbc_gbc_dict.dt[guides_summary, keep:=keep, on="guide"]
guides_summary[,.N,by=keep]

sscale_log10 <- function(x) {
  upper <-  floor(log10(max(x)))
  lower <- ceiling(log10(min(x)))
  if (upper-lower+1>10) {
    scales::trans_breaks("log10", function(x) 10^x)
  } else if (upper-lower>0) {
    scales::trans_breaks("log10", function(x) 10^x, n = upper-lower+1)
  } else {
    c( 10^upper*(1:9)) #10^(upper-1)* (1:9),
  }
}
#
figure(paste0("distribution of number of cells with guide – ", data_set),
       ggplot(guides_summary, aes(x=N))+
         geom_density(color=NA, fill="gray", aes(y=..scaled..)) +
         ylab("scaled density estimate") +
         geom_segment(aes(x=N,xend=N,y=0,yend=0.1,color=keep))+
         xlab('# cells with respective guide')+ expand_limits(x=0) +
         scale_x_log10(breaks = sscale_log10(guides_summary$N)) +
         annotation_logticks(sides="b") +
         true_false_scale +
         theme(legend.justification = c(1, 1), legend.position = c(1, 1)),
       height=3,width=8
)

NUM_GENES_TO_KEEP = 2000

n_cells <- length(unique(counts.dt[,cell_id]))
genes_summary.dt <-
  counts.dt[, .(total_counts=sum(count), cells_with_counts=.N),by=.(gene_id)
            ][, `:=`(keep = rank(-total_counts)<=NUM_GENES_TO_KEEP,
                     frac_cells_with_counts=cells_with_counts/n_cells,
                     mean_counts_per_cell_with_counts=total_counts/cells_with_counts,
                     mean_count=total_counts/n_cells)]

figure(paste0("Distribution fraction of cells with counts over genes – ", data_set),
       ggplot(genes_summary.dt, aes(x=frac_cells_with_counts)) +
         geom_density(color=NA, fill="gray", aes(y=..scaled..)) +
         ylab("scaled density estimate") +
         geom_segment(aes(x=frac_cells_with_counts,xend=frac_cells_with_counts,y=0,yend=0.1,color=keep))+
         scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                       limits=c(genes_summary.dt[,min(frac_cells_with_counts)],1))  +
         annotation_logticks(sides="b") +
         xlab("fraction of cells with at least one count for the gene") +
         true_false_scale +
         theme(legend.justification = c(0, 1), legend.position = c(0.01, 1)),
       height=3, width=8
)
logprint(median(genes_summary.dt[,frac_cells_with_counts]))

figure(paste0("Distribution fraction of cells with counts over kept genes – ", data_set),
       ggplot(genes_summary.dt[keep==TRUE], aes(x=frac_cells_with_counts))+
         geom_density(color=NA, fill="gray", aes(y=..scaled..)) +
         ylab("scaled density estimate") +
         geom_segment(aes(x=frac_cells_with_counts,xend=frac_cells_with_counts,y=0,yend=0.1,color=keep))+
         expand_limits(y=1) +
         xlab("fraction of cells with at least one count for the gene") +
         true_false_scale +
         theme(legend.justification = c(0, 1), legend.position = c(0.01, 1)),
       height=3, width=8
)
logprint(median(genes_summary.dt[keep==TRUE,frac_cells_with_counts]))

figure(paste0("Distribution fraction of cells with counts over genes (zoomed) – ", data_set),
       ggplot(genes_summary.dt, aes(x=frac_cells_with_counts))+
         geom_density(color=NA, fill="gray", aes(y=..scaled..*35)) +
         ylab("scaled density estimate") +
         geom_segment(aes(x=frac_cells_with_counts,xend=frac_cells_with_counts,y=0,yend=0.1,color=keep))+
         expand_limits(y=1) +
         xlab("fraction of cells with at least one count for the gene") +
         true_false_scale +
         coord_cartesian(xlim=genes_summary.dt[keep==TRUE,range(frac_cells_with_counts)],ylim=c(0,1)) +
         theme(legend.justification = c(0.5, 1), legend.position = c(0.5, 1)),
       height=3, width=8
)

figure(paste0("Distribution mean counts per cells with counts over genes – ", data_set),
       ggplot(genes_summary.dt, aes(x=mean_counts_per_cell_with_counts-1))+
         geom_density()+
         scale_x_log10()
)
median(genes_summary.dt[,mean_counts_per_cell_with_counts])


figure(paste0("Distribution mean counts over genes – ", data_set),
       ggplot(genes_summary.dt, aes(x=mean_count))+
         geom_density(color=NA, fill="gray", aes(y=..scaled..)) +
         ylab("scaled density estimate") +
         geom_segment(aes(x=mean_count,xend=mean_count,y=0,yend=0.1,color=keep))+
         scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x)) +
         annotation_logticks(sides="b") +
         xlab("mean count per cell of gene") +
         true_false_scale +
         theme(legend.justification = c(1, 1), legend.position = c(1, 1)),
       height=3, width=8
)
logprint(median(genes_summary.dt[,mean_count]))

figure(paste0("Distribution counts over kept genes – ", data_set),
       ggplot(genes_summary.dt[keep==TRUE], aes(x=mean_count))+
         geom_density(color=NA, fill="gray", aes(y=..scaled..)) +
         ylab("scaled density estimate") +
         geom_segment(aes(x=mean_count,xend=mean_count,y=0,yend=0.1,color=keep))+
         scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x,n=3)) +
         annotation_logticks(sides="b") +
         xlab("mean count per cell of kept gene") +
         ylab("scaled density estimate")+
         true_false_scale, height=3, width=8
)
logprint(median(genes_summary.dt[keep==TRUE,mean_count]))


genes_summary.dt[,.N,by=keep]

cells_summary.dt <-
  counts.dt[,.(total_counts=sum(count),genes_with_counts=.N), by=cell_id
            ][, keep := total_counts>quantile(total_counts, 50/100)]

cells_summary.dt[cbc_gbc_dict.dt[keep==FALSE], keep:=FALSE,on="cell_id"]
cells_summary.dt[,n_guides:=0L]
cells_summary.dt[cbc_gbc_dict.dt[,.(n_guides=.N),by=cell_id], n_guides:=i.n_guides, on="cell_id"]

# cells_summary.dt[n_guides==0, keep := FALSE]

cells_summary.dt[,CDR:=genes_with_counts/n_genes]
cells_summary.dt[,total_counts_scaled:=total_counts/max(total_counts)]


figure(paste0("Distribution of counts over cells – ", data_set),
     ggplot(cells_summary.dt, aes(x=total_counts))+
       geom_density(color=NA, fill="gray", aes(y=..scaled..)) +
       ylab("scaled density estimate") +
       geom_segment(aes(x=total_counts, xend=total_counts, y=0, yend=0.1, color=keep))+
       scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x,n=2)) +
       annotation_logticks(sides="b") +
       xlab("total counts of cell") +
       true_false_scale +
       theme(legend.justification = c(1, 1), legend.position = c(1, 1)),
     height=3, width=8
)
logprint(median(cells_summary.dt[,total_counts]))

cells_summary.dt[,.N,by=keep]


# filtering ------
counts.dt <- counts.dt[gene_id %in% genes_summary.dt[keep==TRUE,gene_id]]
counts.dt <- counts.dt[cell_id %in% cells_summary.dt[keep==TRUE,cell_id]]
cbc_gbc_dict.dt <- cbc_gbc_dict.dt[keep==TRUE]


counts.dt[, gene_id2:=.GRP, keyby=gene_id]
counts.dt[, cell_id2:=.GRP, keyby=cell_id]

# count matrix  ---
count_matrix <- counts.dt[,as.matrix(Matrix::sparseMatrix(cell_id2, gene_id2, x=count))]

genenames.dt <- genenames.dt[unique(counts.dt, by=c("gene_id","gene_id2")), on="gene_id"]
cellnames.dt <- cellnames.dt[unique(counts.dt, by=c("cell_id","cell_id2")), on="cell_id"]

colnames(count_matrix) <- setorder(genenames.dt,gene_id2)[, gene]
rownames(count_matrix) <- setorder(cellnames.dt,cell_id2)[, cell]
# guide matrix  ---
cbc_gbc_dict.dt[, guide_id:=.GRP, by=guide]
guide_matrix <- cbc_gbc_dict.dt[cellnames.dt, as.matrix(Matrix::sparseMatrix(cell_id2, guide_id, x=TRUE)), on="cell_id", nomatch=0]
colnames(guide_matrix) <- cbc_gbc_dict.dt[,.(guide=guide[1]), keyby=guide_id][, guide]
rownames(guide_matrix) <- rownames(count_matrix)
# covariate matrix ---
covariates.dt <- cbc_gbc_dict.dt[,.(guides=paste0(unique(guide),collapse=" & "), n_guides=.N), by=cell_id][cells_summary.dt[keep==TRUE,],on="cell_id"]
covariates.dt[is.na(guides), `:=`(guides="", n_guides=0)]
setorder(covariates.dt,"cell_id")
covariates.dt <- cellnames.dt[,.(cell_id,cell,batch,cell_id2)][covariates.dt,on="cell_id"]
# batch matrix ---
batch_matrix <- model.matrix(~0+batch,covariates.dt)
rownames(batch_matrix) <- rownames(count_matrix)

# export --------------
data_folder <- paste0('data_processed/', data_set)
dir.create(data_folder, showWarnings = FALSE)

fwrite(covariates.dt, file.path(data_folder, "covariates.dt.csv"))

save(count_matrix, file= file.path(data_folder, "count_matrix.RData"))
save(guide_matrix, file= file.path(data_folder, "guide_matrix.RData"))
save(batch_matrix, file= file.path(data_folder, "batch_matrix.RData"))
fwrite(as.data.table(count_matrix, keep.rownames = TRUE), file = file.path(data_folder, "count_matrix.csv"))
fwrite(as.data.table(guide_matrix, keep.rownames = TRUE), file = file.path(data_folder, "guide_matrix.csv"))
fwrite(as.data.table(batch_matrix, keep.rownames = TRUE), file = file.path(data_folder, "batch_matrix.csv"))
feather::write_feather(as.data.table(count_matrix, keep.rownames = TRUE),  file.path(data_folder, "count_matrix.feather"))
feather::write_feather(as.data.table(guide_matrix, keep.rownames = TRUE),  file.path(data_folder, "guide_matrix.feather"))
feather::write_feather(as.data.table(batch_matrix, keep.rownames = TRUE),  file.path(data_folder, "batch_matrix.feather"))
