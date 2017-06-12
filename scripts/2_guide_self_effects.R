
genenames.dt[!is.na(targeted_by)]

counts.dt <- counts.dt[genenames.dt[!is.na(targeted_by)], on="gene_id"]

cells_summary.dt <-
  counts.dt[,.(total_counts=sum(count),genes_with_counts=.N), by=cell_id
            ][, keep := !scater::isOutlier(total_counts, nmads=2, type="lower", log=TRUE)]
cells_summary.dt[,n_guides:=0L]
cells_summary.dt[cbc_gbc_dict.dt[,.(n_guides=.N),by=cell_id], n_guides:=i.n_guides, on="cell_id"]

#cells_summary.dt[n_guides>1, keep := FALSE]

cells_summary.dt[,CDR:=genes_with_counts/n_genes]
cells_summary.dt[,total_counts_scaled:=total_counts/max(total_counts)]


cells_summary.dt[,.N,by=keep]


# filtering
counts.dt <- counts.dt[cell_id %in% cells_summary.dt[keep==TRUE,cell_id]]


counts.dt[, gene_id2:=.GRP, keyby=gene_id]
counts.dt[, cell_id2:=.GRP, keyby=cell_id]


count_matrix <- counts.dt[,as.matrix(Matrix::sparseMatrix(cell_id2, gene_id2, x=count))]

genenames.dt <- genenames.dt[unique(counts.dt, by=c("gene_id","gene_id2")), on="gene_id"]
cellnames.dt <- cellnames.dt[unique(counts.dt, by=c("cell_id","cell_id2")), on="cell_id"]

colnames(count_matrix) <- setorder(genenames.dt,gene_id2)[, gene]
rownames(count_matrix) <- setorder(cellnames.dt,cell_id2)[, cell]



cbc_gbc_dict.dt[,guide_id := as.integer(factor(guide))]

guide_matrix <- cbc_gbc_dict.dt[cells_summary.dt[cellnames.dt, on="cell_id"], on="cell"][!is.na(guide_id), as.matrix(Matrix::sparseMatrix(cell_id2, guide_id, x=TRUE))]

covariates.dt <- cells_summary.dt[cellnames.dt, on="cell_id"]
Y <- stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)


for (null in c("", "null (permuted targeted labels)")){

  dt <- melt(data.table(Y_adj, keep.rownames = TRUE), id.vars="rn")
  setnames(dt, c("cell", "gene", "adjusted expression"))
  dt[,c("ENSG", "gene short") :=tstrsplit(gene,split = "_")]
  dt[, targeted:=FALSE]
  dt[cbc_gbc_dict.dt,targeted := TRUE, on=c("cell","gene short"="target_gene")]

  if (null=="null (permuted targeted labels)") dt[, targeted:=sample(targeted)]

  figure(paste0("box-plot self effect of guides", null),
         ggplot(dt, aes(x=gene,y=`adjusted expression`, color=targeted))+geom_boxplot()+coord_flip()
  )

  dt[,`:=`(k=rank(`adjusted expression`),n=.N), by=.(gene,targeted)]
  dt[,p_self:=k/n]
  dt[,`:=`(p=(rank(`adjusted expression`)-k)/(.N-n)), by=.(gene)]
  dt[(targeted),`:=`(p_targeted=p_self,p_non_targeted=p)]
  dt[(!targeted),`:=`(p_targeted=p,p_non_targeted=p_self)]

  figure(paste0("pp-plot self effect of guides", null),
         ggplot(dt, aes(x=p_non_targeted,y=p_targeted, color=gene))+geom_line()+geom_abline()
  )

  res <- dt[,.(p.value=wilcox.test(`adjusted expression`[targeted],`adjusted expression`[!targeted],alternative="less")$p.value),by=gene]
  print(res[,p.adjust(p.value)])
}
