

genes_summary.dt <-
  counts.dt[, .(total_counts=sum(count), cells_with_counts=.N),by=.(gene_id)
            ][, p := rank(cells_with_counts)/.N]


genes_summary.dt[genenames.dt[!is.na(targeted_by)], on="gene_id"]

counts.dt2 <- counts.dt[genenames.dt[!is.na(targeted_by)], on="gene_id"]

cells_summary.dt <-
  counts.dt2[,.(total_counts=sum(count),genes_with_counts=.N), by=cell_id
            ][, keep := !scater::isOutlier(total_counts, nmads=2, type="lower", log=TRUE)]
cells_summary.dt[,n_guides:=0L]
cells_summary.dt[cbc_gbc_dict.dt[,.(n_guides=.N),by=cell_id], n_guides:=i.n_guides, on="cell_id"]

#cells_summary.dt[n_guides>1, keep := FALSE]

cells_summary.dt[,CDR:=genes_with_counts/n_genes]
cells_summary.dt[,total_counts_scaled:=total_counts/max(total_counts)]


cells_summary.dt[,.N,by=keep]


# filtering
counts.dt2 <- counts.dt2[cell_id %in% cells_summary.dt[keep==TRUE,cell_id]]


counts.dt2[, gene_id2:=.GRP, keyby=gene_id]
counts.dt2[, cell_id2:=.GRP, keyby=cell_id]


count_matrix <- counts.dt2[,as.matrix(Matrix::sparseMatrix(cell_id2, gene_id2, x=count))]

genenames.dt <- genenames.dt[unique(counts.dt2, by=c("gene_id","gene_id2")), on="gene_id"]
cellnames.dt <- cellnames.dt[unique(counts.dt2, by=c("cell_id","cell_id2")), on="cell_id"]

colnames(count_matrix) <- setorder(genenames.dt,gene_id2)[, gene]
rownames(count_matrix) <- setorder(cellnames.dt,cell_id2)[, cell]



cbc_gbc_dict.dt[,guide_id := as.integer(factor(guide))]

guide_matrix <- cbc_gbc_dict.dt[cells_summary.dt[cellnames.dt, on="cell_id"], on="cell"][!is.na(guide_id), as.matrix(Matrix::sparseMatrix(cell_id2, guide_id, x=TRUE))]

covariates.dt <- cells_summary.dt[cellnames.dt, on="cell_id"]
Y <-  sweep(count_matrix, 2, colMeans(count_matrix), "/") #stabilize_Anscombes(count_matrix) #count_matrix
X <- model.matrix(~0+I(CDR^2)+I(CDR^3)+CDR+total_counts_scaled+batch,data=covariates.dt) #+guide

fit <- lm(Y ~ .,as.data.table(X))
Y_adj <- residuals(fit)

null <- ""
for (null in c("", "null (permuted targeted labels)")){

  dt <- melt(data.table(Y_adj, keep.rownames = TRUE), id.vars="rn")
  setnames(dt, c("cell", "gene", "adjusted expression"))
  dt[,c("ENSG", "gene_short") :=tstrsplit(gene,split = "_")]

  dt <- cbc_gbc_dict.dt[dt, on = c("cell","target_gene"="gene_short")]

  dt[, targeted:=!is.na(guide)]

  if (null=="null (permuted targeted labels)") dt[, targeted:=sample(targeted)]

  figure(paste0("box-plot self effect of guides", null),
         ggplot(dt, aes(x=gene,y=`adjusted expression`, color=guide))+geom_boxplot()+coord_flip()+
           geom_hline(yintercept=0, linetype="dashed"),width=10
  )

  dt[,`:=`(k=rank(`adjusted expression`),n=.N), by=.(gene,guide)]
  dt[,p_self:=k/n]
  dt[,`:=`(k2=rank(`adjusted expression`),n2=.N), by=.(gene,targeted)]
  dt[,p_self2:=k2/n2]
  dt[,`:=`(p=(rank(`adjusted expression`)-k2)/(.N-n2)), by=.(gene)]
  dt[(targeted),`:=`(p_targeted=p_self,p_non_targeted=p)]
  dt[(!targeted),`:=`(p_targeted=p,p_non_targeted=p_self)]

  figure(paste0("pp-plot self effect of guides", null),
         ggplot(setorder(dt[(targeted)],p_non_targeted), aes(x=p_non_targeted,y=p_targeted, color=guide))+geom_line()+geom_abline()
  )


  res <- dt[, {
    null_sample<-`adjusted expression`[(!targeted)]
    .SD[(targeted),
        .(p.value=wilcox.test(`adjusted expression`, null_sample, alternative="less")$p.value),
        by=guide]
        },by=gene]
  res[,p.value.adj:=p.adjust(p.value, method="BH")]
  print(res)

  pvalues <- res[,log(p.value)]
  figure1(paste0("hist log p.values", null),
          hist(pvalues),
          sub_title = TRUE)
}
