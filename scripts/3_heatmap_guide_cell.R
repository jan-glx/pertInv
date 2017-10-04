guide_vector <- guide_matrix %*% seq_len(ncol(guide_matrix))

ordering <- order(guide_vector)

lattice::levelplot(Y_adj[ordering[ordering<1000],], scales=list(draw=FALSE),
                   xlab="gene", ylab="cell", col.regions = viridis::viridis(500),
                   at=seq(-1,1,length.out=500))
image(Y_adj[ordering[ordering<1000],1:100], col=viridis::viridis(500))
image(matrix(rep((guide_vector[ordering[ordering<1000]]),100),ncol=100), col=rainbow(max(ordering)))




SS <- Y_tmp[ordering[ordering<1000],1:1000]
annotation_row <- data.frame(row.names = rownames(SS),guide = guide_vector[ordering[ordering<1000]])
pheatmap::pheatmap(SS, cluster_rows = T, cluster_cols = T, annotation_row=annotation_row,
                   show_colnames = F, show_rownames = F)#,annotation_colors=list(guide=sample(rainbow(max(guide_vector))))

