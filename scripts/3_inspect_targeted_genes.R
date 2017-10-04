


genenames.dt[!is.na(targeted_by)]
counts.dt[,.N,by=gene_id][genenames.dt[!is.na(targeted_by)], on="gene_id"][,perc_detected:=100*N/nrow(cellnames.dt)][]

counts.dt[]
genenames.dt[!is.na(targeted_by)][cbc_gbc_dict.dt,on="gene_id"]

cbc_gbc_dict.dt[genenames.dt,on=c("target_gene"="targeted_by"),gene_id:=gene_id,]
counts.dt[cbc_gbc_dict.dt, on=c("gene_id","cell_id"),nomatch=0][,.N,by=c("gene_id","target_gene")][cbc_gbc_dict.dt[,.N,by="gene_id"], on="gene_id"][,perc_detected:=100*N/i.N][]
