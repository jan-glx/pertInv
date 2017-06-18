library(pertInv)

# Load count data -----------------------------------
data_folder <- 'data_raw/'
data_set = "GSM2396859_k562_tfs_13"
# "GSM2396861_k562_ccycle"
# "GSM2396858_k562_tfs_7"
# "GSM2396859_k562_tfs_13"
# "GSM2396860_k562_tfs_highmoi"
# "GSM2396856_dc_3hr"
# "GSM2396857_dc_0hr"
# "GSM2396859_k562_tfs_13"

counts.dt <- fread(paste0('gzip -dc ', data_folder, data_set,'.mtx.txt.gz'))
setnames(counts.dt,c("gene_id","cell_id","count"))
n_genes <- counts.dt[1,gene_id]
n_cells <- counts.dt[1,cell_id]
n_not_no_reads <- counts.dt[1,count]
counts.dt <- counts.dt[-1]

# Load ID-cell and ID-gene mapping
cellnames.dt <- fread(paste0('gzip -dc ', data_folder, data_set, '_cellnames.csv.gz'), header=TRUE)
setnames(cellnames.dt, c("cell_id", "cell"))
cellnames.dt[,cell_id:=cell_id+1] # changing to one_based

cellnames.dt[,batch:=stringr::str_match(cell,".*_(.*_.*)")[,2]]
cellnames.dt[,cell_bc:=stringr::str_match(cell,"(.*)_.*_.*")[,2]]

genenames.dt <- fread(paste0('gzip -dc ', data_folder, data_set, '_genenames.csv.gz'), header=TRUE)
setnames(genenames.dt, c("gene_id", "gene"))
genenames.dt[,gene_id:=gene_id+1] # changing to one_based



# Load cell-guide mapping
cbc_gbc_dict.dt <- fread(paste0('gzip -dc ', data_folder, data_set, '_cbc_gbc_dict_new.csv.gz'), header= FALSE)
setnames(cbc_gbc_dict.dt, c("guide", "cells"))
cbc_gbc_dict.dt[, target_gene:=stringr::str_match(guide,"^p_(?:sg)?((?:(?<=sg).*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]]

cbc_gbc_dict.dt <- cbc_gbc_dict.dt[,.(cell=stringr::str_split(cells, ", ")[[1]]), by=.(guide, target_gene)]
cbc_gbc_dict.dt[,.N,by=cell][,.N,keyby=N]
cbc_gbc_dict.dt[cellnames.dt, cell_id:=cell_id, on="cell"]
## add batch info and subset observed cells
#cbc_gbc_dict.dt <- cbc_gbc_dict.dt[cellnames.dt,on="cell"]

# find targt genes of guides
genenames.dt[, targeted_by:=stringr::str_match(
  gene,
  paste0("_(", paste0(cbc_gbc_dict.dt[,unique(target_gene)],collapse="|"), ")$")
)[,2]]

