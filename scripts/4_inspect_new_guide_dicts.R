# ----------------
data_folder <- 'data_raw/'
data_set = "GSM2396857_dc_0hr"
# promoters_concat_all.csv >=  GSM2396858_k562_tfs_7 #subsampling of guide_cell mappings
# pt2_concat_all.csv >= GSM2396859_k562_tfs_13 # subsampling of cells with exactly one guide
# dc_0hr_concat_all.csv == GSM2396857_dc_0hr_strict
# ph_concat_all.csv == GSM2396860_k562_tfs_highmoi_lenient
cbc_gbc_dict.dt <- fread(paste0('gzip -dc ', data_folder, data_set,'_cbc_gbc_dict_new.csv.gz'), header= FALSE)
setnames(cbc_gbc_dict.dt, c("guide", "cells"))
cbc_gbc_dict.dt[, target_gene:=stringr::str_match(guide,"^p_(?:sg)?((?:(?<=sg).*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]]

cbc_gbc_dict.dt <- cbc_gbc_dict.dt[,.(cell=stringr::str_split(cells, ", ")[[1]]), by=.(guide, target_gene)]
cbc_gbc_dict.dt[,.N,by=cell][,.N,keyby=N]

#ph_concat_all.csv  promoters_concat_all.csv  pt2_concat_all.csv  dc_0hr_concat_all.csv

cbc_gbc_dict.du <-    fread(paste0(data_folder, 'dc_0hr_concat_all.csv'), header= FALSE) #
setnames(cbc_gbc_dict.du, c("guide", "cells"))
cbc_gbc_dict.du[, target_gene:=stringr::str_match(guide,"^p_(?:sg)?((?:(?<=sg).*(?=_))|(?:INTERGENIC))(?:_)?\\d+$")[,2]]

cbc_gbc_dict.du <- cbc_gbc_dict.du[,.(cell=stringr::str_split(cells, ", ")[[1]]), by=.(guide, target_gene)]
cbc_gbc_dict.du[,.N,by=cell][,.N,keyby=N]

cbc_gbc_dict.dt[,.N,by=cell][,.N,keyby=N]
cbc_gbc_dict.du[cbc_gbc_dict.dt,on="cell", nomatch=0][,.N,by=cell][,.N,keyby=N]

cbc_gbc_dict.dt[cbc_gbc_dict.du,on="cell", nomatch=0][,.N,by=cell][,.N,keyby=N]
cbc_gbc_dict.du[!cbc_gbc_dict.dt,on="cell"][,.N,by=cell][,.N,keyby=N]
cbc_gbc_dict.dt[!cbc_gbc_dict.du,on="cell"][,.N,by=cell][,.N,keyby=N]
