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
