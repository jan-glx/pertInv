
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

data_set = "GSM2396857_dc_0hr"
library(rhdf5)
rhdf5::h5createFile("results/data.h5")
file <-  H5Fopen("results/data.h5")
h5createGroup(file, data_set)
save_matrix <- function(mat) {
  name <- deparse(substitute(mat))
  h5createGroup(file, paste(data_set, name, sep="/"))
  h5write(mat, file, paste(data_set, name, "values", sep="/"))
  h5write(rownames(mat), file, paste(data_set, name, "rownames", sep="/"))
  h5write(colnames(mat), file, paste(data_set, name, "colnames", sep="/"))
}
storage.mode(count_matrix) <- "integer"
storage.mode(guide_matrix) <- "integer"
save_matrix(guide_matrix)
save_matrix(count_matrix)

H5close()


a=h5read("results/data.h5", paste(data_set, "count_matrix", "values", sep="/"))
