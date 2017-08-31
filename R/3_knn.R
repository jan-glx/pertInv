
load(file="results/count_matrix.RData")
load(file="results/guide_matrix.RData")

library(FNN)
dim(count_matrix)

one10 <- FNN::get.knnx(count_matrix[1:10000,],count_matrix[1,, drop=F], k=50)$nn.index

ss <- count_matrix[one10,]
library(fitdistrplus)
fitg <- fitdist(ss[,order(colSums(ss))[ sample(ncol(ss),1)]],"pois")
summary(fitg)
plot(fitg)
