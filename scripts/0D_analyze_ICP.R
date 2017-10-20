
files <- list.files(path="results", pattern = "ICP_run_\\d*\\.Rdata",full.names = TRUE)
files <- files[order(as.integer(stringr::str_match(files,"ICP_run_(\\d*)\\.Rdata")[,2]))]

all_res <- unlist(lapply(files, function(file) {load(file=file); res}), recursive = FALSE)
#size of parent set at wich intersection becomes empty
table(sapply(all_res, function(x) x$stop_at_cardinality), useNA="always")

# Variable selection (up to 8)
table(sapply(all_res, function(x) length(x$usedvariables)))

hist(sapply(all_res, function(x) length(x$usedvariables)))


# Variable selection (up to 8)
model_pvals <- sapply(all_res, function(x) (x$bestModel))
hist(model_pvals)

hist(sapply(all_res[is.na(sapply(all_res, function(x) x$stop_at_cardinality))], function(x) (x$bestModel)))


# Variable selection (up to 8)
hist(sapply(all_res, function(x) sum(x$pvalues<0.05)))

cardinality_S_hat <- as.integer((sapply(all_res, function(x) sum(x$pvalues<0.1))))

library(pertInv)

figure("Histogram of cardinality of S_hat",
ggplot(data.table(cardinality=cardinality_S_hat), aes(x=cardinality))+stat_count()
)

figure("Histogram of cardinality of S_hat (zoomed)",
ggplot(data.table(cardinality=cardinality_S_hat)[cardinality>0,], aes(x=cardinality))+stat_count()
)



figure("#selected vs #preselected",
ggplot(data.table(
  n_preselected=sapply(all_res, function(x) length(x$usedvariables)),
  n_selected= as.integer((sapply(all_res, function(x) sum(x$pvalues<0.1)))),
  model_pval =  sapply(all_res, function(x) x$bestModel)),
  aes(x=n_preselected,n_selected,color=-log10(model_pval)))+geom_jitter(size=3,shape=16)+
  viridis::scale_color_viridis(lim=c(0,2),name="neg. log10\nmodel p.value\n(grey: p<0.01)")+
  theme(legend.position = c(0.2, 0.8))
)

all_res.dt <- data.table(
  n_preselected=sapply(all_res, function(x) length(x$usedvariables)),
  n_selected= as.integer((sapply(all_res, function(x) sum(x$pvalues<0.1)))),
  model_pval =  sapply(all_res, function(x) x$bestModel))
