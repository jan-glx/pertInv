devtools::install('.')
source("scripts/1_load_data.R")
source("scripts/2_filter_data.R")

args <- commandArgs(trailingOnly = TRUE)
ii <- as.integer(args[1])
source("scripts/3_ICP.R")
