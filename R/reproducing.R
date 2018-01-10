
#' Download and extract data files from source repositories.
#' @export
setup_data <- function() {
  dir.create("data_raw", showWarnings = FALSE)
  base_path <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90063/suppl/"
  download_data_file <- function(base_path, filename) utils::download.file(paste0(base_path, filename), file.path("data_raw", filename), method = "auto", mode="wb")

  if (!file.exists(file.path("data_raw", "GSM2396858_k562_tfs_7.mtx.txt.gz"))){
    download_data_file(base_path, "GSE90063_RAW.tar")
    utils::untar(file.path("data_raw", "GSE90063_RAW.tar"), exdir="data_raw")
    file.remove(file.path("data_raw", "GSE90063_RAW.tar"))
  }

  download_data_file(base_path, "GSE90063_dc0hr_umi_wt.txt.gz")
  download_data_file(base_path, "GSE90063_dc3hr_umi_wt.txt.gz")
  download_data_file(base_path, "GSE90063_k562_umi_wt.txt.gz")

  base_path <- "https://raw.githubusercontent.com/asncd/MIMOSCA/master/GBC_CBC_pairing/gbc_cbc_dicts/"
  download_and_compress_data_file <- function(base_path, filename, new_filename) {
    in_connection <- file(paste0(base_path, filename), open= "rb", raw=TRUE)
    out_connection <- gzfile(file.path("data_raw", paste0(new_filename, ".gz")), "wb")
    while(length(buf <- readBin(in_connection, "raw", n = 1024*1024))) {
      writeBin(buf, out_connection)
    }
    close(out_connection)
    close(in_connection)
  }

  download_and_compress_data_file(base_path, "dc_0hr_concat_all.csv", "GSM2396857_dc_0hr_cbc_gbc_dict_new.csv")
  download_and_compress_data_file(base_path, "ph_concat_all.csv", "GSM2396860_k562_tfs_highmoi_lenient_cbc_gbc_dict_new.csv")
  download_and_compress_data_file(base_path, "promoters_concat_all.csv", "GSM2396858_k562_tfs_7_cbc_gbc_dict_new.csv")
  download_and_compress_data_file(base_path, "pt2_concat_all.csv", "GSM2396859_k562_tfs_13_cbc_gbc_dict_new.csv")

  dir.create("results", showWarnings = FALSE)
  source(system.file("scripts", "4ab_5ab_9abc_quality_control_preprocessing.R", package = "pertInv"), local=TRUE)
}

ensure_processed_data <- function() {
  if (!file.exists(file.path("data_processed", "GSM2396858_k562_tfs_7", "guide_matrix.RData"))){
    warning("processed data not found - running setup_data ...\n")
    setup_data()
  }
}


#' Reproduce an individual figure
#' @param figure_number character
#' @export
reproduce_figure <- function(figure_number = c("1", "4", "5", "6", "7", "8", "9", "10", "11", "13", "14", "15", "16", "17", "18", "19", "20", "21", "23", "24", "A1")) {
  figure_number <- match.arg(figure_number)
  ensure_processed_data()
  cat("reproducing figure ", figure_number, "...\n")
  switch(figure_number,
         "1" = source(system.file("scripts", "01bc_simulate_example_graph.R", package = "pertInv"), local=TRUE),
         "A1" = ,
         "6" = source(system.file("scripts", "06ab_A1_mean_variance_relationship.R", package = "pertInv"), local=TRUE),
         "7" = source(system.file("scripts", "07_Bayesian_deconvolution.R", package = "pertInv"), local=TRUE),
         "8" = source(system.file("scripts", "08ab_10b_tSNE.R", package = "pertInv"), local=TRUE),
         "10" = {
           source(system.file("scripts", "08ab_10b_tSNE.R", package = "pertInv"), local=TRUE)
           source(system.file("scripts", "10a_13_all_guide_effects.R", package = "pertInv"), local=TRUE)
         },
         "13" = source(system.file("scripts", "10a_13_all_guide_effects.R", package = "pertInv"), local=TRUE),
         "11" = source(system.file("scripts", "11_genewise_variacne_decomposition.R", package = "pertInv"), local=TRUE),
         "14" = ,
         "15" = source(system.file("scripts", "14_15_simulate_effect_of_noise_on_icp.R", package = "pertInv"), local=TRUE),
         "16" = ,
         "17" = ,
         "18" = ,
         "19" = ,
         "20" = source(system.file("scripts", "16_17ab_18_19_20_knockout_inference.R", package = "pertInv"), local=TRUE),
         "4" = ,
         "5" = ,
         "9" = source(system.file("scripts", "4ab_5ab_9abc_quality_control_preprocessing.R", package = "pertInv"), local=TRUE),
         "21" = source(system.file("scripts", "21_additional_advantage_of_inferred_knockout_state.R", package = "pertInv"), local=TRUE),
         "23" = ,
         "24" = source(system.file("scripts", "23ab_24_stan_simulate_and_fit.R", package = "pertInv"), local=TRUE)
  )
  cat("done. Produced figure(s) can be found in folder 'results'.\n")
}

#' Reproduce all figures
#'
#' More efficient than producing each figure individually.
#' @export
reproduce_all <- function() {
  ensure_processed_data()
  cat("reproducing all figures...\n")
  cat("reproducing figure 1...\n")
  source(system.file("scripts", "01bc_simulate_example_graph.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figures 6 and A1...\n")
  source(system.file("scripts", "06ab_A1_mean_variance_relationship.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figure 7...\n")
  source(system.file("scripts", "07_Bayesian_deconvolution.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figures 8 and 10...\n")
  source(system.file("scripts", "08ab_10b_tSNE.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figures 10 and 13...\n")
  source(system.file("scripts", "10a_13_all_guide_effects.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figure 11...\n")
  source(system.file("scripts", "11_genewise_variacne_decomposition.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figures 14 and 15...\n")
  source(system.file("scripts", "14_15_simulate_effect_of_noise_on_icp.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figures 16-20...\n")
  source(system.file("scripts", "16_17ab_18_19_20_knockout_inference.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figure 21...\n")
  source(system.file("scripts", "21_additional_advantage_of_inferred_knockout_state.R", package = "pertInv"), local=TRUE)
  cat("done.\n")
  cat("reproducing figures 23 and 24...\n")
  source(system.file("scripts", "23ab_24_stan_simulate_and_fit.R", package = "pertInv"), local=TRUE)
  cat("done reproducing all figures. Produced figures can be found in folder 'results'.\n")
}

