#' Download and extract data files from source repositories.
#' @export
setup_data <- function() {
  dir.create("data_raw", showWarnings = FALSE)
  base_path <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90063/suppl/"
  download_data_file <- function(base_path, filename) download.file(paste0(base_path, filename), file.path("data_raw", filename), method = "auto", mode="wb")

  if (!file.exists(file.path("data_raw", "GSM2396858_k562_tfs_7.mtx.txt.gz"))){
    download_data_file(base_path, "GSE90063_RAW.tar")
    untar(file.path("data_raw", "GSE90063_RAW.tar"), exdir="data_raw")
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

