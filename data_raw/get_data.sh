wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90063/suppl/GSE90063_RAW.tar
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90063/suppl/GSE90063_dc0hr_umi_wt.txt.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90063/suppl/GSE90063_dc3hr_umi_wt.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90063/suppl/GSE90063_k562_umi_wt.txt.gz
wget https://raw.githubusercontent.com/asncd/MIMOSCA/master/GBC_CBC_pairing/gbc_cbc_dicts/dc_0hr_concat_all.csv
wget https://raw.githubusercontent.com/asncd/MIMOSCA/master/GBC_CBC_pairing/gbc_cbc_dicts/ph_concat_all.csv
wget https://raw.githubusercontent.com/asncd/MIMOSCA/master/GBC_CBC_pairing/gbc_cbc_dicts/promoters_concat_all.csv
wget https://raw.githubusercontent.com/asncd/MIMOSCA/master/GBC_CBC_pairing/gbc_cbc_dicts/pt2_concat_all.csv
tar -xvf GSE90063_RAW.tar
rm GSE90063_RAW.tar
gzip < promoters_concat_all.csv > GSM2396858_k562_tfs_7_cbc_gbc_dict_new.csv.gz
gzip < pt2_concat_all.csv > GSM2396859_k562_tfs_13_cbc_gbc_dict_new.csv.gz
gzip < dc_0hr_concat_all.csv > GSM2396857_dc_0hr_cbc_gbc_dict_new.csv.gz
gzip < ph_concat_all.csv > GSM2396860_k562_tfs_highmoi_lenient_cbc_gbc_dict_new.csv.gz
