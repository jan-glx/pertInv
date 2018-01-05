## pertINV - An R package to reproduce the analysis and figures presented in my master thesis

### Usage
```
git clone git@github.com:jan-glx/pertInv.git
cd pertInv/data_raw 
bash get_data.sh
cd ..
R
> devtools::install(".", dependencies = TRUE)
> source("scripts/000_4ab_5ab_9abc_quality_control_preprocessing.R")
```

This will download this repository, download the data files, install this package to provide key functions for the scripts and run the quality control/preprocessing script, producing figures 4a, 4b, 5a, 5b, and 9 a,b and c. After that other scripts in the script folder can be run to reproduce further figures shown in the thesis.