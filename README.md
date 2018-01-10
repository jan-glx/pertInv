## pertInv - An R package to reproduce the analysis and figures presented in my master thesis


Usage:
```
# source("http://bioconductor.org/biocLite.R")
# install.packages("devtools")
devtools::install_github("jan-glx/pertInv")
pertInv::reproduce_all()
```

This will download and install this package from this repository, download the raw data files from their respective repositories, preprocess the data and generate all figures. It will take less than a few hours to complete. Instead of producing all figures at once, individual figures can be produced using, e.g., `pertInv::reproduce_figure("8")`. The coresponding analysis scripts can be found in the `scripts` folder.
