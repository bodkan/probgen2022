# 2022-03-24 17:47

# set up the slendr package for Binder
install.packages("rgdal")
devtools::install_github("bodkan/slendr")

# configure a dedicated Python environment for slendr
slendr::setup_env(agree = TRUE)

# packages used in the examples
install.packages(c("cowplot", "forcats", "BiocManager"))
BiocManager::install("ggtree")
