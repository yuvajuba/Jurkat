## 1. Function to load or Install Packages
load_or_install <- function(pkgs){
  for(pkg in pkgs){
    if (!require(pkg, character.only = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
}

## 2. List of packages used in this project
List_packages <- c(
  "readxl",    # read xlsx files
  "dplyr",     # data manipulation
  "tibble",    # data manipulation
  "tidyr",     # data manipulation
  "ggplot2",   # nice plotting
  "writexl",   # save files as xlsx
  "ggfortify", # to use the autoplot() function directly for pca  
  "DESeq2",    # differential expression analysis
  "ggrepel"
)

## 3. Load all required packages
load_or_install(List_packages)


#### Cleaning up !
rm(load_or_install, List_packages)