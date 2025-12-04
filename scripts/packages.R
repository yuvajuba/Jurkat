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
  
  ## General packages
  "readxl",         
  "writexl",        
  "dplyr",          
  "tibble",         
  "tidyr",          
  "ggplot2",        
  "ggfortify",      # to use the autoplot() function directly for pca  
  "ggrepel",   
  "DT",
  "forcats",
  "ggvenn",
  
  ## Differential expression analysis
  "DESeq2",         
  "ComplexHeatmap", 
  "circlize",
  
  ## Pathway enrichment analysis
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "igraph",
  "ggraph",
  "tidygraph"
)

## 3. Load all required packages
load_or_install(List_packages)


#### Cleaning up !
rm(load_or_install, List_packages)