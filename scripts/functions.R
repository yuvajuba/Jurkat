## ---------------------------------------------- 1. ggplot scatter plot theme
my_scatter_theme <- function(){
  theme_bw()+
    theme(
      ## Backgrounds:
      plot.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.key = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      plot.title = element_text(colour = "darkcyan", face = "bold", margin = margin(b=0.15, unit = "in"), size = 14),
      legend.title = element_text(size = 12, colour = "darkcyan", face = "bold")
    )
}



## ---------------------------------------------- 2. RNAseq rawcounts filtering and PCA
Pre_process <- function(RCounts,
                        Conditions,
                        Metadata,
                        MaxCount_threshold = 20,
                        CpmCount_threshold = 0.5,
                        MinSample = 4,
                        Biotypes_select = c("protein_coding", "lncRNA")){
  
  ###   RCount    : The raw count data
  ###   Conditions: The experimental data
  
  rawcounts <- RCounts
  Experimental <- Conditions
  
  # filtering by max count threshold :
  raw_filt <- rawcounts[rowMaxs(as.matrix(rawcounts)) >= MaxCount_threshold, ]
  
  # filtering by cpm threshold :
  cpm_counts <- edgeR::cpm(raw_filt)
  raw_filt <- raw_filt[rowSums(cpm_counts > CpmCount_threshold) >= MinSample,]
  
  # filtering by biotypes :
  to_keep <- rownames(Metadata)[Metadata$Gene_biotype %in% Biotypes_select]
  raw_filt <- raw_filt[rownames(raw_filt) %in% to_keep, ]
  
  
  ## 1: build DESeq2 structure
  dds <- DESeqDataSetFromMatrix(countData = raw_filt,
                                colData = Experimental,
                                design = ~ Condition)
  
  ## 2: apply vst
  vsd <- vst(dds, blind = T)
  #### blind = TRUE (gives a neutral VST - not influenced by the design 'condition')
  
  ## 3: extract vst matrix 
  vst_mat <- assay(vsd)
  
  ## 4: run PCA using prcomp()
  pca <- prcomp(t(vst_mat), scale. = TRUE)
  
  var_explained <- round(100*(pca$sdev)^2 / sum(pca$sdev^2), 2) %>% setNames(colnames(pca$x))
  
  ## 5: data to plot
  plt_data <- pca$x %>% 
    as.data.frame() %>% 
    dplyr::mutate(Condition = Experimental$Condition,
                  Replicat = Experimental$Replicat,
                  Sample = rownames(Experimental))
  
  ## 6: plot
  plt <- plt_data %>% 
    ggplot(aes(x= PC1, y= PC2))+
    geom_point(aes(colour = Condition), size = 3)+
    geom_text_repel(aes(label = Sample, colour = Condition), size = 3, show.legend = FALSE)+
    labs(title = "PCA of Jurkat_ct vs Jurkat_H151",
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"))+
    scale_colour_manual(values = c("Jurkat_H151"="brown", "Jurkat_ct"="black"))+
    theme_bw()
  
  
  ## Regrouping the outputs :
  Outs <- list(raw_filt = raw_filt,
               DESeqData = dds,
               PCA_plt = plt)
  
  return(Outs)
}
