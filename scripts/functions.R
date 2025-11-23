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
    theme_bw()
  
  
  ## Regrouping the outputs :
  Outs <- list(raw_filt = raw_filt,
               DESeqData = dds,
               PCA_plt = plt)
  
  return(Outs)
}



## ---------------------------------------------- 3. Volcano plot (DEA)
My_volcano <- function(df, 
                       x= "log2FoldChange", 
                       y= "padj", 
                       lab= "Gene_name",
                       point_size= 2,
                       label_genes= F,
                       lab_size= 3,
                       lab_max_overlap= 20,
                       lab_box.padding= 0.5,
                       lab_colour = "navy", 
                       lab_fill = "white", 
                       lab_force_pull = .5,
                       selected_genes= NULL,
                       logFC_thresh= 1.5,
                       padj_thresh= 1e-5,
                       col_palette= NULL,
                       extend_xlim= 0.5){
  
  if(!(is.data.frame(df))){
    stop("your object isn't a dataframe !")
  }
  
  if(!all(c(x, y) %in% colnames(df))){
    stop(paste0("<",x,"> or <",y,"> isn't found in your dataframe !"))
  }
  
  hline <- -log10(padj_thresh)
  vline <- logFC_thresh
  
  ## Select genes :
  if(is.null(selected_genes)){
    selected_genes <- df %>% 
      dplyr::filter(.data[[y]] < 1e-10,
                    abs(.data[[x]]) > 2) %>% 
      dplyr::pull(.data[[lab]])
  }
  
  ## Color palette :
  if(is.null(col_palette)){
    col_palette <- c("pval & logFC"="darkred",
                     "logFC"="darkgreen",
                     "p-value"="midnightblue",
                     "NS"="lightgray")
  } else {
    if(length(col_palette) < 4){
      warning("The provided palette is of length < 4 ; therefor we're using default palette")
      col_palette <- c("pval & logFC"="darkred",
                       "logFC"="darkgreen",
                       "p-value"="midnightblue",
                       "NS"="lightgray")
      
    } else {
      col_palette <- c("pval & logFC" = col_palette[1],
                       "logFC" = col_palette[2],
                       "p-value" = col_palette[3],
                       "NS" = col_palette[4])
    }
  }
  
  
  ## plot :
  df %>% 
    dplyr::mutate(Col = case_when(abs(.data[[x]]) > logFC_thresh & 
                                    .data[[y]] < padj_thresh ~ "pval & logFC",
                                  abs(.data[[x]]) > logFC_thresh & 
                                    .data[[y]] > padj_thresh ~ "logFC",
                                  abs(.data[[x]]) < logFC_thresh & 
                                    .data[[y]] < padj_thresh ~ "p-value",
                                  abs(.data[[x]]) < logFC_thresh & 
                                    .data[[y]] > padj_thresh ~ "NS",
                                  TRUE ~ "NS"),
                  Sel_genes = ifelse(.data[[lab]] %in% selected_genes,
                                     .data[[lab]],
                                     NA)) -> df.plot
  
  
  
  if(isFALSE(label_genes)){
    
    plt <- df.plot %>% 
      ggplot(aes(x= .data[[x]],
                 y= -log10(.data[[y]])))+
      geom_vline(xintercept = -(vline), linetype = 2, color = "black", linewidth = .4)+
      geom_vline(xintercept = vline, linetype = 2, color = "black", linewidth = .4)+
      geom_hline(yintercept = hline, linetype = 2, color = "black", linewidth = .4)+
      geom_point(aes(colour = Col),
                 size = point_size,
                 alpha = 0.6)+
      labs(title = "Volcano Plot",
           caption = paste0("Number of genes : ", nrow(df)))+
      theme(legend.position = "top",
            legend.margin = margin(b=.05, unit = "in"),
            legend.key = element_rect(colour = "white"),
            panel.background = element_rect(fill = "white",
                                            colour = "darkgray",
                                            linewidth = .5),
            legend.text = element_text(size = 13,
                                       color = "black",
                                       face = "bold"),
            axis.title = element_text(size = 14,
                                      face = "bold",
                                      colour = "darkred"),
            axis.title.x = element_text(margin = margin(t=12)),
            axis.title.y = element_text(margin = margin(r=12)),
            axis.text = element_text(size = 12,
                                     colour = "black"),
            plot.title = element_text(size = 16,
                                      face = "bold",
                                      colour = "darkred",
                                      margin = margin(b=12)),
            plot.caption = element_text(size = 13,
                                        colour = "black"))+
      guides(colour = guide_legend(title = NULL,
                                   override.aes = list(size = 5)))+
      scale_color_manual(values = col_palette)+
      xlim(c((min(df[[x]]) - extend_xlim),(max(df[[x]]) + extend_xlim)))
    
  } else {
    
    plt <- df.plot %>% 
      ggplot(aes(x= .data[[x]],
                 y= -log10(.data[[y]])))+
      geom_vline(xintercept = -(vline), linetype = 2, color = "black", linewidth = .4)+
      geom_vline(xintercept = vline, linetype = 2, color = "black", linewidth = .4)+
      geom_hline(yintercept = hline, linetype = 2, color = "black", linewidth = .4)+
      geom_point(aes(colour = Col),
                 size = point_size,
                 alpha = 0.6)+
      labs(title = "Volcano Plot",
           caption = paste0("Number of genes : ", nrow(df)))+
      theme(legend.position = "top",
            legend.margin = margin(b=.05, unit = "in"),
            legend.key = element_rect(colour = "white"),
            panel.background = element_rect(fill = "white",
                                            colour = "darkgray",
                                            linewidth = .5),
            legend.text = element_text(size = 13,
                                       color = "black",
                                       face = "bold"),
            axis.title = element_text(size = 14,
                                      face = "bold",
                                      colour = "darkred"),
            axis.title.x = element_text(margin = margin(t=12)),
            axis.title.y = element_text(margin = margin(r=12)),
            axis.text = element_text(size = 12,
                                     colour = "black"),
            plot.title = element_text(size = 16,
                                      face = "bold",
                                      colour = "darkred",
                                      margin = margin(b=12)),
            plot.caption = element_text(size = 13,
                                        colour = "black"))+
      guides(colour = guide_legend(title = NULL,
                                   override.aes = list(size = 5)))+
      scale_color_manual(values = col_palette)+
      xlim(c((min(df[[x]]) - extend_xlim),(max(df[[x]]) + extend_xlim)))+
      geom_label_repel(aes(label = .data[["Sel_genes"]]),
                       size = lab_size,
                       max.overlaps = lab_max_overlap, 
                       box.padding = lab_box.padding,
                       label.r = 0.4,
                       label.size = .25, 
                       force_pull = lab_force_pull,
                       colour = lab_colour, 
                       fill = lab_fill)
    
  }
  
  
  return(plt)
}
