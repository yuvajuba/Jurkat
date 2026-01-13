## ---------------------------------------------- 1. ggplot general theme
my_general_theme <- function(){
  theme_bw()+
    theme(
      ## Backgrounds:
      plot.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.key = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      panel.background = element_rect(colour = "gray", fill = "white"),
      ## Titles:
      plot.title = element_text(colour = "darkcyan", face = "bold", margin = margin(b=0.1, unit = "in"), size = 14),
      plot.subtitle = element_text(colour = "darkgreen", face = "bold", margin = margin(b=0.1, unit = "in"), size = 13),
      axis.title = element_text(colour = "darkcyan", face = "bold", size = 13),
      axis.title.x = element_text(margin = margin(t=0.1, unit = "in")),
      axis.title.y = element_text(margin = margin(r=0.1, unit = "in")),
      legend.title = element_text(size = 13, colour = "darkcyan", face = "bold", margin = margin(b=0.1, unit = "in")),
      ## Textes:
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12, face = "bold"),
      plot.caption = element_text(size = 11, face = "bold", margin = margin(t=0.15, unit = "in"))
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
                       label_genes_xthreshold= 2,
                       label_genes_ythreshold= 1e-20,
                       lab_size= 3,
                       lab_max_overlap= 20,
                       lab_box.padding= 0.5,
                       lab_colour = "navy",
                       lab_force = 1,
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
      dplyr::filter(.data[[y]] < label_genes_ythreshold,
                    abs(.data[[x]]) > label_genes_xthreshold) %>% 
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
  
  
  ## plot dataframe:
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
  
  
  ## plot:
  plt <- df.plot %>% 
    ggplot(aes(x= .data[[x]],
               y= -log10(.data[[y]])))+
    geom_vline(xintercept = -(vline), linetype = 2, color = "black", linewidth = .4)+
    geom_vline(xintercept = vline, linetype = 2, color = "black", linewidth = .4)+
    geom_hline(yintercept = hline, linetype = 2, color = "black", linewidth = .4)+
    geom_point(aes(colour = Col),
               size = point_size,
               alpha = 0.6)+
    labs(title = "Differential expressed genes",
         caption = paste0("Number of genes : ", nrow(df)))+
    
    theme(
      plot.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.key = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      panel.background = element_rect(colour = "gray", fill = "white"),
      plot.title = element_text(colour = "darkcyan", face = "bold", margin = margin(b=0.15, unit = "in"), size = 14),
      axis.title = element_text(colour = "darkcyan", face = "bold", size = 13),
      axis.title.x = element_text(margin = margin(t=0.15, unit = "in")),
      axis.title.y = element_text(margin = margin(r=0.15, unit = "in")),
      axis.text = element_text(size = 11),
      legend.text = element_text(size = 12, face = "bold", colour = "black"),
      legend.position = "top",
      legend.margin = margin(b=.05, unit = "in"),
      plot.caption = element_text(size = 12, colour = "black")
    )+
    guides(colour = guide_legend(title = NULL,
                                 override.aes = list(size = 5)))+
    scale_color_manual(values = col_palette)+
    xlim(c((min(df[[x]]) - extend_xlim),(max(df[[x]]) + extend_xlim)))
  
  
  if(isTRUE(label_genes)){
    plt <- plt+
      geom_text_repel(aes(label = .data[["Sel_genes"]]),
                      size = lab_size,
                      max.overlaps = lab_max_overlap,
                      box.padding = lab_box.padding, 
                      force_pull = lab_force_pull, 
                      colour = lab_colour,
                      show.legend = F,
                      force = lab_force)
  }
  
  
  return(plt)
}



## ---------------------------------------------- 4. Venn diagram
Display_Venn <- function(Markers, 
                         colpalette = NULL, 
                         set.names = NULL, 
                         set.name.size = 5,
                         text.size = 5,
                         Padding = 0.04,
                         show_elements = FALSE) {
  
  
  
  # 1. Input checks :
  # --------------- :
  if (!is.list(Markers)) stop("Argument 'Markers' must be a list of vectors.")
  if (!all(sapply(Markers, is.vector))) stop("All elements in the 'Markers' list must be vectors.")
  
  
  n <- min(length(Markers), 4)
  if (n < 2) stop("You need at least 2 sets to make a Venn diagram.")
  if (length(Markers) > 4) warning("Only the first 4 sets will be used.")
  
  
  # 2. Handling sets :
  # ---------------- :
  sets <- Markers[1:n]
  if (!(is.null(set.names))) {
    if (length(set.names) != n) {
      stop("The length of 'set.names' must match the number of sets in 'Markers'")
    }
    names(sets) <- set.names
  } else if (is.null(names(sets))) {
    names(sets) <- paste0("Set", 1:n)
  }
  
  
  # 3. Handling colours :
  # ------------------- :
  default_palette <- c("navy", "red", "darkgreen", "violet")
  colpalette <- if (is.null(colpalette) || length(colpalette) < n) {
    warning("Invalid color palette. Using default.")
    default_palette[1:n]
  } else {
    colpalette[1:n]
  }
  
  
  
  # 4. Calculate intersections :
  # -------------------------- :
  common_genes <- list()
  for (i in 2:n) {
    combs <- combn(sets, i, simplify = FALSE)
    for (x in combs) {
      key <- paste(names(x), collapse = " ∩ ")
      common_genes[[key]] <- purrr::reduce(x, intersect)
    }
  }
  
  
  # 5. Calculate group specific genes :
  # --------------------------------- :
  specific_genes <- list()
  for (i in 2:n) {
    combs <- combn(sets, i, simplify = FALSE)
    for (x in combs) {
      key <- paste(names(x), collapse = " | ")
      specific <- purrr::map2(x, names(x), ~ setdiff(.x, purrr::reduce(x[names(x) != .y], union)))
      specific_genes[[key]] <- specific
    }
  }
  
  # Plot
  p <- ggvenn::ggvenn(
    sets,
    fill_color = colpalette,
    stroke_size = 0.8,
    set_name_color = colpalette,
    set_name_size = set.name.size,
    text_size = text.size,
    text_color = "black",
    stroke_linetype = "solid",
    padding = Padding,
    fill_alpha = 0.4,
    show_stats = "c",
    show_elements = show_elements
  ) +
    theme(plot.background = element_rect(fill = "white", color = "white"))
  
  return(list(plot = p, intersections = common_genes, group_specific = specific_genes))
}



## ---------------------------------------------- 5. Network plot
Net_plot <- function(GOobject,
                     DEAres,
                     condition.name = "condition",
                     select.terms = NULL,
                     nterm = 6,
                     net.layout = "nicely",
                     GOobject.colname.term = "Description",
                     GOobject.colname.gene = "geneID",
                     GOobject.colname.pval = "p.adjust",
                     DEAres.colname.gene = "Gene_name",
                     DEAres.colname.logFC = "log2FoldChange",
                     DEAres.colname.pval = "padj",
                     edge.param.color = NULL,
                     edge.param.width = .8,
                     edge.param.alpha = .5,
                     edge.legend.width = 1.6,
                     edge.legend.position = "bottom",
                     edge.legend.direction = "horizontal",
                     edge.legend.layers = 2,
                     edge.legend.fontsize = 11,
                     node.param.term.color = "black",
                     node.param.term.size = c(8,4),
                     node.param.gene.logFC = T,
                     node.param.gene.color = "darkcyan",
                     node.param.gene.size = 3,
                     node.param.gene.alpha = .7,
                     label.genes = T,
                     label.gene.color = c("red","purple","black"),
                     label.gene.size = 2,
                     label.gene.alpha = 0.9,
                     label.gene.repel = T){
  
  ################ :
  ## 01. Debugging :
  ################ :
  if(!inherits(GOobject, "enrichResult")) {
    stop("You need to provide a GO object from enrichGO or similar.")
  }
  
  if(!(is.data.frame(DEAres))){
    stop("DEAres needs to be a dataframe")
  }
  
  if(!(net.layout %in% c("fr", "kk", "stress", "tree", "mds", "circle", "nicely"))) {
    stop(paste("Please choose a layout within this list : ",
               "fr, kk, stress, tree, mds, circle, gem, drl, nicely",
               sep = "\n"))
  }
  
  
  
  ################ :
  ## 02. Standards :
  ################ :
  GOresults <- GOobject@result %>% 
    dplyr::filter(.data[[GOobject.colname.pval]] < 0.05) %>% 
    dplyr::arrange(.data[[GOobject.colname.pval]]) %>% 
    dplyr::select(Term = GOobject.colname.term,
                  Gene = GOobject.colname.gene)
  
  DEAresults <- DEAres %>% 
    dplyr::filter(.data[[DEAres.colname.pval]] < 0.05) %>% 
    dplyr::select(Gene = DEAres.colname.gene, 
                  logFC = DEAres.colname.logFC)
  
  if(!(is.null(select.terms))){
    selected_term <- select.terms
  } else {
    selected_term <- GOresults$Term[1:nterm]
  }
  
  
  
  ############################ :
  ## 03. Building network data :
  ############################ :
  dat <- GOresults %>% 
    dplyr::mutate(Cluster = condition.name)  %>% 
    tidyr::separate_rows(Gene, sep = "/") %>% 
    left_join(DEAresults, by = "Gene")
  
  dat <- dat[dat$Term %in% selected_term,] %>% 
    dplyr::mutate(Regulation = case_when(logFC > 0 ~ "Up",
                                         logFC < 0 ~ "Down"))
  
  ## Adding group colour to colour shared and distinct genes
  gene_term_map <- dat %>% 
    dplyr::select(Gene, Term) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(n_term = n(),
                     term_for_color = ifelse(n_term == 1, Term, NA_character_),
                     .groups = "drop") %>% 
    dplyr::mutate(color_group = case_when(
      n_term == 1 ~ term_for_color,
      n_term == 2 ~ "shared by two",
      n_term >= 3 ~ "shared by more"
    )) %>% 
    dplyr::select(Gene, color_group)
  
  ## Creating genes info data
  genes_info <- dat %>%
    dplyr::select(Gene, Cluster, logFC, Regulation) %>%
    dplyr::distinct(Gene, .keep_all = TRUE) %>% 
    dplyr::left_join(gene_term_map, by = "Gene")
  
  
  
  ########################### :
  ## 04. Building igraph data :
  ########################### :
  graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Term")], directed = FALSE)
  
  ## Annotating nodes
  V(graph_data)$type <- ifelse(V(graph_data)$name %in% genes_info$Gene, "Gene", "Term")
  V(graph_data)$genes <- case_when(V(graph_data)$type == "Gene" ~ V(graph_data)$name, TRUE ~ NA)
  V(graph_data)$terms <- case_when(V(graph_data)$type == "Term" ~ V(graph_data)$name, TRUE ~ NA)
  V(graph_data)$logFC <- genes_info$logFC[match(V(graph_data)$name, genes_info$Gene)]
  V(graph_data)$reg <- genes_info$Regulation[match(V(graph_data)$name, genes_info$Gene)]
  V(graph_data)$lab_color <- genes_info$color_group[match(V(graph_data)$name, genes_info$Gene)]
  
  ## Annotating edges
  E(graph_data)$term <- dat$Term[match(
    paste0(dat$Gene, dat$Term),
    paste0(ends(graph_data, es = E(graph_data))[, 1],
           ends(graph_data, es = E(graph_data))[, 2])
  )]
  
  ## Final igraph data
  graph_tbl <- tidygraph::as_tbl_graph(graph_data)
  
  
  
  ########################## :
  ## 05. Plotting parameters :
  ########################## :
  MyPalette <- c("#070","#700","#007","#770","#077","#707","#2a0",
                 "#a30","#40a","#ca0","#05c","#b04","#000","#f0f")
  
  nterm <- length(selected_term)
  
  ## Set edge colours
  if(is.null(edge.param.color)){
    if(nterm <= length(MyPalette)){
      e.colors <- MyPalette[1:nterm] %>% setNames(selected_term)
    }
    
  } else if(nterm <= length(edge.param.color)) {
    e.colors <- edge.param.color[1:nterm] %>% setNames(selected_term)
  } else {
    e.colors <- sample(colors(), nterm) %>% setNames(selected_term)
    message("Provided colours don't match the number of terms \nUsing shuffeled colours!")
  }
  
  ## Set edge legend guides
  edge_cols <- switch(
    edge.legend.direction,
    "horizontal" = guide_legend(title = NULL,
                                override.aes = list(edge_width = edge.legend.width),
                                position = edge.legend.position,
                                direction = "horizontal",
                                nrow = edge.legend.layers),
    "vertical" = guide_legend(title = NULL,
                              override.aes = list(edge_width = edge.legend.width),
                              position = edge.legend.position,
                              direction = "vertical",
                              ncol = edge.legend.layers)
  )
  
  ## Set logFC values
  logFC_values <- c(floor(min(genes_info$logFC)), ceiling(max(genes_info$logFC)))
  
  ## Set logFC colour gradient
  if((sign(logFC_values[1]) * sign(logFC_values[2])) == -1){
    gradient_logFC <- scale_colour_gradientn(
      colours = c("#050","#fff","#fff","#500"),
      values = scales::rescale(c(min(genes_info$logFC), -0.3, 0.3, max(genes_info$logFC))),
      limits = c(logFC_values[1], logFC_values[2])
    )
  } else {
    gradient_logFC <- scale_colour_gradient(
      low = ifelse(sign(logFC_values[1]) == -1, "#050", "#ffa"),
      high = ifelse(sign(logFC_values[1]) == -1, "#ffa", "#500")
    )
  }
  
  ## Set term nodes layers size
  node_term <- list(col1 = node.param.term.color[1],
                    col2 = ifelse(length(node.param.term.color) >= 2,
                                  node.param.term.color[2],
                                  node.param.term.color[1]),
                    size1 = node.param.term.size[1],
                    size2 = ifelse(length(node.param.term.size) >= 2,
                                   node.param.term.size[2],
                                   node.param.term.size[1]))
  
  ## Set label genes colours
  gene_lab_cols <- list(`shared by more` = label.gene.color[1],
                        `shared by two` = ifelse(length(label.gene.color) >= 2,
                                                 label.gene.color[2],
                                                 label.gene.color[1]),
                        distinct = ifelse(length(label.gene.color) >= 3,
                                          label.gene.color[3],
                                          "black"))
  
  
  
  ########### :
  ## 06. Plot :
  ########### :
  g <- ggraph(graph_tbl, layout = net.layout)+
    theme_void()
  
  ## Plot edges
  g <- g +
    geom_edge_link(aes(color = term), alpha = edge.param.alpha, edge_width = edge.param.width)+
    scale_edge_color_manual(values = e.colors)+
    guides(edge_color = edge_cols)
  
  ## Plot gene nodes
  if(isFALSE(node.param.gene.logFC)){
    g <- g +
      geom_node_point(colour = node.param.gene.color,
                      alpha = node.param.gene.alpha,
                      size = ifelse(V(graph_tbl)$type == "Gene", node.param.gene.size, 0))
  } else {
    g <- g +
      geom_node_point(aes(colour = logFC),
                      alpha = node.param.gene.alpha,
                      size = ifelse(V(graph_tbl)$type == "Gene", node.param.gene.size, 0))+
      gradient_logFC+
      guides(colour = guide_colourbar(title = "log2FC",
                                      barwidth = unit(.6, "cm"),
                                      barheight = unit(3, "cm")))
  }
  
  ## Plot term nodes
  g <- g +
    #### first layer
    geom_node_point(colour = ifelse(V(graph_tbl)$type == "Term", node_term$col1, "white"),
                    size = ifelse(V(graph_tbl)$type == "Term", node_term$size1, 0),
                    alpha = ifelse(V(graph_tbl)$type == "Term", 0.3, 0),)+
    #### second layer
    geom_node_point(colour = ifelse(V(graph_tbl)$type == "Term", node_term$col2, "white"),
                    size = ifelse(V(graph_tbl)$type == "Term", node_term$size2, 0),
                    alpha = ifelse(V(graph_tbl)$type == "Term", 1, 0),)
  
  ## Plot gene labels
  if(isTRUE(label.genes)){
    g <- g +
      geom_node_text(aes(label = genes),
                     show.legend = F,
                     fontface = "bold.italic",
                     size = label.gene.size,
                     alpha = label.gene.alpha,
                     repel = label.gene.repel,
                     colour = ifelse(V(graph_tbl)$lab_color == "shared by more",
                                     gene_lab_cols$`shared by more`,
                                     ifelse(V(graph_tbl)$lab_color == "shared by two",
                                            gene_lab_cols$`shared by two`,
                                            gene_lab_cols$distinct)))
  }
  
  
  
  ## Setup plot theme
  g <- g +
    theme(
      #### Backgrounds
      plot.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      #### Plot
      plot.margin = margin(b=0.1,t=0.1,r=0.1,l=0.1, unit = "in"),
      plot.caption = element_text(size = 11,
                                  face = "bold.italic",
                                  colour = "#333",
                                  margin = margin(t=0.25, unit = "in")),
      #### Legends
      legend.title = element_text(size = 14,
                                  face = "bold",
                                  colour = "darkcyan",
                                  margin = margin(b=0.15, unit = "in")),
      legend.text = element_text(size = edge.legend.fontsize,
                                 face = "bold"),
      legend.margin = margin(l=0.2,b=0.2,r=0.1,t=0.2, unit = "in")
    )
  
  
  
  ## Setup caption for labels colour info
  shared <- length(genes_info$color_group[genes_info$color_group == "shared by more"])
  shared2 <- length(genes_info$color_group[genes_info$color_group == "shared by two"])
  notshared <- nrow(genes_info)-(shared+shared2)
  
  g <- g +
    labs(caption = paste(
      sprintf("%-30s: %s", 
              paste0("Distinct genes (", gene_lab_cols$distinct,")"),
              notshared),
      sprintf("%-30s: %s", 
              paste0("Genes shared by two (", gene_lab_cols$`shared by two`,")"),
              shared2),
      sprintf("%-30s: %s", 
              paste0("Genes shared by many (", gene_lab_cols$`shared by more`,")"),
              shared),
      sep = "\n"
    ))
  
  return(list(netdata = dat, geneInfoData = genes_info, plot = g))
  
  
}



## ---------------------------------------------- 6. Normalizing & scaling data
nlog_scale_data <- function(RCounts, 
                            Sel_genes){
  
  # -- Checks ----
  if (!is.data.frame(RCounts) & !is.matrix(RCounts)) {
    stop("RCounts must be a data.frame or matrix.")
  }
  if(all(apply(RCounts,2,class) != "numeric")){
    stop("All your variables must be numeric values !")
  }
  
  if(is.null(names(Sel_genes)) || any(names(Sel_genes) == "")){
    stop("Sel_genes must be a named vector with the rownames of your RCounts")
  }
  if(all(!(names(Sel_genes) %in% rownames(RCounts)))){
    stop(paste("There are some unmatched genes between RCounts and Sel_genes",
               "rownames(RCounts) need to match names(Sel_genes)",
               sep = "\n"))
  }
  
  # -- Proceed ----
  ## 1. with scaling :
  n_fact <- estimateSizeFactorsForMatrix(RCounts)
  n_count <- sweep(RCounts, 2, n_fact, FUN = "/")
  n_count <- log10(n_count+1)
  n_count <- n_count[names(Sel_genes),]
  n_count <- t(apply(n_count,1,scale))
  colnames(n_count) <- colnames(RCounts)
  
  ## 2. without scaling
  log_counts <- sweep(RCounts, 2, n_fact, FUN = "/")
  log_counts <- log10(log_counts+1)
  log_counts <- log_counts[names(Sel_genes),]
  
  # -- return ----
  return(list(scaled = n_count,
              unscaled = log_counts))
}



## ---------------------------------------------- 7. Enrichment mapping plot
mapTermsPlot <- function(goObject,
                         goSimData = NULL,
                         pairwiseTermsimMethod = "Wang",
                         selectTerms = 30,
                         similarityCutoff = 0.6,
                         networkLayout = "nicely",
                         nodesColour = NULL,
                         nodesAlpha = 0.9,
                         edgesAlpha = 0.4,
                         edgesColour = "#777",
                         plotTitleSize = 14,
                         legendTitleSize = 13,
                         legendTextSize = 11,
                         labelsColour = "#333",
                         labelsAlpha = 0.9,
                         labelsSize = 3,
                         labelParams = list(distinctTerms = F,
                                            labelBy = c("strength","Count"),
                                            labelsPerCluster = 3,
                                            withTies = T)){
  
  ## ---- 01. Debuging :
  if(!inherits(goObject, "enrichResult")) {
    stop("You need to provide a GO object from enrichGO or similar.")
  }
  
  if(nrow(goObject@termsim) == 0){
    if(is.null(goSimData)){
      stop("Please provide a pre-computed goObject or a goSimdata from GOSemSim::godata()")
    } else {
      goObject <- pairwise_termsim(goObject, 
                                   method = pairwiseTermsimMethod,
                                   semData = goSimData)
    }
  }
  
  if(!(networkLayout %in% c("fr", "kk", "stress", "tree", "mds", "circle", "nicely"))) {
    stop(paste("Please choose a layout within this list : ",
               "fr, kk, stress, tree, mds, circle, gem, drl, nicely",
               sep = "\n"))
  }
  
  if(is.numeric(selectTerms)){
    selected_terms <- goObject@result %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = selectTerms) %>%
      dplyr::pull(Description)
  } 
  
  if(is.character(selectTerms)){
    selected_terms <- selectTerms
    sig_terms <- goObject@result$Description[goObject@result$p.adjust < 0.05]
    
    if(!(all(selected_terms %in% goObject@result$Description[goObject@result$p.adjust < 0.05]))){
      message("Some selected terms are not significant, or don't exist!")
      outs <- dplyr::setdiff(selected_terms,
                             sig_terms)
      message("The intruders are :")
      cat(outs, sep = "\n")
      message("Keeping only significant terms")
      selected_terms <- dplyr::setdiff(selected_terms, outs)
    }
  }
  
  ## ---- 02. Handle Description vs ID in termsim :
  IDs <- FALSE
  selected_terms_ID <- NULL
  
  if(!all(selected_terms %in% rownames(goObject@termsim))){
    selected_terms_ID <- goObject@result$ID[goObject@result$Description %in% selected_terms]
    IDs <- TRUE
    
    if(!all(selected_terms_ID %in% rownames(goObject@termsim))){
      stop("Selected terms not found in the similarity matrix.")
    }
  }
  
  sel <- if (IDs) selected_terms_ID else selected_terms
  goObject@termsim <- goObject@termsim[sel, sel]
  
  
  ## ---- 02. Data preparation :
  #### GO results
  goResult <- goObject@result %>% 
    dplyr::filter(Description %in% selected_terms)
  
  #### Nodes data
  nodeData <- goResult %>% 
    dplyr::select(ID, Description, p.adjust, GeneRatio, Count, geneID)
  
  #### Edges data
  edgeData <- as.data.frame(as.table(goObject@termsim)) %>% 
    dplyr::rename(from = Var1,
                  to   = Var2,
                  weight = Freq)
  
  if(!IDs){
    ## termsim indexed by Description → convert to ID
    edgeData <- edgeData %>% 
      dplyr::left_join(nodeData[,c("ID","Description")], by = c("from" = "Description")) %>%
      dplyr::rename(from_id = ID) %>%
      dplyr::left_join(nodeData[,c("ID","Description")], by = c("to" = "Description")) %>%
      dplyr::filter(from != to, weight >= similarityCutoff) %>%
      dplyr::mutate(pair = paste(pmin(from, to), pmax(from, to))) %>%
      dplyr::distinct(pair, .keep_all = TRUE) %>%
      dplyr::select(from = from_id,
                    to = ID,
                    weight)
  } 
  
  edgeData <- edgeData %>%
    dplyr::filter(from != to, weight >= similarityCutoff) %>%
    dplyr::mutate(pair = paste(pmin(from, to), pmax(from, to))) %>%
    dplyr::distinct(pair, .keep_all = TRUE) %>%
    dplyr::select(-pair)
  
  
  ## ---- 03. Data enrichment :
  #### Terms membership
  tmp_graph <- graph_from_data_frame(d = edgeData,
                                     vertices = nodeData,
                                     directed = FALSE)
  cc <- components(tmp_graph)
  
  #### Enrich nodeData
  nodeData$cluster <- cc$membership[match(nodeData$ID, names(cc$membership))]
  nodeData$cluster <- as.character(nodeData$cluster)
  nodeData$cluster[degree(tmp_graph) == 0] <- "Distinct"
  nodeData$cluster <- factor(nodeData$cluster)
  nodeData$strength <- strength(
    tmp_graph,
    weights = edgeData$n_shared_genes
  )[match(nodeData$ID, names(strength(tmp_graph)))]
  
  nodeData$label <- nodeData$Description
  
  if(!(labelParams$labelBy %in% c("strength", "Count"))){
    stop("labelBy argument should be one of : 'strength' or 'Count' !")
  }
  
  keep_labs <- nodeData %>%
    dplyr::filter(cluster != "Distinct") %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::slice_max(order_by = .data[[labelParams$labelBy]], 
                     n = labelParams$labelsPerCluster,
                     with_ties = labelParams$withTies) %>% 
    dplyr::ungroup() %>% 
    dplyr::pull(label)
  
  nodeData <- nodeData %>% 
    dplyr::mutate(label = dplyr::if_else(label %in% keep_labs,
                                         label,
                                         NA_character_))
  
  if(labelParams$distinctTerms){
    nodeData <- nodeData %>% 
      dplyr::mutate(label = dplyr::if_else(cluster == "Distinct",
                                           Description,
                                           label))
  }
  
  #### Construct term-gene table
  term_genes <- goResult %>%
    dplyr::select(ID, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(genes = list(geneID), .groups = "drop")
  
  #### Enrich edgeData
  edgeData <- edgeData %>%
    dplyr::left_join(term_genes, by = c("from" = "ID")) %>%
    dplyr::rename(genes_from = genes) %>%
    dplyr::left_join(term_genes, by = c("to" = "ID")) %>%
    dplyr::rename(genes_to = genes) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_shared_genes = length(intersect(genes_from, genes_to)),
      shared_genes = list(intersect(genes_from, genes_to))
    ) %>%
    dplyr::ungroup()
  
  
  ## ---- 04. igraph :
  g <- graph_from_data_frame(d = edgeData %>% dplyr::select(from,to,weight,n_shared_genes),
                             vertices = nodeData,
                             directed = FALSE) %>% 
    as_tbl_graph()
  
  
  ## ---- 05. Plot parameters :
  #### Assigning colours
  col_palette <- c("#080","#a80","#a08","#008","#fa0","#a0f","#0af","#f0a",
                   "#af0","#3ff","#700","#fd0","#0fa","#02c","#b02","#000",
                   "#240","#057","#fec","#344","#c70","#c07","#666","#aaa") 
  
  if(!is.null(nodesColour)){
    if(length(nodesColour) >= length(levels(nodeData$cluster))){
      col_palette <- nodesColour
    }
  }
  
  clusters <- levels(nodeData$cluster)
  cluster_cols <- setNames(col_palette[seq_along(clusters)], clusters)
  
  
  ## ---- 06. Plot Basic :
  p <- ggraph(g, layout = networkLayout)+
    geom_edge_link(aes(width = n_shared_genes),
                   alpha = edgesAlpha,
                   colour = edgesColour)+
    geom_node_point(aes(size = Count, colour = cluster),
                    alpha = nodesAlpha)+
    geom_node_text(aes(label = label),
                   colour = labelsColour,
                   alpha = labelsAlpha,
                   size = labelsSize,
                   repel = T)+
    theme_void()
  
  
  ## ---- 07. Final plot :
  p <- p+
    labs(title = paste0("Enriched terms mapping: Cutoff similarity at ",
                        similarityCutoff))+
    scale_colour_manual(values = cluster_cols, 
                        name = "GO term clusters")+
    scale_edge_width(range = c(0.3, 2.5))+
    guides(colour = guide_legend(override.aes = list(size = 4)))+
    theme(
      plot.margin = margin(t=0.05,l=0.05,r=0.05,b=0.05, unit = "in"),
      legend.margin = margin(r=0.1,l=0.05, unit = "in"),
      
      plot.background = element_rect(colour = "white", fill = "white"),
      legend.background = element_rect(colour = "white", fill = "white"),
      
      plot.title = element_text(colour = "darkcyan", 
                                face = "bold", 
                                size = plotTitleSize),
      legend.title = element_text(size = legendTitleSize, 
                                  colour = "darkcyan", 
                                  face = "bold"),
      
      legend.text = element_text(size = legendTextSize, face = "bold")
    )
  
  
  
  return(list(plot = p,
              nodeData = nodeData))
}



## --------------------------------------------------- 8. Heatmap Gene-Terms
heatmapTermGene <- function(nodeData,
                            colnameTerms = 'Description',
                            colnameClusters = 'cluster',
                            colnameGenes = 'geneID',
                            deaData,
                            colnameGenesDea = 'Gene_name',
                            colnameLogFC = 'log2FoldChange',
                            geneList,
                            fillByLogFC = F,
                            plt_gp = list(tilesLineColour = "white",
                                          fillEmptyValues = "whitesmoke",
                                          fillPresence = "darkgreen",
                                          axisTextSize = 12,
                                          axisTextColour = "black",
                                          xAxisAngle = 30,
                                          xAxisMaxCharSize = 60,
                                          legendTitleSize = 13)){
  
  
  ## ---- 1/ debuging:
  if (!all(c(colnameTerms, colnameClusters, colnameGenes) %in% colnames(nodeData))) {
    stop("nodeData must contain the 3 columns! At least 1 is missing!")
  }
  
  if (!all(c(colnameGenesDea, colnameLogFC) %in% colnames(deaData))) {
    stop("deaData must contain the 2 columns! At least 1 is missing!")
  }
  
  
  
  ## ---- 2/ preparing the heatmap data:
  hmData <- nodeData %>% 
    dplyr::select(.data[[colnameTerms]],
                  .data[[colnameClusters]],
                  .data[[colnameGenes]]) %>% 
    tidyr::separate_rows(.data[[colnameGenes]], sep = "/") %>% 
    dplyr::mutate(present = 1) %>% 
    tidyr::pivot_wider(names_from = .data[[colnameGenes]], 
                       values_from = present, 
                       values_fill = 0) %>% 
    tidyr::pivot_longer(-c(.data[[colnameTerms]], .data[[colnameClusters]]), 
                        names_to = "gene", 
                        values_to = "present") %>% 
    dplyr::left_join(deaData[,c(colnameGenesDea, colnameLogFC)], 
                     by = c("gene" = colnameGenesDea)) %>% 
    dplyr::rename(term = .data[[colnameTerms]],
                  log2FC = .data[[colnameLogFC]]) %>% 
    dplyr::mutate(cluster = if_else(cluster == "Distinct", term, cluster)) %>% 
    dplyr::filter(gene %in% geneList) %>% 
    dplyr::group_by(gene, cluster) %>% 
    dplyr::summarise(present = max(present),
                     log2FC = dplyr::first(log2FC), .groups = "drop") %>% 
    dplyr::mutate(log2FC = if_else(present != 0, log2FC, NA))
  
  #### pull only clusters with genes!
  clustersWithGenes <- hmData %>% 
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n_genes = sum(present), .groups = "drop") %>%
    dplyr::filter(n_genes > 0) %>%
    dplyr::arrange(desc(n_genes)) %>%
    dplyr::pull(cluster)
  
  #### Final filtering:
  hmData <- hmData %>% 
    dplyr::filter(cluster %in% clustersWithGenes) %>% 
    dplyr::mutate(gene = factor(gene, levels = rev(unique(gene))),
                  cluster = factor(cluster, levels = clustersWithGenes))
  
  
  
  ## ---- 3/ setup plot parameters:
  logFCValues <- hmData$log2FC[!is.na(hmData$log2FC)]
  logFCValues <- c(floor(min(logFCValues)), ceiling(max(logFCValues)))
  
  #### set fill gradient colours!
  if ((sign(logFCValues[1]) * sign(logFCValues[2])) == -1) {
    gradient_logFC <- scale_fill_gradientn(
      colours = c("#050","#fff","#fff","#500"),
      values = scales::rescale(c(min(hmData$log2FC), -0.3, 0.3, max(hmData$log2FC))),
      limits = c(logFCValues[1], logFCValues[2]), 
      na.value = plt_gp$fillEmptyValues
    )
  } else {
    gradient_logFC <- scale_fill_gradient(
      low = ifelse(sign(logFCValues[1]) == -1, "#050", "#ffa"),
      high = ifelse(sign(logFCValues[1]) == -1, "#ffa", "#500"), 
      na.value = plt_gp$fillEmptyValues
    )
  }
  
  
  
  ## ---- 4/ plot:
  p <- hmData %>% 
    ggplot(aes(x= cluster, y= gene, fill = factor(present)))+
    geom_tile(color = plt_gp$tilesLineColour, linewidth = 1.4)+
    theme_bw()+
    labs(y= "", x= "", fill= "")+
    guides(fill= "none")+
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, plt_gp$xAxisMaxCharSize),
                     guide = guide_axis(angle = plt_gp$xAxisAngle))+
    scale_fill_manual(values = c("0" = plt_gp$fillEmptyValues, "1" = plt_gp$fillPresence))+
    theme(
      axis.text = element_text(size = plt_gp$axisTextSize, 
                               colour = plt_gp$axisTextColour, 
                               face = "bold"),
      plot.background = element_rect(colour = "white", fill = "white"),
      legend.background = element_rect(colour = "white", fill = "white")
    )
  
  #### case with logFC values
  if (fillByLogFC) {
    p <- hmData %>% 
      ggplot(aes(x= cluster, y= gene, fill = log2FC))+
      geom_tile(color = plt_gp$tilesLineColour, linewidth = 1.4)+
      theme_bw()+
      labs(y= "", x= "", fill= "Log2FC")+
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, plt_gp$xAxisMaxCharSize),
                       guide = guide_axis(angle = plt_gp$xAxisAngle))+
      gradient_logFC+
      theme(
        axis.text = element_text(size = plt_gp$axisTextSize, 
                                 colour = plt_gp$axisTextColour, 
                                 face = "bold"),
        legend.title = element_text(size = plt_gp$legendTitleSize,
                                    colour = "darkcyan",
                                    face = "bold",
                                    margin = margin(b=0.15, unit = "in")),
        plot.background = element_rect(fill = "whitesmoke", colour = "whitesmoke"),
        legend.background = element_rect(fill = "whitesmoke", colour = "whitesmoke")
      )
  }
  
  
  return(p)
  
}





## --------------------------------------------------- 9. Network Gene-Clusters

Net_clust_plot <- function(clustData,
                           colnameCluster = "cluster",
                           colnameGene = "geneID",
                           deaData,
                           colnameGeneDea = "Gene_name",
                           colnameLogFC = "log2FoldChange",
                           net.layout = "nicely",
                           edge.param.color = NULL,
                           edge.param.width = .8,
                           edge.param.alpha = .5,
                           edge.legend.width = 1.6,
                           edge.legend.position = "bottom",
                           edge.legend.direction = "horizontal",
                           edge.legend.layers = 2,
                           edge.legend.fontsize = 11,
                           node.param.term.color = "black",
                           node.param.term.size = c(8,4),
                           node.param.gene.logFC = T,
                           node.param.gene.color = "darkcyan",
                           node.param.gene.size = 3,
                           node.param.gene.alpha = .7,
                           label.genes = T,
                           label.gene.color = c("red","purple","black"),
                           label.gene.size = 2,
                           label.gene.alpha = 0.9,
                           label.gene.repel = T){
  
  
  ## ---- 01/ Debuging :
  if (!all(c(colnameCluster, colnameGene) %in% colnames(clustData))) {
    stop("clustData must contain the 2 mentionned columns! At least 1 is missing!")
  }
  
  if (!all(c(colnameGeneDea, colnameLogFC) %in% colnames(deaData))) {
    stop("deaData must contain the 2 mentionned columns! At least 1 is missing!")
  }
  
  if(!(net.layout %in% c("fr", "kk", "stress", "tree", "mds", "circle", "nicely"))) {
    stop(paste("Please choose a layout within this list : ",
               "fr, kk, stress, tree, mds, circle, gem, drl, nicely",
               sep = "\n"))
  }
  
  
  
  
  
  ## ---- 02/ Building the data :
  clustData <- clustData %>% 
    dplyr::select(Cluster = colnameCluster,
                  Gene = colnameGene)
  
  deaData <- deaData %>% 
    dplyr::select(Gene = colnameGeneDea,
                  logFC = colnameLogFC)
  
  ##
  dat <- clustData %>% 
    tidyr::separate_rows(Gene, sep = "/") %>% 
    dplyr::left_join(deaData, by = "Gene") %>% 
    dplyr::mutate(Regulation = case_when(logFC > 0 ~ "Up",
                                         logFC < 0 ~ "Down"))
  
  ##
  gene_term_map <- dat %>% 
    dplyr::select(Gene, Cluster) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(Gene) %>% 
    dplyr::summarise(n_term = n(),
                     term_for_color = ifelse(n_term == 1, Cluster, NA_character_),
                     .groups = "drop") %>% 
    dplyr::mutate(color_group = case_when(
      n_term == 1 ~ term_for_color,
      n_term == 2 ~ "shared by two",
      n_term >= 3 ~ "shared by more"
    )) %>% 
    dplyr::select(Gene, color_group)
  
  ##
  genes_info <- dat %>% 
    dplyr::select(Gene, logFC, Regulation) %>%
    dplyr::distinct(Gene, .keep_all = TRUE) %>% 
    dplyr::left_join(gene_term_map, by = "Gene")
  
  
  
  
  
  ## ---- 03/ Building igraph object :
  graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Cluster")], directed = FALSE)
  
  ## Annotating nodes
  V(graph_data)$type <- ifelse(V(graph_data)$name %in% genes_info$Gene, "Gene", "Cluster")
  V(graph_data)$genes <- case_when(V(graph_data)$type == "Gene" ~ V(graph_data)$name, TRUE ~ NA)
  V(graph_data)$clusters <- case_when(V(graph_data)$type == "Cluster" ~ V(graph_data)$name, TRUE ~ NA)
  V(graph_data)$logFC <- genes_info$logFC[match(V(graph_data)$name, genes_info$Gene)]
  V(graph_data)$reg <- genes_info$Regulation[match(V(graph_data)$name, genes_info$Gene)]
  V(graph_data)$lab_color <- genes_info$color_group[match(V(graph_data)$name, genes_info$Gene)]
  
  ## Annotating edges
  E(graph_data)$cluster <- dat$Cluster[match(
    paste0(dat$Gene, dat$Cluster),
    paste0(ends(graph_data, es = E(graph_data))[, 1],
           ends(graph_data, es = E(graph_data))[, 2])
  )]
  
  ## Final igraph data
  graph_tbl <- tidygraph::as_tbl_graph(graph_data)
  
  
  
  
  
  ## ---- 04/ Plotting parameters :
  MyPalette <- c("#070","#700","#007","#770","#077","#707","#2a0",
                 "#a30","#40a","#ca0","#05c","#b04","#000","#f0f")
  
  clusters <- sort(unique(dat$Cluster))
  n_clust  <- length(clusters)
  
  if (is.null(edge.param.color)) {
    if (n_clust <= length(MyPalette)) {
      e.colors <- setNames(MyPalette[seq_len(n_clust)], clusters)
    } else {
      e.colors <- setNames(
        c(MyPalette, sample(colors(), n_clust - length(MyPalette))),
        clusters
      )
      message("Palette extended with random colors.")
    }
  } else {
    if (length(edge.param.color) < n_clust) {
      stop("edge.param.color does not contain enough colors for the number of clusters.")
    }
    e.colors <- setNames(edge.param.color[seq_len(n_clust)], clusters)
  }
  
  
  ## Set edge legend guides
  edge_cols <- switch(
    edge.legend.direction,
    "horizontal" = guide_legend(title = NULL,
                                override.aes = list(edge_width = edge.legend.width),
                                position = edge.legend.position,
                                direction = "horizontal",
                                nrow = edge.legend.layers),
    "vertical" = guide_legend(title = NULL,
                              override.aes = list(edge_width = edge.legend.width),
                              position = edge.legend.position,
                              direction = "vertical",
                              ncol = edge.legend.layers)
  )
  
  ## Set logFC values
  logFC_values <- c(floor(min(genes_info$logFC)), ceiling(max(genes_info$logFC)))
  
  ## Set logFC colour gradient
  if((sign(logFC_values[1]) * sign(logFC_values[2])) == -1){
    gradient_logFC <- scale_colour_gradientn(
      colours = c("#050","#fff","#fff","#500"),
      values = scales::rescale(c(min(genes_info$logFC), -0.3, 0.3, max(genes_info$logFC))),
      limits = c(logFC_values[1], logFC_values[2])
    )
  } else {
    gradient_logFC <- scale_colour_gradient(
      low = ifelse(sign(logFC_values[1]) == -1, "#050", "#ffa"),
      high = ifelse(sign(logFC_values[1]) == -1, "#ffa", "#500")
    )
  }
  
  ## Set cluster nodes layers size
  node_term <- list(col1 = node.param.term.color[1],
                    col2 = ifelse(length(node.param.term.color) >= 2,
                                  node.param.term.color[2],
                                  node.param.term.color[1]),
                    size1 = node.param.term.size[1],
                    size2 = ifelse(length(node.param.term.size) >= 2,
                                   node.param.term.size[2],
                                   node.param.term.size[1]))
  
  ## Set label genes colours
  gene_lab_cols <- list(`shared by more` = label.gene.color[1],
                        `shared by two` = ifelse(length(label.gene.color) >= 2,
                                                 label.gene.color[2],
                                                 label.gene.color[1]),
                        distinct = ifelse(length(label.gene.color) >= 3,
                                          label.gene.color[3],
                                          "black"))
  
  
  
  
  ## ---- 05/ plotting :
  g <- ggraph(graph_tbl, layout = net.layout)+
    theme_void()
  
  ## Plot edges
  g <- g +
    geom_edge_link(aes(color = cluster), alpha = edge.param.alpha, edge_width = edge.param.width)+
    scale_edge_color_manual(values = e.colors)+
    guides(edge_color = edge_cols)
  
  ## Plot gene nodes
  if(isFALSE(node.param.gene.logFC)){
    g <- g +
      geom_node_point(colour = node.param.gene.color,
                      alpha = node.param.gene.alpha,
                      size = ifelse(V(graph_tbl)$type == "Gene", node.param.gene.size, 0))
  } else {
    g <- g +
      geom_node_point(aes(colour = logFC),
                      alpha = node.param.gene.alpha,
                      size = ifelse(V(graph_tbl)$type == "Gene", node.param.gene.size, 0))+
      gradient_logFC+
      guides(colour = guide_colourbar(title = "log2FC",
                                      barwidth = unit(.6, "cm"),
                                      barheight = unit(3, "cm")))
  }
  
  ## Plot cluster nodes
  g <- g +
    #### first layer
    geom_node_point(colour = ifelse(V(graph_tbl)$type == "Cluster", node_term$col1, "white"),
                    size = ifelse(V(graph_tbl)$type == "Cluster", node_term$size1, 0),
                    alpha = ifelse(V(graph_tbl)$type == "Cluster", 0.3, 0),)+
    #### second layer
    geom_node_point(colour = ifelse(V(graph_tbl)$type == "Cluster", node_term$col2, "white"),
                    size = ifelse(V(graph_tbl)$type == "Cluster", node_term$size2, 0),
                    alpha = ifelse(V(graph_tbl)$type == "Cluster", 1, 0),)
  
  ## Plot gene labels
  if(isTRUE(label.genes)){
    g <- g +
      geom_node_text(aes(label = genes),
                     show.legend = F,
                     fontface = "bold.italic",
                     size = label.gene.size,
                     alpha = label.gene.alpha,
                     repel = label.gene.repel,
                     colour = ifelse(V(graph_tbl)$lab_color == "shared by more",
                                     gene_lab_cols$`shared by more`,
                                     ifelse(V(graph_tbl)$lab_color == "shared by two",
                                            gene_lab_cols$`shared by two`,
                                            gene_lab_cols$distinct)))
  }
  
  ## Setup plot theme
  g <- g +
    theme(
      #### Backgrounds
      plot.background = element_rect(colour = "white", fill = "white"),
      legend.background = element_rect(colour = "white", fill = "white"),
      #### Plot
      plot.margin = margin(b=0.1,t=0.1,r=0.1,l=0.1, unit = "in"),
      plot.caption = element_text(size = 11,
                                  face = "bold.italic",
                                  colour = "#333",
                                  margin = margin(t=0.25, unit = "in")),
      #### Legends
      legend.title = element_text(size = 14,
                                  face = "bold",
                                  colour = "darkcyan",
                                  margin = margin(b=0.15, unit = "in")),
      legend.text = element_text(size = edge.legend.fontsize,
                                 face = "bold"),
      legend.margin = margin(l=0.2,b=0.2,r=0.1,t=0.2, unit = "in")
    )
  
  ## Setup caption for labels colour info
  shared <- length(genes_info$color_group[genes_info$color_group == "shared by more"])
  shared2 <- length(genes_info$color_group[genes_info$color_group == "shared by two"])
  notshared <- nrow(genes_info)-(shared+shared2)
  
  g <- g +
    labs(caption = paste(
      sprintf("%-30s: %s", 
              paste0("Distinct genes (", gene_lab_cols$distinct,")"),
              notshared),
      sprintf("%-30s: %s", 
              paste0("Genes shared by two (", gene_lab_cols$`shared by two`,")"),
              shared2),
      sprintf("%-30s: %s", 
              paste0("Genes shared by many (", gene_lab_cols$`shared by more`,")"),
              shared),
      sep = "\n"
    ))
  
  return(list(netdata = dat, geneInfoData = genes_info, plot = g))
}






## --------------------------------------------------- 10. Get GSEA gene sets with gmt-files

getGmtPaths <- function(gmtFile, genes){
  gmt <- gmtPathways(gmtFile)
  
  tmp <- unique(unlist(gmt)) %>% 
    intersect(genes)
  
  for(i in 1:length(gmt)){
    gmt[[i]] <- intersect(gmt[[i]], tmp)
  }
  
  return(gmt)
}








## ------------------------------------------------------- 11. GSEA plotEnrichment

gseaEnrichPlot <- function(gseaResults,
                           rankedGenes,
                           pathwayName,
                           pathwayGenes,
                           genesToHighlight = NULL,
                           segmentsSizeScale = 0.07,
                           segmentLineWidth = 0.5,
                           lineWidth = 0.05,
                           lineColour = "green"){
  
  
  ## ---- 1/ Debuging :
  if(!is.null(genesToHighlight)){
    if(!all(genesToHighlight %in% names(rankedGenes))){
      stop("genesToHighlight and rankedGenes don't match!")
    }
  }
  
  if(!(pathwayName %in% gseaResults$pathway)){
    stop(paste0(pathwayName, " doesn't exist in the results!"))
  }
  
  
  
  ## ---- 2/ Prepare the data :
  ## 01.
  dat <- data.frame(
    rank_pos = seq_along(rankedGenes),
    rank_val = rankedGenes,
    genes = names(rankedGenes),
    genes_highlighted = NA_character_
  )
  
  ## 02.
  pw_idx <- match(pathwayGenes, names(rankedGenes))
  pw_idx <- pw_idx[!is.na(pw_idx)] %>% sort()
  gseaStats <- fgsea:::calcGseaStat(
    stats = rankedGenes,
    selectedStats = pw_idx,
    returnAllExtremes = T
  )
  
  tops <- gseaStats$tops
  bottoms <- gseaStats$bottoms
  
  diff <- (max(tops) - min(bottoms)) 
  
  ## 03.
  N <- length(rankedGenes)
  Xs <- as.vector(rbind(pw_idx - 1, pw_idx))
  Ys <- as.vector(rbind(bottoms, tops))
  
  toPlot <- data.frame(x=c(0, Xs, N + 1), 
                       y=c(0, Ys, 0))
  
  
  ## 04.
  df_segments <- data.frame(
    g = names(rankedGenes)[pw_idx],
    idx = pw_idx
  ) %>% 
    dplyr::mutate(
      deg = case_when(g %in% genesToHighlight ~ "is Deg",
                      TRUE ~ "not Deg")
    )
  
  
  ## 05. Plot
  g <- toPlot %>% ggplot(aes(x = x, y = y)) + 
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(colour = "darkcyan",
                                    margin = margin(b=0.2, unit = "in"),
                                    face = "bold"),
          axis.title = element_text(face = "bold",
                                    colour = "darkcyan"),
          axis.title.x = element_text(margin = margin(t=0.15, unit = "in")),
          axis.title.y = element_text(margin = margin(r=0.15, unit = "in"))) +
    labs(x = "Rank", 
         y = "Enrichment score", 
         title = pathwayName, 
         colour = "", 
         alpha = "")
  
  g <- g +
    geom_hline(yintercept = 0, colour = "gray45", linetype = "dashed")+
    geom_hline(yintercept = max(tops), color = "red", linetype = "dashed")+
    geom_hline(yintercept = min(bottoms), color = "red", linetype = "dashed")+
    geom_line(linewidth = lineWidth, 
              colour = lineColour,
              alpha = 0.7) +
    geom_segment(data = df_segments,
                 mapping = aes(x = idx, 
                               y = -diff*segmentsSizeScale, 
                               xend = idx, 
                               yend = diff*segmentsSizeScale,
                               colour = deg,
                               alpha = deg),
                 linewidth = segmentLineWidth)
  
  if(is.null(genesToHighlight)){
    g <- g +
      guides(colour = "none",
             alpha = "none")
    
  } else {
    g <- g +
      scale_colour_manual(values = c("is Deg" = "darkorange",
                                     "not Deg" = "gray20"))+
      scale_alpha_manual(values = c("is Deg" = 1,
                                    "not Deg" = 0.3)) +
      guides(colour = guide_legend(override.aes = list(linewidth = 1.2)),
             alpha = guide_legend(override.aes = list(linewidth = 1.2)))
    
  }
  
  
  
  
  
  return(g)
  
}










