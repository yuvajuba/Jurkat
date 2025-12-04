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
    theme(plot.background = element_rect(fill = "whitesmoke", color = "whitesmoke"))
  
  return(list(plot = p, intersections = common_genes, group_specific = specific_genes))
}



## ---------------------------------------------- 5. Network plot
Net_plot <- function(GOobject,
                     DEAres,
                     condition.name = "condition",
                     go.obj.colnames = list(term = "Description",
                                            gene = "geneID",
                                            pval = "p.adjust"),
                     dea.res.colnames = list(gene = "Gene_name",
                                             logFC = "log2FoldChange",
                                             pval = "padj"),
                     select.terms = NULL,
                     net.layout = "nicely",
                     edge.params = list(color = NULL,
                                        width = .8,
                                        alpha = .5),
                     node.params = list(term.col = "black",
                                        term.size = c(8,4),
                                        gene.col = "logFC",
                                        gene.shape = F,
                                        gene.size = 3,
                                        gene.alpha = .7),
                     label.params = list(gene.col = c("navy", "black"),
                                         gene.size = 2)){
  
  ## parameters :
  net.layout.options <- c(
    "fr", "kk", "stress", "tree", "mds", "circle", "gem", "drl", "nicely"
  )
  
  MyPalette <- c("purple", "orangered", "darkgreen", "navy", "darkred", "darkcyan",
                 "gold", "violet", "black", "brown", "khaki3")
  
  
  ## debuging :
  if (!inherits(GOobject, "enrichResult")) {
    stop("You need to provide a GO object from enrichGO or similar.")
  }
  
  if(!(is.data.frame(DEAres))){
    stop("DEAres needs to be a dataframe")
  }
  
  if(!(net.layout %in% net.layout.options)) {
    stop(paste("Please choose a layout within this list : ",
               "fr, kk, stress, tree, mds, circle, gem, drl, nicely",
               sep = "\n"))
  }
  
  #### Build Network data :
  dat <- GOobject@result %>% 
    dplyr::filter(.data[[go.obj.colnames$pval]] < 0.05) %>% 
    dplyr::arrange(p.adjust) %>% 
    dplyr::select(Term = go.obj.colnames$term, 
                  Gene = go.obj.colnames$gene)
  
  default_terms <- dat %>% head(9) %>% dplyr::pull(Term)
  
  dea <- DEAres %>% 
    dplyr::filter(.data[[dea.res.colnames$pval]] < 0.05) %>% 
    dplyr::select(Gene = dea.res.colnames$gene, 
                  logFC = dea.res.colnames$logFC)
  
  dat <- dat %>% 
    dplyr::mutate(Cluster = condition.name)  %>% 
    tidyr::separate_rows(Gene, sep = "/") %>% 
    left_join(dea, by = "Gene")
  
  ### selecting terms
  selected_terms <- if(is.null(select.terms)) default_terms else select.terms
  
  ### updating netdata
  dat <- dat[dat$Term %in% selected_terms,] %>% 
    dplyr::mutate(Regulation = case_when(logFC > 0 ~ "Up",
                                         logFC < 0 ~ "Down"))
  
  ### adding group colour to colour shared and distinct genes
  gene_term_map <- dat %>%
    dplyr::select(Gene, Term) %>%
    distinct() %>%
    group_by(Gene) %>%
    summarise(n_terms = n(),
              term_for_color = ifelse(n_terms == 1, Term, NA_character_),
              .groups = "drop") %>%
    mutate(
      color_group = ifelse(n_terms == 1, term_for_color, "Shared")
    ) %>% 
    dplyr::select(Gene, color_group)
  
  ### creating genes info data
  genes_info <- dat %>%
    dplyr::select(Gene, Cluster, logFC, Regulation) %>%
    dplyr::distinct(Gene, .keep_all = TRUE) %>% 
    dplyr::left_join(gene_term_map, by = "Gene")
  
  
  
  #### built igraph data :
  ## 1. Initiate igraph object
  graph_data <- igraph::graph_from_data_frame(dat[, c("Gene", "Term")], directed = FALSE)
  
  ## 2. Annotating nodes (Verticies "V") [type, genes, terms, logFC, reg]
  V(graph_data)$type <- ifelse(V(graph_data)$name %in% genes_info$Gene, "Gene", "Term")
  V(graph_data)$genes <- case_when(V(graph_data)$type == "Gene" ~ V(graph_data)$name, TRUE ~ NA)
  V(graph_data)$terms <- case_when(V(graph_data)$type == "Term" ~ V(graph_data)$name, TRUE ~ NA)
  V(graph_data)$logFC <- genes_info$logFC[match(V(graph_data)$name, genes_info$Gene)]
  V(graph_data)$reg <- genes_info$Regulation[match(V(graph_data)$name, genes_info$Gene)]
  V(graph_data)$lab_color <- genes_info$color_group[match(V(graph_data)$name, genes_info$Gene)]
  
  ## 3. Annotating edges (Edges "E") [term]
  E(graph_data)$term <- dat$Term[match(
    paste0(dat$Gene, dat$Term),
    paste0(ends(graph_data, es = E(graph_data))[, 1],
           ends(graph_data, es = E(graph_data))[, 2])
  )]
  
  ## 4. Get the igraph table
  graph_tbl <- tidygraph::as_tbl_graph(graph_data)
  
  
  #### ggraph plot :
  ## 1. params :
  nterm <- length(selected_terms)
  
  if(is.null(edge.params$color)){
    if(nterm <= length(MyPalette)){
      e.colors <- MyPalette[1:nterm]
      e.colors <- setNames(e.colors, selected_terms)
    } else {
      e.colors <- sample(colors(), nterm)
      e.colors <- setNames(e.colors, selected_terms)
    }
    
  } else {
    e.colors <- edge.params$color[1:nterm]
    e.colors <- setNames(e.colors, selected_terms)
  }
  
  ## 2. plot :
  ## 2.1 Basic + edges
  g <- ggraph(graph_tbl, layout = net.layout)+
    theme_void()+
    geom_edge_link(aes(color = term), alpha = edge.params$alpha, edge_width = edge.params$width)+
    scale_edge_color_manual(values = e.colors)
  
  if(nterm > 9){
    g <- g +
      guides(edge_color = guide_legend(title = NULL,
                                       override.aes = list(edge_width = 1.6),
                                       position = "right",
                                       direction = "vertical",
                                       ncol = ifelse(nterm <= 12, 1, 2)))
  } else {
    g <- g +
      guides(edge_color = guide_legend(title = NULL,
                                       override.aes = list(edge_width = 1.6),
                                       position = "bottom",
                                       direction = "horizontal",
                                       nrow = ifelse(nterm <= 3, 1, 
                                                     ifelse(nterm <= 6, 2, 3))))
  }
  
  
  ## 2.2 node - genes
  logFC_values <- c(floor(min(genes_info$logFC)), ceiling(max(genes_info$logFC)))
  
  if(length(node.params$term.size) == 2){
    node.term.size <- node.params$term.size
  } else {
    node.term.size <- c(node.params$term.size, node.params$term.size)
  }
  
  if(node.params$gene.col == "logFC"){
    if(isTRUE(node.params$gene.shape)){
      g <- g +
        geom_node_point(aes(colour = logFC,
                            shape = reg),
                        alpha = node.params$gene.alpha,
                        size = ifelse(V(graph_tbl)$type == "Term", 
                                      node.term.size[2],
                                      node.params$gene.size))+
        scale_colour_gradientn(colours = c("#060", "#fff", "#fff", "#600"),
                               values = scales::rescale(c(logFC_values[1], -0.5, 0.5, logFC_values[2])),
                               limits = c(logFC_values[1], logFC_values[2]))+
        scale_shape_manual(values = c("Down"=19, "Up"=17), na.translate = F)+
        guides(colour = guide_colourbar(title = "log2FC",
                                        barwidth = unit(.6, "cm"),
                                        barheight = unit(3, "cm")),
               shape = guide_legend(title = "Regulated genes",
                                    override.aes = list(size = 3)))
    } else {
      g <- g +
        geom_node_point(aes(colour = logFC),
                        alpha = node.params$gene.alpha,
                        size = ifelse(V(graph_tbl)$type == "Term", 
                                      node.term.size[2],
                                      node.params$gene.size))+
        scale_colour_gradientn(colours = c("#060", "#fff", "#fff", "#600"),
                               values = scales::rescale(c(logFC_values[1], -0.5, 0.5, logFC_values[2])),
                               limits = c(logFC_values[1], logFC_values[2]))+
        guides(colour = guide_colourbar(title = "log2FC",
                                        barwidth = unit(.6, "cm"),
                                        barheight = unit(3, "cm")))
    }
    
  } else {
    if(!(node.params$gene.col %in% colors())){
      if(!(grepl("^#([A-Fa-f0-9]{6})$", node.params$gene.col))){
        stop(paste("Error in node.params$gene.color",
                   "You need to provide a valide color name ",
                   "a hexadecimal code or 'logFC' to color by logFC values",
                   sep = "\n"))
      }
    }
    
    if(isTRUE(node.params$gene.shape)){
      g <- g +
        geom_node_point(aes(shape = reg),
                        colour = node.params$gene.col,
                        alpha = node.params$gene.alpha,
                        size = ifelse(V(graph_tbl)$type == "Term", 
                                      node.term.size[2],
                                      node.params$gene.size))+
        scale_shape_manual(values = c("Down"=19, "Up"=17), na.translate = F)+
        guides(shape = guide_legend(title = "Regulated genes",
                                    override.aes = list(size = 3)))
    } else {
      g <- g +
        geom_node_point(aes(colour = reg),
                        alpha = node.params$gene.alpha,
                        size = ifelse(V(graph_tbl)$type == "Term", 
                                      node.term.size[2],
                                      node.params$gene.size))+
        scale_colour_manual(values = c("black", node.params$gene.col), na.translate = F)+
        guides(colour = guide_legend(title = "Regulated genes",
                                     override.aes = list(size = 3)))
    }
    
    
    
  }
  
  
  ## 2.3 node - terms
  g <- g+
    ## first layer node term
    geom_node_point(colour = ifelse(!(is.na(V(graph_tbl)$terms)), node.params$term.col, "white"),
                    alpha = ifelse(!(is.na(V(graph_tbl)$terms)), .3, 0),
                    size = ifelse(!(is.na(V(graph_tbl)$terms)), node.term.size[1], 0))+
    ## second layer node term
    geom_node_point(colour = ifelse(!(is.na(V(graph_tbl)$terms)), node.params$term.col, "white"),
                    alpha = ifelse(!(is.na(V(graph_tbl)$terms)), 1, 0),
                    size = ifelse(!(is.na(V(graph_tbl)$terms)), node.term.size[2], 0))
  
  
  ## 2.4 node - label genes
  g <- g +
    geom_node_text(aes(label = genes),
                   colour = ifelse(V(graph_tbl)$lab_color == "Shared", 
                                   label.params$gene.col[1], 
                                   label.params$gene.col[2]),
                   show.legend = F,
                   fontface = "bold.italic",
                   repel = T,
                   size = label.params$gene.size)
  
  ## 2.5 plot theme
  g <- g +
    theme(
      legend.margin = margin(t= 0.2, l= 0.2, r=0.1, b=0.1, unit = "in"),
      legend.title = element_text(size = 13, 
                                  face = "bold", 
                                  colour = "darkcyan", 
                                  margin = margin(b=0.1, t=0.2, unit = "in")),
      legend.text = element_text(size = 10,
                                 face = "bold"),
      plot.background = element_rect(colour = "whitesmoke", fill = "whitesmoke"),
      legend.background = element_rect(colour = "whitesmoke", fill = "whitesmoke")
    )
  
  return(list(netdata = dat, plot = g))
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


