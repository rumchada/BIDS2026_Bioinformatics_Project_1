# Project 1 Part 1 and Part 3 Functions
#QC functions and dim_reduction functions.

pc_var_association <- function(bulk_ds){
  library(pheatmap)
  vst_mat <- assay(bulk_ds, "var_stable")
  #Select top variable gene
  top_var_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), 2000)
  vst_top <- vst_mat[top_var_genes, ]
  # Run PCA via R
  # Remember PCA is taking the right unitary martrix (VT) (feature space) of SVD
  # VT tells us the principal directions of the data set
  # multiply the the original matrix by the right unitary which is equivalent to the
  # U time signma (left unitary) * Singular values (Size of the PCs)
  pca_res <- prcomp(t(vst_top))
  
  top_pcs <- pca_res$x[, 1:10]
  
  #creating a list of possible covariate to compare with PCs
  meta_data <- as.data.frame(bulk_ds@colData)
  
  
  covars <- split(meta_data, 
                  rep(1:ncol(meta_data),
                      each = nrow(meta_data)))
  
  covars <- lapply(seq_len(ncol(meta_data)),
                   function(x) meta_data[ , x])
  
  names(covars) <- colnames(meta_data)
  
  
  
  
  get_p_value <- function(pc_vector, var) {
    
    valid <- !is.na(pc_vector) & !is.na(var)
    pc_vector <- pc_vector[valid]
    var_vector <- var[valid]
    
    
    if (is.numeric(var_vector)) {
      # Continuous variable: use linear regression/correlation test
      return(summary(lm(pc_vector ~ var_vector))$coefficients[2, 4])
    } else {
      # Categorical variable: use ANOVA
      fit <- aov(pc_vector ~ as.factor(var_vector))
      
      anova_sum <- summary(fit)[[1]]
      p_val <- anova_sum[1, "Pr(>F)"]
      
      if (is.null(p_val) || length(p_val) != 1) {
        return(NA)
      }
      return(p_val)
    }
  }
  
  
  p_matrix <- matrix(NA, 
                     nrow = length(covars), 
                     ncol = ncol(top_pcs),
                     dimnames = list(names(covars), colnames(top_pcs)))
  
  
  for(vars in names(covars)){
    for(pcs in colnames(top_pcs)){
      
      current_var <- covars[[vars]]
      current_pc <- top_pcs[, pcs]
      
      
      p_matrix[vars, pcs] <- get_p_value(current_pc, current_var)
    }
  }
  
  log_p_matrix <- -log10(p_matrix) %>% as.data.frame() %>% na.omit()
  
  # 5. Plot the heatmap
  p_val_map <-pheatmap(log_p_matrix,  cluster_rows=F, cluster_cols=F, main = "Correlation between of MetaVar X within PCX") 
  return(p_val_map)
}

baseline_dimreduction <- function(bulk_dataset, 
                                  grouping = NA){
  
  library(umap)
  library(ggplot2)
  
  if (!is.na(grouping) && is.character(grouping)) {
    if (!grouping %in% colnames(colData(bulk_dataset))) {
      stop(paste("The grouping variable", grouping, "was not found in your colData."))
    }
    anno_variable <- c("condition", grouping)
    color_variable <- grouping
  } else {
    anno_variable <- c("condition")
    color_variable <- "condition"
  }
  
  
  #prioritize variance stable assay within the bulk_dataset
  assays(bulk_dataset) <- assays(bulk_dataset)[c("var_stable", "counts", "log_counts")]
  # takes the euclidan distance of a variance stabilized dataset
  # row by row euclidean distance against each other
  # plots the visual distance between each sample onto heatmap
  distance_plot <- plot_sample_clustering(bulk_dataset, anno_vars = anno_variable, distance = "euclidean")
  
  
  #____Extracting Top Variable Gene per sample -----#
  # Setting Variance Stable Counts as Priority Assay
  vst_mat <- assay(bulk_dataset, "var_stable")
  #Select top variable gene
  top_var_genes <- head(order(rowVars(vst_mat), decreasing = TRUE), 2000)
  vst_top <- vst_mat[top_var_genes, ]
  # Run PCA via R
  # Remember PCA is taking the right unitary martrix (VT) (feature space) of SVD
  # VT tells us the principal directions of the data set
  # multiply the the original matrix by the right unitary which is equivalent to the
  # U time signma (left unitary) * Singular values (Size of the PCs)
  pca_res <- prcomp(t(vst_top))
  
  
  # Plotting first two PCAs
  pca_plot <- plot_pca(bulk_dataset, PC_x = 1, PC_y = 2, color_by = color_variable)
  
  
  
  # Calculating Percentage of Variance Per PC for the Scree Plot
  
  # Calculate percentage variance explained
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
  
  # Create a data frame for plotting (first 10 PCs)
  scree_data <- data.frame(
    PC = paste0("PC", 1:10),
    Variance = var_explained[1:10]
  )
  
  # Ensure PC is treated as a factor in order
  scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)
  
  # Scree Plot ensures that we are plotting PCs that capture the most amount of variability within our data
  scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance, group = 1)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(title = "Scree Plot: Variance Explained per PC",
         x = "Principal Component",
         y = "Percentage of Variance (%)") +
    theme_minimal()
  
  
  #Input the Left Unitary * Singular Values
  # 1. Run PCA
  pca_res <- prcomp(t(assay(bulk_dataset, "var_stable")), rank. = 50)
  
  # 2. Use the 'x' slot (the PC scores) as input for UMAP
  # We only take the number of PCs identified by your scree plot (e.g., 1:15)
  umap_out <- umap(pca_res$x[, 1:15])
  
  
  umap_df <- data.frame(
    Sample_IDS = rownames(pca_res$x),
    UMAP_1 = umap_out$layout[, 1],
    UMAP_2 = umap_out$layout[, 2],
    ColorGroup = as.factor(colData(bulk_dataset)[[color_variable]]))
  
  umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = ColorGroup)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "UMAP of COVID-19 Samples", subtitle = "Input: Top 15 PCs from VST Data")
  
  
  
  dimreduction_visuals <- c(distance_plot, scree_plot, pca_plot$plot, umap_plot)
  return(dimreduction_visuals)
}



# Part 2 Intra-Dataset Differntial Expression Functions
edgeR_diffexp <- function(dds_object, condition_col, ref_group_name) {
  require(edgeR)
  #For each condition combination, runs an indepedent differential expression analysis
  #0-Intercept
  
  #Converting to a DGEList
  raw_counts <- dds_object@assays@data$counts
  metadata <- dds_object@colData
  y <- DGEList(counts = raw_counts , group = metadata$condition)
  
  #Filter by low counts
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  #Performs TMM Normalization, Corrects for systemic biases between samples.
  y <- calcNormFactors(y)
  
  # reaffirm healthy as level reference
  y$samples$condition <- relevel(
    factor(y$samples$group),
    ref = ref_group_name
  )
  
  # ~0+ not to include an intercept column and instead to include a column for each group
  # Always use a model intecept term because then the data will no fit to the glm line as closely
  # A model without an intercept term would only be recommended in cases where there is a strong biological reason why a zero covariate should be associated with a zero expression value (Law, et.al., F10000Research, 2020, PMID: 33604029)
  
  design <- model.matrix(~0+condition, data = y$samples)
  conds <- gsub("^condition", "", colnames(design))
  
  contrast_strings <- combn(conds, 2, FUN = function(x) paste0(x[1], "-", x[2]))
  
  
  contrast_strings
  colnames(design) <- gsub(glue("^{condition_col}"), "", colnames(design))
  
  my.contrasts <- makeContrasts(
    contrasts = contrast_strings,
    levels = colnames(design))
  
  #Test Using a GLM
  results <- list()
  
  # Run Dispersion on the Data
  y <- estimateDisp(y, design)
  
  #fit model once
  fit <- glmQLFit(y, design)
  
  # run the test for each contrast
  for(i in seq_len(ncol(my.contrasts))) {
    
    qlf <- glmQLFTest( fit, contrast = my.contrasts[, i])
    
    results[[colnames(my.contrasts)[i]]] <-
      topTags(qlf, n = nrow(y), p.value = 0.05, adjust.method = 'bonferroni')$table %>%
      as.data.frame() %>%
      tibble::rownames_to_column("geneid") %>%
      dplyr::rename_with( ~ "log2foldchange", .cols = "logFC") %>%
      dplyr::rename_with( ~ "p.val_adj", .cols = "FWER" )
  }
  
  return(results)
}


# Volcano Plot Visualizations
volcano_plot <- function(diffexp_df, 
                         #data-alterations
                         unique = FALSE,
                         log2fc_thresh = 0,
                         p.val_thresh = 0.05,
                         #visual parts
                         lower_xlim = -5, 
                         upper_xlim = 5, 
                         step = 1) {
  
  library(ggrepel)
  library(ggplot2)
  
volcano_df <- diffexp_df %>%
    mutate(
      color = case_when(
        p.val_adj < p.val_thresh & log2foldchange >  log2fc_thresh  ~ "Upregulated",
        p.val_adj < p.val_thresh & log2foldchange < -log2fc_thresh  ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  
  volplot <- ggplot(volcano_df) +
    
    aes(x = log2foldchange,y = -log10(p.val_adj), color = color, label = geneid) +
    
    geom_point(size = 2) +
    
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray",linetype = "dashed"
    ) +
    
    geom_hline(yintercept = -log10(0.05), col = "gray",linetype = "dashed") +
    
    theme_classic(base_size = 12) +
    
    theme(axis.title.y = element_text(margin = margin(0, 20, 0, 0),size = rel(1.1),color = "black"
    ),
    
    axis.title.x = element_text(hjust = 0.5, margin = margin(20, 0, 0, 0),size = rel(1.1),
                                color = "black"
    ),
    
    plot.title = element_text(hjust = 0.5)
    ) +
    
    scale_color_manual(values = c("Upregulated" = "red","Downregulated" = "blue","Not Significant" = "gray"
    )
    ) +
    
    geom_text(nudge_y = 0.5, check_overlap = TRUE) +
    
    coord_cartesian(xlim = c(lower_xlim, upper_xlim)) +
    
    scale_x_continuous(
      breaks = seq(lower_xlim, upper_xlim, step)
    )
  
  return(list(volcano_df, volplot))
}


#unpacks enrichGO results, converts them to the user's desired gene IDs, then aligns them to their log2fc results

enrichgo_unpack_ver2 <- function(intial_table = filtered_results$`healthy-covid19`, 
                                 key_table = enrichgo_results$`healthy-covid19`,
                                 chosen_desc,
                                 #assume user does not want to convert geneids
                                 #arguments below are for the conversion helperfunction
                                 convert = TRUE,
                                 from_type = "ensembl",
                                 to_type = "symbol",
                                 ensembl_dataset = "hsapiens_gene_ensembl"){
  
  
  
  # Loop will be placed here.
  cluster <- tibble(key@result) %>%
    dplyr::filter(str_detect(Description, chosen_desc)) %>%
    dplyr::mutate(Description = as.character(Description))
  cluster
  
  
  ### helper functionf or converting geneidA to geneidB
  ### user
  gene_id_converter <- function(vector, from_type, to_type, ensembl_dataset){
    library(biomaRt)
    # geneid pulling for correct attribute
    id_map <- c(
      ensembl = "ensembl_gene_id",
      entrez  = "entrezgene_id",
      symbol  = "external_gene_name"
    )
    # p
    from_attr <- id_map[[from_type]]
    to_attr <-id_map[[to_type]]
    
    ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)
    
    return_df <- getBM( attributes = c(from_attr, to_attr), filters = from_attr, values = vector, mart = ensembl
    )
    return(return_df)
    
  }
  
  ###creates a helper function to unlist and strsplits the enriched Gene IDs per term
  unlist_converter <- function(str) {
    if (grepl("^ENSG|\\d", str)){
      bruh <- as.character(unlist(strsplit(str, "/")))
    }else{
      bruh <- as.integer(unlist(strsplit(str, "/")))
    }
    return(bruh)
  }
  
  #apply this function to the geneID column
  extract_result <- lapply(cluster$geneID, unlist_converter)
  #Next check the length of the new vector to see if it matches the Count column in the GSEA DF.
  length(extract_result[[1]])
  
  #create an empty vector
  path_vector <- c()
  #for each index in the first column of numbers
  for (i in seq_along(cluster$Description)) {
    #replicate the name of the Description for the length of the GO term in column 1
    path_vector <- c(path_vector, rep(cluster$Description[i], length(extract_result[[i]])))
  }
  
  # make new tibble that pairs the name of the Pathway from the GOTerm and enriched Entrez ID from the "post-processed" GeneID column
  pathway_tibble <- tibble(gene_desc = path_vector, geneid = unlist(extract_result))
  pathway_tibble
  
  if(convert == FALSE){
    pathway_genenames <- pathway_tibble
  }else{
    pathway_genenames <- gene_id_converter(pathway_tibble$geneid, 
                                           from_type = from_type,
                                           to_type = to_type,
                                           ensembl_dataset = ensembl_dataset)
    
    new_colname <- paste0(to_type,"_", "geneid")
    pathway_genenames <- pathway_genenames %>%
      rename(geneid = colnames(pathway_genenames)[1]) %>%
      rename(new_colname = colnames(pathway_genenames)[2])
  }
  
  
  final_path_tibble <- inner_join(pathway_tibble, pathway_genenames, by = "geneid")
  final_path_tibble <- final_path_tibble %>% mutate(geneid = as.character(geneid))
  final_path_tibble
  
  
  final_path_fc_tibble <- inner_join(final_path_tibble, initial_table, by = "geneid")
  
  fpfc_tibble <- final_path_fc_tibble %>%
    relocate(new_colname, .after = geneid)
  
  return(fpfc_tibble)
}


### helper functionf or converting geneidA to geneidB
### user
gene_id_converter <- function(vector, from_type, to_type, ensembl_dataset){
  library(biomaRt)
  # geneid pulling for correct attribute
  id_map <- c(
    ensembl = "ensembl_gene_id",
    entrez  = "entrezgene_id",
    symbol  = "external_gene_name"
  )
  # p
  from_attr <- id_map[[from_type]]
  to_attr <-id_map[[to_type]]
  
  ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)
  
  return_df <- getBM( attributes = c(from_attr, to_attr), filters = from_attr, values = vector, mart = ensembl
  )
  return(return_df)
  
}

heatmap_function <- function(initial_table, table_list) {
  
  require(colorRamp2)
  require(ComplexHeatmap)
  
  # Initialize the list to store heatmap objects
  heatmap_list <- list()
  
  # Extract the raw count matrix once outside the loop for efficiency
  # Assumes a SingleCellExperiment or SummarizedExperiment-like structure
  expr_raw <- as.matrix(initial_table@assays@data$counts)
  
  # Loop through the list of filtered differential expression results
  for (i in seq_along(table_list)) {
    
    # Grab the ith table
    filtered_diff_exp_results <- table_list[[i]]
    
    # Split the names of the ctrl-treatment comparison
    comparison_name <- names(table_list)[i]
    ctrl_treatment <- strsplit(comparison_name, "-")[[1]]
    #itemizing the names of the
    ctrl_name      <- ctrl_treatment[1]
    treatment_name <- ctrl_treatment[2]
    
    # Message tracking
    message(glue::glue("Processing: {ctrl_name} vs {treatment_name}"))
    
    # Pull the differentially expressed gene IDs
    diffexp_geneids <- as.character(filtered_diff_exp_results$geneid)
    
    # Filter the matrix for the target genes and matching sample columns
    # Note: Added anchors to ensure exact suffix matching
    col_pattern <- glue::glue("_{ctrl_name}$|_{treatment_name}$")
    matching_cols <- grepl(col_pattern, colnames(expr_raw))
    
    expr <- expr_raw[diffexp_geneids, matching_cols, drop = FALSE]
    
    # Scaling and z-score normalizing the counts matrix
    expr_z <- t(scale(t(expr)))
    expr_z[is.na(expr_z)] <- 0
    rownames(expr_z) <- rownames(expr)
    colnames(expr_z) <- colnames(expr)
    
    # Gene ID conversion
    conversion_table <- gene_id_converter(
      rownames(expr_z), 
      "ensembl", 
      "symbol", 
      "hsapiens_gene_ensembl"
    )
    
    gene_map <- setNames(
      conversion_table$external_gene_name,
      conversion_table$ensembl_gene_id
    )
    
    new_names <- gene_map[rownames(expr_z)]
    
    # Fallback to Ensembl ID if no symbol is found
    rownames(expr_z) <- ifelse(is.na(new_names), rownames(expr_z), new_names)
    
    # Clustering data (Optional: ComplexHeatmap does this natively if cluster_rows = TRUE, 
    # but kept if you intend to pass hclust objects explicitly)
    gene_dist <- dist(expr_z)
    gene_hclust <- hclust(gene_dist)
    
    # Define color palette
    col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    
    # Generate Heatmap object
    p <- Heatmap(
      expr_z, 
      name            = "Z-score",
      col             = col_fun,
      cluster_columns = TRUE,     
      cluster_rows    = TRUE, # To use your manual cluster, change to: cluster_rows = gene_hclust     
      row_names_gp    = gpar(fontsize = 12),
      column_title    = glue::glue("{ctrl_name} - {treatment_name} Diff Exp Genes")
    )
    
    # Store the plot object into the initialized list
    heatmap_list[[comparison_name]] <- p
  }
  
  return(heatmap_list)
}




