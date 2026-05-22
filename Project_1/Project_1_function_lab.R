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


