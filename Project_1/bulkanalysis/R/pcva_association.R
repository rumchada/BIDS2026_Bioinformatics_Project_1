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

  meta_data[] <- lapply(meta_data, factor)

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

      if (length(p_val) != 1 || is.na(p_val) || is.nan(p_val) || is.infinite(p_val)) {
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

    current_var <- covars[[vars]]

    if (is.factor(current_var)) {
      n_levels <- nlevels(current_var)
    } else {
      n_levels <- length(unique(current_var))
    }
    if (n_levels < 2) next

    for (pcs in colnames(top_pcs)) {

      current_pc <- top_pcs[, pcs]

      p_matrix[vars, pcs] <- get_p_value(current_pc, current_var)
    }
  }

  log_p_matrix <- -log10(p_matrix)

  max_finite <- max(log_p_matrix[is.finite(log_p_matrix)])

  log_p_matrix[is.infinite(log_p_matrix)] <- max_finite + 1

  log_p_matrix <- log_p_matrix[!apply(is.na(p_matrix), 1 , all), ]

  # 5. Plot the heatmap
  p_val_map <- pheatmap(log_p_matrix,
                        legend_labels = "log10(pvalue)",
                        cluster_rows=F,
                        cluster_cols=F,
                        main = "PC-Metadata Associations (-log10(p-value))")
  return(p_val_map)
}
