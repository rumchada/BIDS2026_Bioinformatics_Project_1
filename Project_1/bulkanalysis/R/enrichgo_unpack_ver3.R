enrichgo_unpack_ver3 <- function(
    # Pass the whole list so we can grab comparisons
  initial_table_list = filtered_results,
  key = ora_output_down$enrichgo_results,
  pval_adj_threshold = 0.001,

  # Assume user does not want to convert geneids
  convert = FALSE,
  from_type = "ensembl",
  to_type = "symbol",
  ensembl_dataset = "hsapiens_gene_ensembl"){
  require(tibble)
  require(dplyr)

  ### helper functionf or converting geneidA to geneidB
  ### user
  gene_id_converter <- function(vector, from_type, to_type, ensembl_dataset){
    require(biomaRt)
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


  # FIX 2: Corrected the spelling typo here
  unpacking_results <- list()

  # FIX 1: Loop over names(key) so 'comparison' is a valid character string identifier
  for (comparison in names(key)) {

    # Safeguard check: If the comparison doesn't exist in the cluster results, skip
    if (is.null(key[[comparison]])) next

    cluster <- tibble::as_tibble(key[[comparison]]@result) %>%
      dplyr::mutate(Description = as.character(Description))

    if (!is.null(pval_adj_threshold)) {
      significant_clusters <- cluster %>%
        dplyr::filter(p.adjust < pval_adj_threshold)
    } else {
      significant_clusters <- cluster
    }

    labeled_gene_lists <- list()

    # FIX 3: Pull the specific comparison dataframe from your filtered results list
    current_initial_table <- initial_table_list[[comparison]] %>%
      dplyr::mutate(geneid = as.character(geneid))

    # If there are no significant clusters, jump to saving the empty list or next comparison
    if (nrow(significant_clusters) == 0) {
      unpacking_results[[comparison]] <- labeled_gene_lists
      next
    }

    for (i in seq_len(nrow(significant_clusters))) {

      go_term <- significant_clusters$Description[i]

      # Extract genes for current GO term
      genes <- unlist_converter(significant_clusters$geneID[i])

      pathway_tibble <- tibble::tibble(
        gene_desc = rep(go_term, length(genes)),
        geneid = as.character(genes)
      )

      if (convert == FALSE) {
        pathway_genenames <- pathway_tibble
      } else {
        pathway_genenames <- gene_id_converter(
          pathway_tibble$geneid,
          from_type = from_type,
          to_type = to_type,
          ensembl_dataset = ensembl_dataset
        )

        new_colname <- paste0(to_type, "_geneid")

        # FIX 4: Cleaner, bulletproof column renaming by position
        colnames(pathway_genenames)[1] <- "geneid"
        colnames(pathway_genenames)[2] <- new_colname
      }

      # Join converted names
      final_path_tibble <- dplyr::inner_join(pathway_tibble, pathway_genenames, by = "geneid") %>%
        dplyr::mutate(geneid = as.character(geneid))

      # Join fold changes / statistics using our dynamically updated initial table
      final_path_fc_tibble <- dplyr::inner_join(final_path_tibble, current_initial_table, by = "geneid")

      # Relocate converted gene column
      if (convert == TRUE) {
        fpfc_tibble <- final_path_fc_tibble %>%
          dplyr::relocate(dplyr::all_of(new_colname), .after = geneid)
      } else {
        fpfc_tibble <- final_path_fc_tibble
      }

      # Store result under the GO Term name
      labeled_gene_lists[[go_term]] <- fpfc_tibble
    }

    # Store the list of terms into the master list under the Comparison name
    unpacking_results[[comparison]] <- labeled_gene_lists
  }
  return(unpacking_results)
}
