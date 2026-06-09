ora_dotplot <- function(ora_results, gspval_cutoff, org_by, num_disp, graph_title, file_name){
  if(org_by %in% 'Count'){
    #inside if
    dotplot <- ora_results%>%
      dplyr::arrange(desc(Count))%>%
      utils::head(num_disp)%>%
      dplyr::mutate(as.factor(Description))%>%
      ggplot() +
      aes(y = Count, x = fct_reorder(Description, Count)) +
      scale_colour_gradient(low = "red", high = "blue") +
      geom_point(aes(size = GeneRatio, color = pvalue)) +
      labs(title = graph_title, x = "Enriched Term", y = "Gene Count") +
      coord_flip() +
      theme(text = element_text(face = "bold") , plot.title = element_text(hjust = 1))

  }
  if(org_by %in% 'pvalue'){
    #inside if
    dotplot <- ora_results %>%
      dplyr::filter(pvalue < gspval_cutoff) %>%
      utils::head(num_disp) %>%
      dplyr::mutate(as.factor(Description)) %>%
      dplyr::mutate(pvalue = -1 * log(pvalue)) %>%
      ggplot() +
      aes(x = pvalue, y = fct_reorder(Description, pvalue))+
      scale_colour_gradient(low = "red", high = "blue") +
      geom_point(aes(size = GeneRatio, color = Count)) +
      labs(title = graph_title, x = "-log10(P)", y = "Enriched Term") +
      theme(text = element_text(face = "bold"), plot.title = element_text(hjust = 1))


  }
  return(list(dotplot, ggsave(file_name, device = "png", width = 8, height = 6, units = "in")))
}


ora_enrichgo <- function(filtered_results,
                         direction = c("up", "down"),
                         log2fc_thresh = 1,
                         pval_thresh = 0.05,
                         top_terms = 15,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL",
                         ont = "BP"){

  direction <- match.arg(direction)

  # store enrichment results
  enrichgo_results <- vector("list", length(filtered_results))
  names(enrichgo_results) <- names(filtered_results)

  # run enrichGO
  for(i in seq_along(filtered_results)){

    if(direction == "up"){

      genes <- filtered_results[[i]]$geneid[
        filtered_results[[i]]$log2foldchange > log2fc_thresh
      ]

    } else if(direction == "down"){

      genes <- filtered_results[[i]]$geneid[
        filtered_results[[i]]$log2foldchange < -log2fc_thresh
      ]

    }

    enrichgo_results[[i]] <- clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont,
      pAdjustMethod = "BH"
    )
  }

  # store visualization outputs
  enrichgo_visuals <- vector("list", length(enrichgo_results))
  names(enrichgo_visuals) <- names(enrichgo_results)

  # generate plots
  for(i in seq_along(enrichgo_results)){

    ora_dotplot_count <- ora_dotplot(
      enrichgo_results[[i]],
      pval_thresh,
      "Count",
      top_terms,
      glue::glue(
        "ORA Results Ranked by Count: {names(enrichgo_results)[i]}"
      ),
      glue::glue(
        "ORA Results Ranked by Count: {names(enrichgo_results)[i]}.png"
      )
    )

    ora_dotplot_pvalue <- ora_dotplot(
      enrichgo_results[[i]],
      pval_thresh,
      "pvalue",
      top_terms,
      glue::glue(
        "ORA Results Ranked by PValue: {names(enrichgo_results)[i]}"
      ),
      glue::glue(
        "ORA Results Ranked by PValue: {names(enrichgo_results)[i]}.png"
      )
    )

    enrichgo_visuals[[i]] <- list(
      ora_count_ranked = ora_dotplot_count,
      ora_pvalue_ranked = ora_dotplot_pvalue
    )
  }

  return(list(
    enrichgo_results = enrichgo_results,
    enrichgo_visuals = enrichgo_visuals
  ))
}
