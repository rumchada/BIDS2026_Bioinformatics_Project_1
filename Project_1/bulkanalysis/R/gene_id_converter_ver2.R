gene_id_converter_ver2 <- function(vector, from_type, to_type, ensembl_datset){
  require(biomaRt)
  # geneid pulling for correct attribute
  id_map <- c(
    ensembl = "ensembl_gene_id",
    entrez  = "entrezgene_id",
    symbol  = "external_gene_name"
  )
  # Pass the the ID mapping index over to the new attribute string
  from_attr <- id_map[[from_type]]
  to_attr <-id_map[[to_type]]

  #initialize the return and parameters
  #baseline NULL intialization
  return_df <- NULL
  max_attempts <- 2
  attempt <- 1
  org_db <- org.Hs.eg.db


  #Fix 1: There is a huge API crash almost 50% of the time that I use getBM() function to connect to the API
  #Implementing a while loop/tryCatch function to loop API calls just in case that the getBM work
  #Solution:After 10 attempts,

  #while the current attemp is less than or equal to the max attempts and while the return DF is null
  while(attempt <= max_attempts && is.null(return_df)){
    return_df <- tryCatch({

      message(sprintf("biomaRt attempt %d of %d...", attempt, max_attempts))

      #call the getBM ensemble API
      ensembl <- useMart(biomart = "ensembl", dataset = ensembl_dataset)
      getBM( attributes = c(from_attr, to_attr), filters = from_attr, values = vector, mart = ensembl)

      #If there is any kind of error
    }, error=function(e){
      #message for that error
      message(sprintf("  -> biomaRt attempt %d failed: %s", attempt, trimws(e$message)))
      #Return a NULL DF
      return(NULL)
    })

    # During that current attempt, If the return_df is still null
    if (is.null(return_df)) {
      # add on the attempt counter
      attempt <- attempt + 1
      # If the attempt counter is still less than the max attempts
      if (attempt <= max_attempts) {
        Sys.sleep(2) # Brief pause to prevent hammering the server
      }
    }# end of second if statement
  }#end of while loop

  if(is.null(return_df)){
    dbi_map <- c(
      ensembl = "ENSEMBL",
      entrez  = "ENTREZID",
      symbol  = "SYMBOL"
    )

    dbi_from <- dbi_map[[from_type]]
    dbi_to   <- dbi_map[[to_type]]

    return_df <- AnnotationDbi::mapIds(org_db,
                                       keys = vector,
                                       keytype = dbi_from,
                                       column =   dbi_to,
                                       multiVals = "first")

    return_df <- data.frame(
      from_attr = names(return_df),
      to_attr = unname(return_df))

    colnames(return_df) <- c(from_attr, to_attr)
  }
  return(return_df)
}#end of whole function
