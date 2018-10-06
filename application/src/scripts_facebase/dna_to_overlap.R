#' @param file_list a list of file names for dc gene results.
meta_analysis_of_dna_results <- function(file_list, top = NULL) {
  #########################################################
  #
  # Helper functions
  #
  #########################################################
  
  #' Inport dna output as data.frame
  #' @param file_name the full name (including path) for saved dna results.
  #' @value A data.frame containing dna results.
  get_df_from_dc_results <- function(file_name, top = NULL) {
    df_string <- readChar(file_name, file.info(file_name)$size)
    df_edited <- gsub("(& )|( \\\\\\\\)", "", df_string)
    df_edited <- gsub("( \t)", "\t", df_edited)
    
    df <- read.table(text = df_edited, sep = "\t", header = TRUE, 
                     stringsAsFactors = FALSE)
    
    if(!is.null(top)) df <- df[1:min(nrow(df), top), ]
    return(df)
  }
  
  #' Find overlapping genes accross data.frames.
  #' @param df_list a list of data.frames to compare. The gene names are assumed
  #' to be in the first column.
  #' @value A vector containing the names of genes in the overlap. If the
  #' intersection is empty, NULL is returned.
  get_overlapping_genes <- function(df_list) {
    gene_names <- df_list[[1]][, 1]
    overlap <- sapply(gene_names, function(gene) {
      ifelse(all(sapply(df_list, function(df) gene %in% df[, 1])), TRUE, FALSE)
      })
    if(sum(overlap) == 0) return(NULL)
    return(gene_names[overlap])
  }
  
  
  df_list <- lapply(file_list, get_df_from_dc_results, top = top)
  overlap <- get_overlapping_genes(df_list)
  if(is.null(overlap)) return(NULL)
  
  # Subset the rows of each data.frame in the list
  df_list <- lapply(df_list, function(df) df[df[, 1] %in% overlap, ])
  
  # Use the first data.frame to use as a summary.
  df <- df_list[[1]]
  # Update d_genes as the average of d_genes across data.frames. 
  # d_genes is assumed to be in the second column.
  if(length(overlap) == 1) {
    df[, 2] <- mean(sapply(df_list, function(df) df[, 2]))
  } else {
    df[, 2] <- apply(sapply(df_list, function(df) df[, 2]), 1, mean)
  }
  
  df <- df[order(df[, 2], decreasing = TRUE), ]
  return(df)
}