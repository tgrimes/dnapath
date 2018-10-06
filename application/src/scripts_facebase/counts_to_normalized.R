# Tyler Grimes
# January 13, 2017
#
# Differential network analysis

#' Compute association scores from counts
#' 
#' @param dir_input Directory containing the count data as .csv files.
#' @param dir_output Directory where normalized .csv files will be saved.
counts_to_normalized <- function(dir_input_counts = "./data/facebase/counts/",
                                 dir_input_lengths = "./data/facebase/gene_length/",
                                 dir_output_counts = "./data/facebase/normalized/") {
  cat("Running: `counts_to_normalized()` with params:\n")
  cat("\t dir_input_counts =", dir_input_counts, "\n")
  cat("\t dir_input_lengths =", dir_input_lengths, "\n")
  cat("\t dir_output_counts =", dir_output_counts, "\n")
  
  if(!dir.exists(dir_output_counts)) {
    dir.create(dir_output_counts, recursive = TRUE)
  }
  
  load_file <- paste0(dir_input_lengths, list.files(dir_input_lengths)[[1]])
  cat("Loading gene lengths from", load_file, "\n")
  gene_length <- read.csv(load_file)
  
  # Obtain a list of files in input directory.
  file_names <- list.files(path = dir_input_counts, recursive = FALSE)
  file_names <- file_names[grepl(".csv", file_names)] # Subset on .csv files.
  
  # Load .csv file, normalize the counts, then save to output directory.
  normalize_counts <- function(csv_file_name) {
    file_name <- paste(dir_input_counts, csv_file_name, sep = "")
    cat("\tLoading .csv file from", file_name, "\n")
    counts <- read.csv(file_name)
    
    cat("\t - applying TPM normalization and log2(1 + x) transformation.\n")
    # Use TPM normalization on counts.
    counts <- tpm(counts, gene_length)
    # Apply log2(1 + x) transform to counts.
    counts[, -1] <- log(1 + counts[, -1], 2)
    
    
    # Export .csv file containing counts for each entrez gene.
    save_file <- paste(dir_output_counts, csv_file_name, sep = "")
    cat("\tSaving normalized .csv file to", save_file, "\n")
    write.csv(counts, 
              file = save_file, 
              row.names = FALSE)
    return(save_file)
  }
  
  lapply(file_names, normalize_counts)
}

#' Transcripts per kilobase million (TPM) normalization
#' 
#' @param counts a p by (n + 1) matrix of read counts. The first column is 
#' assumed to contain identifiers for each gene; the name for this column 
#' should be specified by the parameter 'id'. The remaining columns give the 
#' counts for each sample.
#' @param gene_length a p by 2 matrix containing gene lengths. The first column  
#' is assumed to contain identifiers for each gene; the name for this column 
#' should be specified by the parameter 'id'. The second column should give the
#' length of each gene.
#' @param id the name of the column containing gene identifiers. 
tpm <- function(counts, gene_length, id = "entrezgene") {
  df <- right_join(gene_length, counts, by = id)
  # Remove identifier column. First column of df will contain lengths.
  df <- dplyr::select(df, -contains(id)) 
  
  # 1. Divide read counts by gene length (in kilobases).
  normalized <- t(apply(df, 1, function(x) x[-1] / x[1] * 1000))
  # 2. Divide RPK values by total number of RPK in sample per million. 
  normalized <- apply(normalized, 2, function(x) x / sum(x) * 10^6)
  
  normalized <- cbind(dplyr::select(counts, id), normalized)
  return(normalized)
}
