# Apr. 25, 2017
#
# The .bam files located in /bam_files are converted to count data 
# and stored in a .csv file.
# library(dplyr)
# library(GenomicFeatures) # For makeTxDbFromUCSC() and exonsBy().
# library(biomaRt) # For useMart() and getBM().
# library(Rsamtools) # For BamFileList() and seqinfo().
# library(GenomicAlignments) # For summarizeOverlaps().

#' Converts .bam files to count data
#' 
#' @param dir_input Directory containing folders, which themselves contain 
#' the .bam files.
#' @param dir_output_counts Directory where .csv count files are saved to.
#' @param dir_output_lengths Directory where .csv file containing gene lengths 
#' are saved to.
bam_to_counts <- function(dir_input = "./data/facebase/raw/bam/",
                          dir_output_counts = "./data/facebase/counts/",
                          dir_output_lengths = "./data/facebase/gene_length/") {
  cat("Running script: `script_bam_to_counts` with params:\n")
  cat("\t dir_input =", dir_input, "\n")
  cat("\t dir_output_counts =", dir_output_counts, "\n")
  cat("\t dir_output_lengths =", dir_output_lengths, "\n")
  
  if(!dir.exists(dir_output_counts)) {
    dir.create(dir_output_counts, recursive = TRUE)
  }
  if(!dir.exists(dir_output_lengths)) {
    dir.create(dir_output_lengths, recursive = TRUE)
  }
  
  # Obtain reference genome (txdb) through UCSC and GRangesList of exons. 
  txdb <- makeTxDbFromUCSC(genome = "mm9")
  exons <- exonsBy(txdb, by = "gene") # Type of gene ID: Entrez Gene ID
  
  # Obtain length of each gene and save as .csv file.
  gene_lengths <- sapply(exons, function(exon) sum(exon@ranges@width))
  df <- tibble(entrezgene = as.integer(names(gene_lengths)), 
               length = gene_lengths)
  save_file <- paste0(dir_output_lengths, "gene_length.csv")
  cat("Saving gene lengths to", save_file, "\n")
  write.csv(df, save_file, row.names = FALSE)
  rm(df, gene_lengths) # Delete gene length information; no longer needed.
  gc()
  
  cat("Obtaining MGI gene symbols and other info from biomart.")
  mart <- useMart(biomart = "ensembl", 
                  dataset = "mmusculus_gene_ensembl")
  gene_info <- getBM(attributes = c("entrezgene", "mgi_symbol", "chromosome_name"), 
                     mart = mart)
  # Remove any rows without an entrezgene ID. 
  gene_info <- dplyr::filter(gene_info, !is.na(entrezgene))
  # Remove any rows with mgi_symbol missing.
  gene_info <- dplyr::filter(gene_info, !(is.na(mgi_symbol) | mgi_symbol == "")) # 21585 rows.
  cat("\t", nrow(gene_info), "genes from biomart after removing rows without an mgi symbol.\n")
  
  # Only keep genes in chromosomes 1 - 19 or X. 
  gene_info <- dplyr::filter(gene_info, chromosome_name %in% c(paste(1:19), "X")) # 21094 rows.
  cat("\t", nrow(gene_info), "genes after removing those not on chromosomes 1 - 19 or X.\n")
  
  gene_info <- gene_info %>%
    dplyr::group_by(entrezgene) %>%
    dplyr::summarise(mgi_symbol = first(mgi_symbol))
  cat("\t", nrow(gene_info), "genes after handling duplicate mgi symbols.\n")
  
  # .bam files are stored in different folders. Obtain a list of folder names.
  folders <- list.files(path = dir_input, recursive = FALSE)
  
  # Create a .csv file of counts for each set of .bam files.
  make_csv_of_counts <- function(bam_folder_name, dir_input, dir_output) {
    current_path <- paste(dir_input, bam_folder_name, sep = "")
    cat("Loading .bam files in", current_path, "\n")
    filenames <- file.path(current_path, 
                           list.files(current_path, pattern = "*.bam"))
    
    # Using Rsamtools, create reference to .bam files.
    bamfiles <- BamFileList(filenames)
    seqinfo(bamfiles) # Print out some information about the sequence.
    
    # Align reads to reference genome and obtain counts.
    se <- summarizeOverlaps(features = exons, 
                            reads = bamfiles,
                            mode = "Union")
    counts <- as.data.frame(assay(se))
    sample_names <- paste("sample", 1:length(filenames), sep = "_")
    colnames(counts) <- sample_names
    
    # Join counts to filtered genes obtained from biomart.
    counts <- data.frame(entrezgene = as.integer(rownames(counts)), 
                         counts, stringsAsFactors = FALSE)
    counts <- dplyr::inner_join(gene_info, counts, by = "entrezgene")
    
    # Remove mgi symbol column.
    counts <- dplyr::select(counts, -mgi_symbol)
    
    # Export .csv file containing counts for each entrez gene.
    save_file <- paste(dir_output, bam_folder_name, ".csv", sep = "")
    cat("Saving .csv file to", save_file, "\n")
    write.csv(counts, 
              file = save_file, 
              row.names = FALSE)
    return(save_file)
  }
  
  lapply(folders, make_csv_of_counts, dir_input, dir_output_counts)
}