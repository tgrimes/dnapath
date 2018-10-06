# Tyler Grimes
# December 5, 2017
#
# Differential network analysis
source("src/init_required_packages_and_files.R")

# Create network inference function for dna from association function in SeqNet.
fn_network_inference <- function(fn = pcor_shrinkC) {
  function(x) {
    # We will only work with genes that have nonzero variation. 
    index <- which(apply(x, 2, sd) > 0)
    if(length(index) < 2) return(NA)
    
    # Center and scale the gene expression profile by column.
    x[, index] <- scale(x[, index])
    
    # Estimate the association network, then standardize the matrix of scores.
    scores <- fn(x[, index])
    diag(scores) <- 0
    
    # Complete the association matrix A.
    p <- ncol(x)
    if(length(index) == p) {
      A <- scores
    } else {
      A <- matrix(0, p, p)
      A[index, index] <- scores
    }
    colnames(A) <- colnames(x)
    return(A)
  }
}

# Use .bam files to compute differential expression scores.
dir_output <- "./output/facebase/" # Starting directory for output. A subdirectory will be
# when the data are filtered. Each DNA method will be
# given a subdirectory in the filtered data directory.
dir_bam <- "./data/facebase/raw/bam/"   # Directory containing the folders of .bam files.
dir_counts <- "./data/facebase/counts/" # Directory where raw counts will be stored.
dir_lengths <- "./data/facebase/gene_length/" # Directory where raw counts will be stored.
dir_normalized <- "./data/facebase/normalized/" # Directory where normalized counts will be stored.
dir_log <- "./logs/facebase/"

source("src/scripts_facebase/bam_to_counts.R")
source("src/scripts_facebase/counts_to_normalized.R")
source("src/scripts_facebase/normalized_to_dna.R")
source("src/scripts_facebase/dna_to_summary.R")


# First convert .bam files to counts. 
# Annotation is done using mm9 genome, with genes identified by Entrez Gene ID.
#
# open_log("1_bam_to_counts", dir_log)
# bam_to_counts(dir_input = dir_bam,
#               dir_output_counts = dir_counts,
#               dir_output_lengths = dir_lengths)
# close_log()

open_log("2_counts_to_normalized", dir_log)
counts_to_normalized(dir_input_counts = dir_counts,
                     dir_input_lengths = dir_lengths,
                     dir_output_counts = dir_normalized)
close_log()

#
# Filter raw counts based on low mean expression, low variance, etc.
#
open_log("3_filter_counts", dir_log)
counts_files_list <- paste0(dir_counts, list.files(dir_counts))
#Load and filter counts.
counts_list <- lapply(counts_files_list, function(csv_file) {
  create_counts(csv_file) %>%
    apply_filter_percent_zero_at_most(0.34) %>%
    (function(counts) {
      cat("\t\t", nrow(counts), "genes remaining\n")
      return(counts)
    })
})
names(counts_list) <- gsub(".csv", "", list.files(dir_counts))
filtered_counts_name <- save_counts(counts_list,
                                    file_name = NULL,
                                    dir_output = dir_output,
                                    make_subdir = TRUE)
close_log()

# The filtered data subdirectory will be the root for further output.
dir_filtered <- paste0(dir_output, filtered_counts_name, "/")



#
# 4. Perform dna on normalized counts.
#
filtered_counts_name <- "zeroes0.34_normalized"
dir_filtered <- paste0(dir_output, filtered_counts_name, "/")
lp_set <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
signif <- 0.1

open_log("4_dna_corr", dir_log)
normalized_to_dna(dir_filtered,
                  network_inference = fn_network_inference(corC),
                  lp_set = lp_set,
                  dir_label = "corr")
close_log()
gc()

open_log("4_dna_pcor", dir_log)
normalized_to_dna(dir_filtered,
                  network_inference = fn_network_inference(pcor_shrinkC),
                  lp_set = lp_set,
                  dir_label = "pcor")
close_log()
gc()



#
# 5. Summarize dna results.
# 
filtered_counts_name <- "zeroes0.34_normalized"
dir_filtered <- paste0(dir_output, filtered_counts_name, "/")
signif <- 0.1

open_log("5_summary_corr", dir_log)
lapply(c("Lateral", "Medial", "Nasal", "Oral"), function(pair_name) {
  summarize_results_from_dna(dataset = "facebase",
                             measure = "corr",
                             lp = 2,
                             pair_name = pair_name,
                             dir_output = dir_filtered,
                             signif = signif)
})
close_log()
gc()

open_log("5_summary_pcor", dir_log)
lapply(c("Lateral", "Medial", "Nasal", "Oral"), function(pair_name) {
  summarize_results_from_dna(dataset = "facebase",
                             measure = "pcor",
                             lp = 2,
                             pair_name = pair_name,
                             dir_output = dir_filtered,
                             signif = signif)
})
close_log()
gc()
