normalized_to_dna <- function(dir_filtered,
                              network_inference, 
                              lp_set = c(1, 2), 
                              dir_label = "New results") {
  if(is.null(dir_label) | dir_label == "") {
    stop("dir_label cannot be null nor an empty string.")
  }
  
  # Load in counts list.
  counts_pair <- readRDS(paste0(dir_filtered, filtered_counts_name))
  
  # Get list of reactome pathways. Subset on those of interest.
  # pathway_list <- get_reactome_pathways("mouse")
  # pathway_list <- pathway_list[sapply(pathway_list, length) %in% 10:100]
  # pathway_list <- remove_redundant_pathways(pathway_list, threshold = 0.9)
  # saveRDS(pathway_list, "data/reactome/Mus_musculus_10to100_nonredundant0.9.rds")
  pathway_list <- get_reactome_pathways(file_name = "Homo_sapiens_10to100_nonredundant0.9.rds")
  cat(length(pathway_list), "Reactome pathways found.\n")
  
  # Create directory for output tables (DC genes and pathways)
  # If this main directory exists, then it's assumed the analysis has already
  # been performed; we stop and return an error. 
  dir_results <- paste0(dir_filtered, "dna_results/", dir_label, "/")
  if(dir.exists(dir_results)) stop("Result for this analysis already exist.")
  dir.create(dir_results, recursive = TRUE)
  
  # Save the network_inference function. This is the key difference between
  # different analyses.
  save_file <- paste0(dir_results, "network_inference")
  cat("... saving network_inference function to", save_file, "\n")
  saveRDS(network_inference, file = save_file)
  
  # Perform differential network analysis on each pathway.
  results_by_lp <-  dna(counts_pair, pathway_list, n_perm = 100, 
                        network_inference = network_inference, lp_set = lp_set)
  
  for(k in 1:length(results_by_lp)) {
    results_list <- results_by_lp[[k]]
    
    results_list <- results_list[!sapply(results_list, is.null)]
    results_list <- results_list[!sapply(results_list, function(x) is.na(x$p_value_path))]
    
    save_file <- paste0(dir_results, "lp", lp_set[k], "_neuro.rds")
    cat("Saving dna results to", save_file, "\n")
    saveRDS(results_list, save_file)
  } # End for k in lp_set
}