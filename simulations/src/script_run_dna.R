if(!exists("results_file")) {
  stop("Variable `results_file` does not exist in global environment.")
}
if(!exists("seed_network")) {
  stop("Variable `seed_network` needs to be specified.")
}
if(!exists("network_number")) {
  stop("Variable `network_number` needs to be specified.")
}
if(!exists("output_dir")) {
  stop("Variable `output_dir` needs to be specified.")
}
if(!exists("n_iter_set")) {
  stop("Variable `n_iter_set` does not exist in global environment.")
}
if(!exists("n_set")) {
  stop("Variable `n_set` does not exist in global environment.")
}
if(!exists("lp_set")) {
  stop("Variable `lp_set` does not exist in global environment.")
}
if(!exists("measure_set")) {
  stop("Variable `measure_set` does not exist in global environment.")
}
if(!exists("signal_set")) {
  stop("Variable `signal_set` does not exist in global environment.")
}
if(!exists("n_perm")) {
  stop("Variable `n_perm` does not exist in global environment.")
}
if(!exists("signif_set")) {
  stop("Variable `signif_set` does not exist in global environment.")
}

#
# Load the two networks
#
network <- get_networks()
pathway_list <- get_pathway_list(network[[1]])

dnw <- abs(get_adj_matrix_from_network(network[[1]]) -
             get_adj_matrix_from_network(network[[2]]))
truth_edges <- dnw[lower.tri(dnw)]
truth_genes <- 1 * (apply(dnw, 2, sum) > 0)
truth_paths <- sapply(pathway_list, function(path) 1 * (sum(dnw[path, path]) > 0))

#
# Load the generated x data
#
gen <- get_generated_data()

counts_joined <- rbind(gen[[1]]$x, gen[[2]]$x)
colnames(counts_joined) <- colnames(gen[[1]]$x) 
p <- ncol(counts_joined)

rm(gen, dnw, network)
gc()


add_noise_to_pathway <- function(nodes, p, signal) {
  if(signal == 1) return(nodes)
  if(p <= length(nodes)) warning("p expected to be greater than number of nodes.")
  # Replace (1 - signal)*100% of nodes with random nodes in {1:p}\nodes.
  pert_path <- sample(nodes, round(signal * length(nodes)))
  remaining_nodes <- 1:p
  remaining_nodes <- remaining_nodes[!(remaining_nodes %in% nodes)]
  pert_path <- c(pert_path, sample(remaining_nodes, 
                                   round((1 - signal) * length(nodes))))
  return(pert_path)
}


cat("... saving results to", results_file, "\n")
# If results file does not exist, create it and add the header.
if(!file.exists(results_file)) {
  file.create(results_file)
  sink(results_file, append = TRUE)
  header <- c("time", "sample_number", "network_number", "network_seed", "component", 
              "n", "p", "lp", "signal", "perm", "signif", "measure", "monotonized",
              "TP", "FN", "FP", "TN")
  cat(header, sep = "\t")
  cat("\n")
} else {
  sink(results_file, append = TRUE)
}


N <- max(n_set)
M <- iter_max * N
# Begin simulation over simulation parameters.
tryCatch({
  for(n in n_set) {
    for(i in n_iter_set) {
      
      df <- counts_joined[c((N * (i - 1) + 1):(N * (i - 1) + n),
                            (N * (i - 1) + 1):(N * (i - 1) + n) + M), ]
      
      # Seed the permutation samples (function of sample number and sample size).
      set.seed(i + 1000 * n)
      permutations <- sapply(1:n_perm, function(k) sample(1:(2 * n), n))
      
      for(j in 1:length(measure_set)) {
        for(signal in signal_set) {
          
          if(is.na(signal)) {
            pathways <- list(no_pathway_info = 1:p)
          } else {
            pathways <- pathway_list
            pathways <- lapply(pathways, add_noise_to_pathway, 
                               p = p, signal = signal)
          }
          
          # Perform the differential network analysis.
          dna_results_by_lp <- dna(df, pathways, permutations = permutations,
                                   network_inference = measure_set[[j]], 
                                   lp_set = lp_set)
          
          for(k in 1:length(lp_set)) {
            dna_results <- dna_results_by_lp[[k]]
            
            # For each alpha level, record performance on genes, edges, and pathways.
            for(signif in signif_set[signif_set >= (1 / (n_perm + 1))]) {
              for(monotonized in c(FALSE, TRUE)) {
                #
                # Evaluate performance on genes.
                #
                dc_genes <- summarize_genes_over_pathways(dna_results, 
                                                          gene_list = 1:p, 
                                                          signif = signif,
                                                          monotonized = monotonized)
                perf <- get_confusion_matrix(truth_genes, dc_genes)
                
                cat(gsub("(.* )", "", Sys.time()), "\t")
                results <- c(i, network_number, seed_network, "genes", n, p, 
                             lp_set[k], signal, n_perm, signif, names(measure_set)[j], 
                             monotonized, perf$TP, perf$FN, perf$FP, perf$TN)
                cat(results, sep = "\t")
                cat("\n")
                
                #
                # Evaluate performance on edges.
                #
                dc_edges <- summarize_dnw_over_pathways(dna_results, 
                                                        gene_list = 1:p, 
                                                        signif = signif,
                                                        monotonized = monotonized)
                dc_edges <- dc_edges[lower.tri(dc_edges)]
                perf <- get_confusion_matrix(truth_edges, dc_edges)
                
                cat(gsub("(.* )", "", Sys.time()), "\t")
                results <- c(i, network_number, seed_network, "edges", n, p, 
                             lp_set[k], signal, n_perm, signif, names(measure_set)[j], 
                             monotonized, perf$TP, perf$FN, perf$FP, perf$TN)
                cat(results, sep = "\t")
                cat("\n")
              } # End for monotonized in T/F
              
              #
              # Evaluate performance on pathways (if pathway info provided).
              #
              # Note: monotonized (step-down) p-values do not apply here.
              if(!is.na(signal)) {
                dc_paths <- summarize_pathways(dna_results, 
                                               signif = signif,
                                               monotonized = FALSE)
                perf <- get_confusion_matrix(truth_paths, dc_paths)
                
                cat(gsub("(.* )", "", Sys.time()), "\t")
                results <- c(i, network_number, seed_network, "paths", n, p, 
                             lp_set[k], signal, n_perm, signif, names(measure_set)[j], 
                             FALSE, perf$TP, perf$FN, perf$FP, perf$TN)
                cat(results, sep = "\t")
                cat("\n")
              }
            } #End for signal in signal_set
          } # End for lp in lp_set
        } # End for signal in signal_set
      } # End for measure in measure_set
    } # End for i in n_iter_set
  } # End for n in n_set
  
  # Close the connection to the results file.
  sink()
  
}, error = function(e) {
  # Close the connection to the results file and print error.
  sink()
  print(e)
})
