#' Differential network analysis
#' 
#' @param counts_pair a list of two p by n data.frames of gene expression counts. 
#' The first column of each data.frame should be named entrezgene and contain
#' the entrezgene id's for each gene.
#' @param n_perm the number of random permutations to perform during 
#' permutation testing. If n_perm == 1, the permutation tests are not performed.
#' @param pathway_name optional name for labeling the pathway associated with
#' the given counts.
#' @param network_inference a function used to infer the pathway network. Should
#' return a p by p matrix of association scores.
#' @param diff_connectivity a function used to compute differential connectivity
#' scores of each gene. Should return a vector of length p.
#' @return a list containing results for each pathway. The p-values for 
#' differential connectivity of each gene are given, along with the overall
#' differential connectivity of the pathway. Pathways without any significance
#' are excluded from the results.
dna <- function(counts_joined, pathway_list, permutations, network_inference, 
                lp_set = 1) {
  if(is.null(colnames(counts_joined)))
    colnames(counts_joined) <- 1:ncol(counts_joined)
  
  if(is.null(pathway_list) || is.na(pathway_list))
    pathway_list = list(no_pathway_info = colnames(counts_joined))
  
  if(is.list(pathway_list)) {
    pathway_list <- lapply(pathway_list, sort)
    results_by_pathway <- lapply(pathway_list, function(pathway) {
      dna(counts_joined[, pathway], pathway, permutations, network_inference, lp_set)
    })
    
    # Re-arrange results by lp, then by pathway.
    results_by_lp <- vector("list", length(lp_set))
    for(i in 1:length(lp_set)) {
      results_by_lp[[i]] <- lapply(results_by_pathway, function(results) results[[i]])
    }
    
    return(results_by_lp)
  }
  
  
  n <- nrow(counts_joined) / 2 # Number of observations for each network.
  p <- ncol(counts_joined) # Total number of genes in union.
  n_perm <- ncol(permutations)
  
  # Store differential connectivity scores from each permutation.
  scores_1 <- network_inference(counts_joined[1:n, ])
  scores_2 <- network_inference(counts_joined[-(1:n), ])
  
  d_path_list <- lapply(lp_set, function(lp) d_pathwayC(scores_1, scores_2, lp = lp))
  d_gene_list <- lapply(lp_set, function(lp) d_genesC(scores_1, scores_2, lp = lp))
  d_edge_list <- lapply(lp_set, function(lp) d_edgesC(scores_1, scores_2, lp = lp))
  
  # Obtain ranks of original DC scores for monotonization.
  d_path_ranks <- lapply(d_path_list, function(d_path) order(d_path, decreasing = TRUE))
  d_gene_ranks <- lapply(d_gene_list, function(d_gene) order(d_gene, decreasing = TRUE))
  d_edge_ranks <- lapply(d_edge_list, function(d_edge) order(d_edge, decreasing = TRUE))
  
  pval_path_list <- lapply(lp_set, function(lp) 1)
  pval_gene_list <- lapply(lp_set, function(lp) rep(1, p))
  pval_edge_list <- lapply(lp_set, function(lp) rep(1, p * (p - 1) / 2))
  
  pval_path_list_mono <- lapply(lp_set, function(lp) 1)
  pval_gene_list_mono <- lapply(lp_set, function(lp) rep(1, p))
  pval_edge_list_mono <- lapply(lp_set, function(lp) rep(1, p * (p - 1) / 2))
  
  if(!is.null(n_perm) && n_perm > 1) {
    for(i in 1:n_perm) {
      network_1 <- permutations[, i]
      
      scores_1 <- network_inference(counts_joined[network_1, ])
      scores_2 <- network_inference(counts_joined[-network_1, ])
      
      for(k in 1:length(lp_set)) {
        pval_path_list[[k]] <- pval_path_list[[k]] +
          (d_pathwayC(scores_1, scores_2, lp = lp_set[[k]]) >= d_path_list[[k]])
        pval_gene_list[[k]] <- pval_gene_list[[k]] +
          (d_genesC(scores_1, scores_2, lp = lp_set[[k]]) >= d_gene_list[[k]])
        pval_edge_list[[k]] <- pval_edge_list[[k]] +
          (d_edgesC(scores_1, scores_2, lp = lp_set[[k]]) >= d_edge_list[[k]])
        
        # Compute step-down p-values.
        # Rank index for original scores:
        index <- d_path_ranks[[k]]
        n_index <- length(index)
        # Calculate scores of permuted data.
        d_path_perm <- d_pathwayC(scores_1, scores_2, lp = lp_set[[k]])
        # Take the cumulative max, in rank order, starting from the min.
        d_path_perm <- cummax(d_path_perm[rev(index)])[n_index:1][order(index)]
        # Accumulate number of permuted scores that are larger than original.
        pval_path_list_mono[[k]] <- pval_path_list_mono[[k]] +
          (d_path_perm >= d_path_list[[k]])
        
        # Repeat step-down procedure for genes and edges.
        index <- d_gene_ranks[[k]]
        n_index <- length(index)
        d_gene_perm <- d_genesC(scores_1, scores_2, lp = lp_set[[k]])
        d_gene_perm <- cummax(d_gene_perm[rev(index)])[n_index:1][order(index)]
        pval_gene_list_mono[[k]] <- pval_gene_list_mono[[k]] +
          (d_gene_perm >= d_gene_list[[k]])
        
        index <- d_edge_ranks[[k]]
        n_index <- length(index)
        d_edge_perm <- d_edgesC(scores_1, scores_2, lp = lp_set[[k]])[index]
        d_edge_perm <- cummax(d_edge_perm[rev(index)])[n_index:1][order(index)]
        pval_edge_list_mono[[k]] <- pval_edge_list_mono[[k]] +
          (d_edge_perm >= d_edge_list[[k]])
      }
    }
    
    pval_path_list <- lapply(pval_path_list, function(pval) pval / (n_perm + 1))
    pval_gene_list <- lapply(pval_gene_list, function(pval) pval / (n_perm + 1))
    pval_edge_list <- lapply(pval_edge_list, function(pval) pval / (n_perm + 1))
    
    for(k in 1:length(lp_set)) {
      index <- d_path_ranks[[k]]
      pval_path_list_mono[[k]] <- 
        cummax(pval_path_list_mono[[k]][index] / (n_perm + 1))[order(index)]
      index <- d_gene_ranks[[k]]
      pval_gene_list_mono[[k]] <- 
        cummax(pval_gene_list_mono[[k]][index] / (n_perm + 1))[order(index)]
      index <- d_edge_ranks[[k]]
      pval_edge_list_mono[[k]] <- 
        cummax(pval_edge_list_mono[[k]][index] / (n_perm + 1))[order(index)]
    }
  }
  
  results <- vector("list", length(lp_set))
  for(i in 1:length(lp_set)) {
    results[[i]] <- list(gene_names = as.numeric(colnames(counts_joined)),
                         p_value_path = pval_path_list[[i]],
                         p_value_genes = pval_gene_list[[i]],
                         p_value_edges = pval_edge_list[[i]],
                         p_value_path_mono = pval_path_list_mono[[i]],
                         p_value_genes_mono = pval_gene_list_mono[[i]],
                         p_value_edges_mono = pval_edge_list_mono[[i]],
                         d_path = d_path_list[[i]],
                         d_genes = d_gene_list[[i]],
                         d_edges = d_edge_list[[i]],
                         lp = lp_set[[i]])
  }
  
  return(results)
}

#' Summarize differential connectivity scores over different pathways
#' 
#' @param results_list a list of results from performing [dna()] on several pathways.
#' @param gene_list a set of entrezgene IDs to be summarized.
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used. 
#' @return The p by p differential network, with nonzero entries containing
#' the significant DC scores.
summarize_dnw_over_pathways <- function(results_list, gene_list = NULL, 
                                        signif = 0.05, monotonized = FALSE) {
  if(is.null(gene_list)) {
    gene_list <- unique(unlist(lapply(results_list, 
                                      function(results) results$gene_names)))
  }
  
  p <- length(gene_list)
  dnw <- matrix(0, p, p)
  colnames(dnw) <- gene_list
  
  for(results in results_list) {
    # Compute the differential network for given pathway.
    dnw_p <- summarize_dnw_for_pathway(results, 
                                       signif = signif, 
                                       monotonized = monotonized)
    
    # Index the columns so they can be aligned with the overall network.
    index <- sapply(results$gene_names, 
                    function(gene) which(gene_list == gene)[1])
    
    # Keep the largest dc score for each gene.
    dnw[index, index] <- ifelse(dnw_p >= dnw[index, index], 
                                dnw_p, 
                                dnw[index, index])
  }
  
  return(dnw)
}

#' Summarize differential connectivity scores for a pathway
#' 
#' @param results The results from performing [dna()] on a single pathway.
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used. 
#' @return The p by p differential network, with nonzero entries containing
#' the significant DC scores.
summarize_dnw_for_pathway <- function(results, signif = 0.05, monotonized = FALSE) {
  genes_in_dnw <- results$gene_names
  p <- length(genes_in_dnw)
  
  if(p == 0) 
    stop("results pathway does not contain any genes.")
  
  dnw <- matrix(0, p, p)
  colnames(dnw) <- genes_in_dnw
  
  d_edges <- results$d_edges
  
  # If p-value is above signif, set to zero.
  if(monotonized) {
    d_edges[results$p_value_edges_mono > signif] <- 0
  } else {
    d_edges[results$p_value_edges > signif] <- 0
  }
  
  dnw[lower.tri(dnw)] <- d_edges
  dnw <- dnw + t(dnw)
  
  return(dnw)
}


#' Summarize differential connectivity of genes over different pathways
#' 
#' @param results_list a list of results from performing [dna()] on several pathways.
#' @param gene_list a set of entrezgene IDs to be summarized.
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used. 
#' @return A vector of DC gene scores, with nonzero entries containing
#' the significant DC scores.
summarize_genes_over_pathways <- function(results_list, gene_list = NULL, 
                                          signif = 0.05, monotonized = FALSE) {
  if(is.null(gene_list)) {
    gene_list <- unique(unlist(lapply(results_list, 
                                      function(results) results$gene_names)))
  }
  
  p <- length(gene_list)
  d_genes <- rep(0, p)
  names(d_genes) <- gene_list
  
  for(results in results_list) {
    d_genes_p <- summarize_genes_for_pathway(results,
                                             signif = signif,
                                             monotonized = monotonized)
    index <- sapply(results$gene_names, 
                    function(gene) which(gene_list == gene)[1])
    
    d_genes[index] <- ifelse(d_genes_p >= d_genes[index], 
                             d_genes_p, 
                             d_genes[index])
  }
  
  return(d_genes)
}

#' Summarize differential connectivity scores for a pathway
#' 
#' @param results The results from performing [dna()] on a single pathway.
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used.
#' @return The p by p differential network, with nonzero entries containing
#' the significant DC scores.
summarize_genes_for_pathway <- function(results, signif = 0.05, monotonized = FALSE) {
  gene_names <- unique(results$gene_names)
  
  d_genes <- drop(results$d_genes)
  names(d_genes) <- gene_names
  
  if(monotonized) {
    d_genes[drop(results$p_value_genes_mono) > signif] <- 0
  } else {
    d_genes[drop(results$p_value_genes) > signif] <- 0
  }
  
  return(d_genes)
}

#' Summarize differential connectivity scores from different pathways
#' 
#' @param results_list a list of results from performing dna() on several pathways.
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used.
#' @return list of results
summarize_pathways <- function(results_list, signif = 0.05, monotonized = FALSE) {
  pathway_names <- names(results_list)
  
  if(monotonized) {
    p_values <- sapply(results_list, function(results) results$p_value_path_mono)
  } else {
    p_values <- sapply(results_list, function(results) results$p_value_path)
  }
  
  d_pathway <- sapply(results_list, function(results) results$d_path)
  
  d_pathway[p_values > signif] <- 0
  
  return(d_pathway)
}
