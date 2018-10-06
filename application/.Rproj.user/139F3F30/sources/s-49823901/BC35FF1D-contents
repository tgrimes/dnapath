
#' Differential network analysis
#' 
#' @param counts_pair a list of two p by n data.frames of gene expression counts. 
#' The first column of each data.frame should be named entrezgene and contain
#' the entrezgene id's for each gene.
#' @param pathway a vector containing entrezgene IDs that are in a pathway.
#' (Can also be a list of vectors corresponding to a set of pathways).
#' These are used to subset the data.frames in counts_pair.
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
dna <- function(counts_pair, pathway_list, n_perm = 100, pathway_name = NULL,
                network_inference = function(x) {
                  # We will only work with genes that have nonzero variation. 
                  index <- which(apply(x, 2, sd) > 0)
                  if(length(index) < 2) return(NA)
                  
                  # Center and scale the gene expression profile by column.
                  x[, index] <- apply(x[, index], 2, standardize_scores,
                                      ignore_zeroes = FALSE, robust = FALSE)
                  
                  # Estimate the association network. Center and scale the scores.
                  scores <- run_pcor(x[, index])$scores
                  scores <- standardize_scores(scores, ignore_zeroes = TRUE, robust = FALSE)
                  
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
                }, 
                lp_set = 2, seed = NULL) {
  
  if(!is.list(counts_pair))  stop("counts_pair should be a list.")
  if(length(counts_pair) != 2) stop("counts_pair should be a list of two data.frames.")
  if(colnames(counts_pair[[1]])[1] != "entrezgene") 
    stop("entrezgene is not in the first column of the counts_pair[[1]].")
  if(colnames(counts_pair[[2]])[1] != "entrezgene") 
    stop("entrezgene is not in the first column of the counts_pair[[2]].")
  
  # If a list of pathways are given, perform dna over each pathway.
  if(!is.null(pathway_list) & is.list(pathway_list)) {
    results_by_pathway <- lapply(pathway_list, function(pathway) {
      dna(counts_pair, pathway, n_perm, pathway_name, network_inference,
          lp_set, seed)
    })
    
    # Re-arrange results by lp, then by pathway.
    if(length(lp_set) > 1) {
      results_by_lp <- vector("list", length(lp_set))
      for(i in 1:length(lp_set)) {
        results_by_lp[[i]] <- lapply(results_by_pathway, 
                                     function(results) results[[i]])
      }
    } else {
      results_by_lp <- results_by_pathway
    }
    
    names(results_by_lp) <- lp_set
    return(results_by_lp)
    
  } else {
    pathway <- pathway_list
  }
  
  # Set seed for reproducibility of permutation test.
  if(!is.null(seed)) {
    set.seed(seed)
  } else {
    seed <- runif(1) * 10^8
    set.seed(seed)
  }
  
  # If a pathway is provided, subset counts.
  if(!is.null(pathway)) {
    mean_expr_1 <- mean(as.matrix(dplyr::select(counts_pair[[1]], -entrezgene)))
    mean_expr_2 <-  mean(as.matrix(dplyr::select(counts_pair[[2]], -entrezgene)))
    counts_pair <- lapply(counts_pair, function(counts) {
      counts <- dplyr::filter(counts, entrezgene %in% pathway)
      return(counts)
    })
  }
  if(any(sapply(counts_pair, nrow) < 2)) {
    warning("Fewer than two genes are expressed in one group. Returning NULL.")
    return(NULL)
  }
  
  ## Join the two datasets and relabel samples
  counts_joined <- dplyr::full_join(counts_pair[[1]], counts_pair[[2]], "entrezgene")
  n <- c(ncol(counts_pair[[1]]) - 1, ncol(counts_pair[[2]]) - 1)
  colnames(counts_joined) <- c("entrezgene", 
                               paste(names(counts_pair)[1], 1:n[1], sep = "_"),
                               paste(names(counts_pair)[2], 1:n[2], sep = "_"))
  
  # Format data as n by p matrix. Store entrezgene id as column names.
  gene_names <- drop(counts_joined$entrezgene)
  counts_joined <- t(dplyr::select(counts_joined, -entrezgene))
  colnames(counts_joined) <- gene_names
  # Fill in missing values with 0
  counts_joined[which(is.na(counts_joined))] <- 0
  
  # Center counts in each group (used for t-tests).
  t_counts <- counts_joined
  t_counts[1:n[1], ] <- t_counts[1:n[1], ] - mean_expr_1
  t_counts[-(1:n[1]), ] <- t_counts[-(1:n[1]), ] - mean_expr_2
  
  # Calculate mean expression of the pathway in each group
  mean_expression_1 <- mean(counts_joined[1:n[1], ])
  mean_expression_2 <- mean(counts_joined[-(1:n[1]), ])
  mean_expression <- mean(counts_joined)
  t_test_pathway <- t.test(apply(t_counts[1:n[1], ], 1, mean),
                           apply(t_counts[-(1:n[1]), ], 1, mean))$p.val
  
  # Calculate mean expression of each gene in each group.
  mean_expression_genes_1 <- apply(counts_joined[1:n[1], ], 2, mean)
  mean_expression_genes_2 <- apply(counts_joined[-(1:n[1]), ], 2, mean)
  mean_expression_genes <- apply(counts_joined, 2, mean)
  t_counts <- counts_joined
  t_counts[1:n[1], ] <- t_counts[1:n[1], ] - mean(t_counts[1:n[1], ])
  t_counts[-(1:n[1]), ] <- t_counts[-(1:n[1]), ] - mean(t_counts[-(1:n[1]), ])
  t_test_genes <- apply(t_counts, 2, function(x) {
    y = x[-(1:n[1])]
    x = x[1:n[1]]
    if(sd(x) == 0 | sd(y) == 0) return(NA)
    t.test(x, y)$p.val
  })
  
  N <- nrow(counts_joined) # Total number of observations across both samples.
  p <- ncol(counts_joined) # Total number of genes in union.
  n_lp <- length(lp_set)
  
  # If total number of observations is small, use exact set of permutations. 
  # Order of samples doesn't matter, so keep only first half if sample sizes
  # are the same. If sample sizes are unequal, all permutations are used.
  if(n[1] == n[2] & choose(N, n[1]) / 2 <= n_perm) {
    permutations <- combn(1:N, n[1])
    permutations <- permutations[, 1:(ncol(permutations) / 2)] 
    permutations <- permutations[, -1]
    # cat("Iterating over all", ncol(permutations), "permutations.\n")
  } else if(n[1] != n[2] & choose(N, n[1]) <= n_perm) {
    permutations <- combn(1:N, n[1])
    permutations <- permutations[, -1]
    # cat("Iterating over all", ncol(permutations), "permutations.\n")
  } else {
    permutations <- cbind(sapply(1:n_perm, function(x) sample(1:N, n[1])))
    # cat("Iterating over", ncol(permutations), "random permutations.\n")
  }
  
  n_perm <- ncol(permutations)
  
  # Initialize lists to store DC scores and p-values for each lp in lp_set.
  # Scores are initialized to 0 and p-values to 1.
  scores_1 <- network_inference(counts_joined[1:n[1], ])
  scores_2 <- network_inference(counts_joined[-(1:n[1]), ])
  
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
  
  # Save inferred networks.
  nw1 <- scores_1
  nw2 <- scores_2
  
  # If inferred networks are empty, end early without performing perm test.
  if(all(scores_1 == 0) && all(scores_2 == 0)) {
    cat("Inferred network is empty. Permutation test not performed.\n")
  } else {
    for(i in 1:n_perm) {
      network_1 <- permutations[, i]
      
      scores_1 <- network_inference(counts_joined[network_1, ])
      scores_1[which(is.na(scores_1))] <- 0
      scores_2 <- network_inference(counts_joined[-network_1, ])
      scores_2[which(is.na(scores_2))] <- 0
      
      for(k in 1:n_lp) {
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
  }
  
  # If no permutation testing is done, set p-values to NA.
  if(n_perm == 1) {
    pval_path_list <- lapply(pval_path_list, function(pval) NA)
    pval_gene_list <- lapply(pval_gene_list, function(pval) NA)
    pval_edge_list <- lapply(pval_edge_list, function(pval) NA)
    
    pval_path_list_mono <- lapply(pval_path_list_mono, function(pval) NA)
    pval_gene_list_mono <- lapply(pval_gene_list_mono, function(pval) NA)
    pval_edge_list_mono <- lapply(pval_edge_list_mono, function(pval) NA)
  } else {
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
  
  results <- vector("list", n_lp)
  for(i in 1:n_lp) {
    results[[i]] <- list(entrezgene = gene_names,
                         p_value_path = pval_path_list[[i]],
                         p_value_genes = pval_gene_list[[i]],
                         p_value_edges = pval_edge_list[[i]],
                         p_value_path_mono = pval_path_list_mono[[i]],
                         p_value_genes_mono = pval_gene_list_mono[[i]],
                         p_value_edges_mono = pval_edge_list_mono[[i]],
                         d_pathway = d_path_list[[i]],
                         d_genes = d_gene_list[[i]],
                         d_edges = d_edge_list[[i]],
                         n_perm = n_perm,
                         n_genes_in_pathway = length(unique(pathway)),
                         n_genes_in_counts = p,
                         mean_expression = mean_expression,
                         mean_expression_1 = mean_expression_1,
                         mean_expression_2 = mean_expression_2,
                         t_test_pathway = t_test_pathway,
                         mean_expression_genes = mean_expression_genes,
                         mean_expression_genes_1 = mean_expression_genes_1,
                         mean_expression_genes_2 = mean_expression_genes_2,
                         t_test_genes = t_test_genes,
                         pathway = pathway_name,
                         lp = lp_set[i],
                         seed = seed)
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
                                        signif = 0.1, monotonized = FALSE) {
  if(is.null(gene_list)) {
    gene_list <- unique(unlist(lapply(results_list, 
                                      function(results) results$entrezgene)))
    gene_list <- sort(gene_list)
  }
  
  # signif <- get_signif_level(pathway_list, signif_overall = signif...
  
  p <- length(gene_list)
  dnw <- matrix(0, p, p)
  colnames(dnw) <- gene_list
  
  for(results in results_list) {
    # Compute the differential network for given pathway.
    dnw_p <- summarize_dnw_for_pathway(results, 
                                       signif = signif, 
                                       monotonized = monotonized)
    
    # Index the columns so they can be aligned with the overall network.
    index <- sapply(results$entrezgene, 
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
summarize_dnw_for_pathway <- function(results, signif = 0.1, monotonized = FALSE) {
  genes_in_dnw <- results$entrezgene
  p <- length(genes_in_dnw)
  
  if(p == 0) {
    warning("results pathway does not contain any genes.")
    return(NULL)
  }
  
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
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used. 
#' @return A vector of DC gene scores, with nonzero entries containing
#' the significant DC scores.
summarize_genes_over_pathways <- function(results_list, gene_list = NULL, 
                                          signif = 0.1, monotonized = FALSE) {
  if(is.null(gene_list)) {
    gene_list <- unique(unlist(lapply(results_list, 
                                      function(results) results$entrezgene)))
    gene_list <- sort(gene_list)
  }
  
  p <- length(gene_list)
  
  df <- tibble(entrezgene = gene_list,
               d_genes = 0,
               p_value_genes = 1,
               p_value_genes_mono = 1,
               n_pathways = 0,
               n_significant = 0,
               n_significant_mono = 0,
               mean_expression = 0,
               mean_expression_1 = 0,
               mean_expression_2 = 0,
               t_test = 1,
               pathways = "")
  
  # If pathway names are not provided in results_list, create them.
  if(is.null(results_list[[1]]$pathway)) {
    for(i in 1:length(results_list)) {
      results_list[[i]]$pathway <- names(results_list)[i]
    }
  }
  
  for(results in results_list) {
    df_p <- summarize_genes_for_pathway(results,
                                        signif = signif,
                                        monotonized = monotonized)
    index <- sapply(df_p$entrezgene, function(gene) which(gene_list == gene)[1])
    
    # Keep largest DC score for each gene. 
    df[index, ]$d_genes <- ifelse(df_p$d_genes >= df[index, ]$d_genes, 
                                  df_p$d_genes, 
                                  df[index, ]$d_genes)
    # Keep smallest p-value.
    df[index, ]$p_value_genes <- ifelse(df_p$p_value_genes <= 
                                          df[index, ]$p_value_genes, 
                                         df_p$p_value_genes, 
                                         df[index, ]$p_value_genes)
    # Keep smallest p-value.
    df[index, ]$p_value_genes_mono <- ifelse(df_p$p_value_genes_mono <= 
                                               df[index, ]$p_value_genes_mono, 
                                        df_p$p_value_genes_mono, 
                                        df[index, ]$p_value_genes_mono)
    
    # Keep track of how many pathways the gene appears in.
    df[index, ]$n_pathways <- df[index, ]$n_pathways + 1
    
    # Keep track of how often (# of pathways) gene is significantly DC.
    df[index, ]$n_significant <- df[index, ]$n_significant + 
      (df_p$p_value_genes <= signif)
    df[index, ]$n_significant_mono <- df[index, ]$n_significant_mono + 
      (df_p$p_value_genes_mono <= signif)
    
    df[index, ]$pathways <- paste0(df[index, ]$pathways, "; ", results$pathway)
    
    df[index, ]$mean_expression <- df_p$mean_expression
    df[index, ]$mean_expression_1 <- df_p$mean_expression_1
    df[index, ]$mean_expression_2 <- df_p$mean_expression_2
    
    df[index, ]$t_test <- df_p$t_test
  }
  
  return(df)
}

#' Summarize differential connectivity scores for a pathway
#' 
#' @param results The results from performing [dna()] on a single pathway.
#' @param signif The desired significance level. DC scores with p-values below
#' this threshold will be set to zero.
#' @param monotonized If true, monotonized (i.e. step-down) p-values from the 
#' permutation test will be used.
#' @return A vector of p DC scores, with nonzero entries containing
#' the significant DC scores.
summarize_genes_for_pathway <- function(results, signif = 0.1, monotonized = FALSE) {
  gene_names <- results$entrezgene
  d_genes <- drop(results$d_genes)
  p_value_genes <- drop(results$p_value_genes)
  p_value_genes_mono <- drop(results$p_value_genes_mono)
  
  mean_expression <- drop(results$mean_expression_genes)
  mean_expression_1 <- drop(results$mean_expression_genes_1)
  mean_expression_2 <- drop(results$mean_expression_genes_2)
  
  t_test <- drop(results$t_test_genes)
  
  return(tibble(entrezgene = gene_names,
                d_genes = d_genes,
                p_value_genes = p_value_genes,
                p_value_genes_mono = p_value_genes_mono,
                mean_expression = mean_expression,
                mean_expression_1 = mean_expression_1,
                mean_expression_2 = mean_expression_2,
                t_test = t_test))
}



#' Summarize differential connectivity scores from different pathways
#' 
#' @param results_list a list of results from performing dna() on several pathways.
#' @param gene_list a set of entrezgene IDs to be summarized.
#' @param signif a significance threshold for differential connectivity.
#' @return list of results
summarize_pathways <- function(results_list, signif = 0.1) {
  # Remove any NULL results from the list.
  results_list <- results_list[!sapply(results_list, is.null)]
  
  pathway_names <- names(results_list)
  p_values <- 
    sapply(results_list, function(results) results$p_value_path)
  d_pathway <- 
    sapply(results_list, function(results) results$d_pathway)
  n_genes <- 
    sapply(results_list, function(results) results$n_genes_in_pathway)
  n_genes_expressed <- 
    sapply(results_list, function(results) results$n_genes_in_counts)
  n_genes_significant <- 
    sapply(results_list, function(results) sum(results$p_value_genes_mono <= signif))
  mean_expression <- 
    sapply(results_list, function(results) results$mean_expression)
  mean_expression_1 <- 
    sapply(results_list, function(results) results$mean_expression_1)
  mean_expression_2 <- 
    sapply(results_list, function(results) results$mean_expression_2)
  t_test <-
    sapply(results_list, function(results) results$t_test_pathway)
  
  df <- tibble(pathway = pathway_names,
               p_values = p_values,
               d_pathway = d_pathway,
               n_genes = n_genes,
               n_genes_expressed = n_genes_expressed,
               n_genes_significant = n_genes_significant,
               prop_genes_expressed = n_genes_expressed / n_genes,
               prop_genes_significant = n_genes_significant / n_genes,
               mean_expression = mean_expression,
               mean_expression_1 = mean_expression_1,
               mean_expression_2 = mean_expression_2,
               t_test = t_test)
  
  return(df)
}



#' Summarize differential connectivity scores within a pathway
#' 
#' @param results the result from dna().
#' @param gene_list a set of entrezgene IDs to be summarized.
#' @param signif a significance threshold for differential connectivity.
#' @return list of results
summarize_genes <- function(results, gene_list = NULL, signif = 0.1) {
  if(is.null(gene_list)) {
    gene_list <- results$entrezgene
  }
  
  d_genes <- results$d_genes
  p_values <-results$p_value_genes
  n_perm <- results$n_perm
  
  df <- tibble(entrezgene = gene_list,
               d_genes = d_genes,
               p_values = p_values,
               n_perm = n_perm)
  
  return(df)
}