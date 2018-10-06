source("src/init_required_packages_and_files.R")
source("src/scripts_neuroblastoma/plot_network_diff.R")

dir_output <- "./output/neuroblastoma/" 
filtered_counts_name <- "zeroes0.8"
dir_filtered <- paste0(dir_output, filtered_counts_name, "/")

measure <- "pcor"
lp <- 2
dir_output <- dir_filtered
signif <- 0.01
monotonized <- FALSE

load_file <- paste0(dir_output, "dna_results/", measure, "/lp", lp, "_",
                    "neuro.rds")
results_list <- readRDS(load_file)
# results_list <- results_list[summarize_pathways(results_list, signif = signif)$prop_genes_expressed >= 0.8]

index_res <- which(sapply(results_list, function(r) {
  r$p_value_path <= signif &
    any(r$p_value_edges_mono <= signif) &
    all(c("9475") %in% r$entrezgene)
  }))

cat("Pathways:\n")
pathway_name <- "InlA-mediated entry of Listeria monocytogenes into host cells"
results <- results_list[[pathway_name]]

# ii <- 1
# results <- results_list[[index_res[ii]]] # 5,
# pathway_name <- names(results_list)[index_res[ii]]
results
df <- tibble(entrezgene = results$entrezgene,
             d_gene = drop(results$d_gene),
             p_value = drop(results$p_value_genes),
             p_value_mono = drop(results$p_value_genes_mono))
df <- entrez_to_symbol(df, "human")

dnw <- summarize_dnw_for_pathway(results, monotonized = monotonized)
colnames(dnw) <- df$hgnc_symbol

df_genes <- summarize_genes_for_pathway(results) %>%
  entrez_to_symbol("human") %>%
  entrez_to_description("human")
df_genes %>%
  dplyr::mutate(DC_score = round(d_genes, 2),
                DC_pval = round(p_value_genes, 2),
                DC_pval_mono = round(p_value_genes_mono, 2),
                DE_pval = round(t_test, 2),
                ant = round(mean_expression_1, 2),
                post = round(mean_expression_2, 2),
                FC = round(ifelse(ant > post, ant / post, post / ant), 2)) %>%
  dplyr::select(hgnc_symbol, DC_score, DC_pval, DC_pval_mono, DE_pval, ant, post, FC, description)
## View()


counts_list <- readRDS(paste0(dir_filtered, "zeroes0.8"))
counts_list <- lapply(counts_list, function(counts) {
  hgnc <- entrez_to_symbol(counts, "human")$hgnc_symbol
  index <- which(hgnc %in% df$hgnc_symbol)
  counts <- t(counts[index, -1])
  colnames(counts) <- hgnc[index]
  return(counts)
})

nw <- vector("list", 2)
for(i in 1:2) {
  nw[[i]] <- matrix(0, nrow(dnw), ncol(dnw))
  index <- which(df_genes$hgnc_symbol %in% colnames(counts_list[[i]]))
  nw[[i]][index, index] <- pcor_shrinkC(counts_list[[i]])
  colnames(nw[[i]]) <- df_genes$hgnc_symbol
}
nw1 <- nw[[1]]
nw2 <- nw[[2]]
nw1[dnw == 0] <- 0
nw2[dnw == 0] <- 0

more_active_in_nw1 <- which((nw1 - nw2) > 0)
more_active_in_nw2 <- which((nw2 - nw1) > 0)
nw1[more_active_in_nw2] <- 0
nw2[more_active_in_nw1] <- 0


set.seed(5)
g <- plot_network(nw1 + nw2, main = "DNW")
coords <- igraph::layout.fruchterman.reingold(g)
node_color <- c(adjustcolor("orange", 0.5),
                adjustcolor("red", 0.5))[1 + (df_genes$mean_expression_2 > df_genes$mean_expression_1)]
node_weights <- apply(rbind(df_genes$mean_expression_1 / df_genes$mean_expression_2,
                            df_genes$mean_expression_2 / df_genes$mean_expression_1),
                      2, max) * (df_genes$t_test < 0.1)
node_weights[node_weights == 0] <- 1
node_weights[is.na(node_weights)] <- 1
node_weights[node_weights == Inf] <- 1
node_weights[is.nan(node_weights)] <- 1
node_weights <- node_weights^2 + 1
edge_weights <- dnw[lower.tri(dnw)][dnw[lower.tri(dnw)] != 0]
edge_weights <- edge_weights / max(edge_weights)


png("InlA_mediated_entry.png",
    height = 1080, width = 1080, res = 300)
par(mar = c(0, 0, 2.5, 0))
plot_network_diff(nw2, nw1, coords = coords, 
                  node_color = node_color, 
                  node_scale = 2, node_weights = node_weights,
                  edge_scale = 3, edge_weights = edge_weights,
                  main = paste0("InlA-mediated entry of Listeria\nmonocytogenes into host cells"),
                  as_subgraph = FALSE, include_vertex_labels = TRUE)
dev.off()
