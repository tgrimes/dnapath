source("src/init_required_packages_and_files.R")
source("src/scripts_facebase/plot_network_diff.R")

measure = "pcor"
lp <- 2
pair <- 1
pairs <- list(c(1, 5), c(2, 6), c(3, 7), c(4, 8))[[pair]]
pair_name <- c("Lateral", "Medial", "Nasal", "Oral")[pair]
dir_filtered <- "output/facebase/zeroes0.34_normalized/"
dir_output <- dir_filtered
signif = 0.1
monotonized <- FALSE

load_file <- paste0(dir_output, "dna_results/", measure, "/lp", lp, "_",
                    pair_name, ".rds")
results_list <- readRDS(load_file)

cat("Pathways:\n")
pathway_name <- "Signaling by NOTCH1"
results <- results_list[[pathway_name]]
results
df <- tibble(entrezgene = results$entrezgene,
             d_gene = drop(results$d_gene),
             p_value = drop(results$p_value_genes),
             p_value_mono = drop(results$p_value_genes_mono))
df <- entrez_to_symbol(df, "mouse")

dnw <- summarize_dnw_for_pathway(results, monotonized = monotonized)
colnames(dnw) <- df$mgi_symbol

df_genes <- summarize_genes_for_pathway(results) %>%
  entrez_to_symbol("mouse") %>%
  entrez_to_description("mouse")
df_genes %>%
  dplyr::mutate(DC_score = round(d_genes, 2),
                DC_pval = round(p_value_genes, 2),
                DC_pval_mono = round(p_value_genes_mono, 2),
                DE_pval = round(t_test, 2),
                ant = round(mean_expression_1, 2),
                post = round(mean_expression_2, 2),
                FC = round(ifelse(ant > post, ant / post, post / ant), 2)) %>%
  dplyr::select(mgi_symbol, DC_score, DC_pval, DC_pval_mono, DE_pval, ant, post, FC, description)


counts_list <- readRDS(paste0(dir_filtered, "zeroes0.34_normalized"))
counts_list <- lapply(counts_list, function(counts) {
  mgi <- entrez_to_symbol(counts, "mouse")$mgi_symbol
  index <- which(mgi %in% df$mgi_symbol)
  counts <- t(counts[index, -1])
  colnames(counts) <- mgi[index]
  return(counts)
})

nw <- vector("list", 2)
for(i in 1:2) {
  nw[[i]] <- matrix(0, nrow(dnw), ncol(dnw))
  index <- which(df_genes$mgi_symbol %in% colnames(counts_list[[pairs[i]]]))
  nw[[i]][index, index] <- pcor_shrinkC(counts_list[[pairs[i]]])
  colnames(nw[[i]]) <- df_genes$mgi_symbol
}
nw1 <- nw[[1]]
nw2 <- nw[[2]]
nw1[dnw == 0] <- 0
nw2[dnw == 0] <- 0

more_active_in_nw1 <- which((nw1 - nw2) > 0)
more_active_in_nw2 <- which((nw2 - nw1) > 0)
nw1[more_active_in_nw2] <- 0
nw2[more_active_in_nw1] <- 0

set.seed(361473)
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

png("signaling_by_notch1.png", 
    height = 1080, width = 1080, res = 300)
par(mar = c(0, 0, 1.5, 0))
plot_network_diff(nw1, nw2, coords = coords, 
                  node_color = node_color, 
                  node_scale = 2, node_weights = node_weights,
                  edge_scale = 3, edge_weights = edge_weights,
                  main = "Signaling by NOTCH1",
                  as_subgraph = FALSE, include_vertex_labels = TRUE)
dev.off()
