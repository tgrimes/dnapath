i <- 3
network_number <- 1
seed_network <- 10000000 + 1000 * i + network_number
network <- get_networks()

nw1 <- network$`1`
nw1_adj <- get_adj_matrix_from_network(nw1)
paths <- get_pathway_list(nw1)

nw2 <- network$`2`
nw2_adj <- get_adj_matrix_from_network(nw2)
paths <- get_pathway_list(nw2)

set.seed(1)
png(paste0("network_diff.png"), height = 680, width = 680, res = 100)
par(mar = c(0, 0, 1, 0))
plot_network_diff(nw1, nw2, main = "Differential network", edge_scale = 1.5,
                  as_subgraph = FALSE, include_vertex_labels = TRUE)
dev.off()

for(i in 1:length(paths)) {
  png(paste0("network_diff_pathway_", i, ".png"), res = 180)
  par(mar = c(0, 0, 1, 0))
  plot_network_diff(nw1_adj[paths[[i]], paths[[i]]],
                    nw2_adj[paths[[i]], paths[[i]]],
                    as_subgraph = FALSE,
                    include_vertex_labels = TRUE, edge_scale = 1.5,
                    main = paste("Pathway", i))
  dev.off()
}
