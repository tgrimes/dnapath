if(!exists("seed_network")) {
  stop("Variable `seed_network` needs to be specified.")
}
if(!exists("output_dir")) {
  stop("Variable `output_dir` needs to be specified.")
}

save_file <- paste0(output_dir, "figures/network_diff_", seed_network, ".png")
if(file.exists(save_file)) {
  stop(paste("plot already exists for this network seed at", save_file))
}

# Plot the two networks.
network <- get_networks()

cat("... saving network differential plot to", save_file, "\n")
png(filename = save_file,
    width = 1680,
    height = 1680,
    res = 300)
par(mfrow = c(1, 1),
    mar = c(0, 0, 1.5, 0))
plot_network_diff(network[[1]], network[[2]], include_vertex_labels = TRUE,
                  main = "Graphical representation of the differential networks", 
                  as_subgraph = FALSE, node_scale = 10, edge_scale = 1.2,
                  label_scale = 0.65)
dev.off()
