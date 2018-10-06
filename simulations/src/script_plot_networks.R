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
    width = 1280,
    height = 1280,
    res = 256)
par(mfrow = c(1, 1))
plot_network_diff(network[[1]], network[[2]], as_subgraph = TRUE,
                  main = "Network 1 vs Network 2", node_scale = 10)
dev.off()
