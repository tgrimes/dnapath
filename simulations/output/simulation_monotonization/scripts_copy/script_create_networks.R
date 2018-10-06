if(!exists("network_params")) {
  stop("Variable `network_params` needs to be specified.")
}
if(!exists("output_dir")) {
  stop("Variable `output_dir` needs to be specified.")
}
if(!exists("seed_network")) {
  stop("Variable `seed_network` needs to be specified.")
}

save_file <- paste0(output_dir, "generated_networks/",
                    "network_seed_", seed_network, ".rds")
if(file.exists(save_file)) {
  stop(paste("network already exists for this seed at", save_file))
}


set.seed(seed_network)
p <- network_params$p
n_pathways <- network_params$n_pathways
mu <- network_params$mu
sigma <- network_params$sigma
size <- mu^2 / (sigma^2 - mu)
pathways <- lapply(rnbinom(n_pathways, mu = mu, size = size), function(x) {
  sample(1:p, ifelse(x < 2, 2, x))
})

# Create network structures for the simulation.
cat("Creating network with p =", p, "nodes. Seed =", seed_network, "\n")
network <- create_network(p = p,
                          modules = pathways,
                          module_neig = 1,
                          module_prob = 0.15)
network$seed <- seed_network
pathway_nodes <- unique(unlist(sapply(network$modules, function(x) x$nodes)))

# Turn on random hub nodes:
n_hubs <- c(3, 6, 9)[which(n_pathways == c(4, 10, 20))]
index <- sample(pathway_nodes, n_hubs)
n_hubs <- length(index)
for(i in index) {
  network <- rewire_connections_to_node(network, i, 0.5)
}
for(i in 1:length(network$modules)) {
  network <- remove_small_components_from_module(network, i)
}

# sapply(network$modules, function(x) {length(x$nodes)})
# plot(network)
# 
# par(mfrow = c(ceiling(sqrt(n_pathways)), 
#               ceiling(n_pathways / ceiling(sqrt(n_pathways)))))
# sapply(1:n_pathways, function(i) plot_network(network$modules[[i]]$struct, main = paste("Pathway", i)))
# par(mfrow = c(1, 1))



# cat("Turning off one hub in network.\n")
network2 <- network
# Modify connections to various nodes:
off_hubs <- sample(index, round(n_hubs / 3)) # Turn off 1/3 of hubs.
rewire_hubs <- sample(setdiff(index, off_hubs), round(n_hubs / 3)) # Rewire 1/3.
rewire_other <- sample(setdiff(pathway_nodes, index), round(p * 0.025)) 
for(i in off_hubs) {
  network2 <- remove_connections_to_node(network2, i)
}
for(i in rewire_hubs) {
  network2 <- rewire_connections_to_node(network2, i, 0.5)
}
for(i in rewire_other) {
  network2 <- rewire_connections_to_node(network2, i, 0.02)
}


network_list <- list(`1` = network,
                     `2` = network2)
cat("... saving list of networks to", save_file, "\n")
saveRDS(network_list, save_file)

# plot(network2)
# plot_network_diff(network, network2)
# 
# par(mfrow = c(ceiling(sqrt(n_pathways)), 
#               ceiling(n_pathways / ceiling(sqrt(n_pathways)))))
# sapply(1:n_pathways, function(i) {
#   plot_network_diff(network$modules[[i]]$struct, network2$modules[[i]]$struct, 
#                     main = paste("Pathway", i))
# })
# par(mfrow = c(1, 1))
# 
# sort(apply(get_adj_matrix_from_network(network) != get_adj_matrix_from_network(network2), 2, sum),
#      decreasing = TRUE)
