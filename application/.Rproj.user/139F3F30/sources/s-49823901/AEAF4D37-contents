# Modified function from SeqNet package:
plot_network_diff <- function (network_1, network_2, as_subgraph = TRUE, 
          node_color = adjustcolor("orange", 0.5), node_scale = 5, node_weights = NULL,
          edge_weights = NULL, edge_scale = 1, 
          coords = NULL, main = "Untitled", include_vertex_labels = FALSE, 
          ...) 
{
  if (class(network_1) == "network") {
    adj_matrix_1 <- get_adj_matrix_from_network(network_1)
  }
  else if (is.matrix(network_1)) {
    if (!all(network_1 %in% c(0, 1))) {
      network_1[network_1 != 0] <- 1
    }
    adj_matrix_1 <- network_1
  }
  else {
    stop("network_1 should be either a matrix or a network class object.")
  }
  if (class(network_2) == "network") {
    adj_matrix_2 <- get_adj_matrix_from_network(network_2)
  }
  else if (is.matrix(network_2)) {
    if (!all(network_2 %in% c(0, 1))) {
      network_2[network_2 != 0] <- 1
    }
    adj_matrix_2 <- network_2
  }
  else {
    stop("network_2 should be either a matrix or a network class object.")
  }
  if (is.null(colnames(adj_matrix_1))) {
    colnames(adj_matrix_1) <- 1:ncol(adj_matrix_1)
  }
  if (is.null(colnames(adj_matrix_2))) {
    colnames(adj_matrix_2) <- 1:ncol(adj_matrix_2)
  }
  if (all(adj_matrix_1 == 0) && all(adj_matrix_2 == 0)) {
    plot(0, xaxt = "n", yaxt = "n", bty = "n", pch = "", 
         ylab = "", xlab = "", main = main)
    return(NULL)
  }
  adj_matrix_both <- (adj_matrix_1 | adj_matrix_2) * 1
  if (as_subgraph) {
    degree_zero <- apply(adj_matrix_both, 2, sum) == 0
    adj_matrix_both <- adj_matrix_both[!degree_zero, !degree_zero]
    adj_matrix_1 <- adj_matrix_1[!degree_zero, !degree_zero]
    adj_matrix_2 <- adj_matrix_2[!degree_zero, !degree_zero]
  }
  g <- igraph::graph_from_adjacency_matrix(adj_matrix_both, 
                                           mode = "undirected", 
                                           weighted = NULL)
  
  if(is.null(node_weights)) {
    igraph::V(g)$size <- log(apply(abs(adj_matrix_1 - adj_matrix_2), 
                                   2, sum) + 1) + 1
    igraph::V(g)$size <- igraph::V(g)$size/max(igraph::V(g)$size) * 
      node_scale
  } else {
    igraph::V(g)$size <- node_weights * node_scale
  }
  
  igraph::V(g)$frame.color <- "white"
  
  if(is.null(edge_weights)) {
    edge_weights <- rep(edge_scale, length(igraph::E(g)))
  } else {
    edge_weights <- edge_weights * edge_scale
  }
  g_1 <- igraph::graph_from_adjacency_matrix(adj_matrix_1, 
                                             mode = "undirected", 
                                             weighted = NULL)
  g_2 <- igraph::graph_from_adjacency_matrix(adj_matrix_2, 
                                             mode = "undirected", 
                                             weighted = NULL)
  if (include_vertex_labels) {
    vertex.label.color <- rgb(0.1, 0.1, 0.1, 0.8)
  }
  else {
    vertex.label.color <- rgb(0, 0, 0, 0)
  }
  edge.color <- vector("character", length(igraph::E(g)))
  subset_1 <- which(attr(igraph::E(g), "vnames") %in% attr(igraph::E(g_1), 
                                                           "vnames"))
  edge.color[subset_1] <- "orange"
  subset_2 <- which(attr(igraph::E(g), "vnames") %in% attr(igraph::E(g_2), 
                                                           "vnames"))
  edge.color[subset_2] <- "red"
  edge.color[intersect(subset_1, subset_2)] <- "black"
  if (is.null(coords)) {
    coords <- igraph::layout.fruchterman.reingold(g_1)
  }
  
  plot(g, vertex.color = node_color, vertex.label.font = 2, 
       vertex.label.color = vertex.label.color, vertex.label.cex = 0.7, 
       edge.color = edge.color, edge.width = edge_weights,
       layout = coords, main = main, 
       ...)
  return(g)
}
