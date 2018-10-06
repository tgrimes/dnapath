entrez_to_reactome <- function(df, 
                               species = c("human", "mouse"),
                               dir_annotation = "data/annotations/",
                               G = NULL) {
  ##############################################################
  ##############################################################
  warning("entrez_to_reactome is under construction. Returning df.")
  return(df)
  ##############################################################
  ##############################################################
  
  if(!("entrezgene" %in% colnames(df)))
    stop("df must contain a column called `entrezgene` containing entrez gene id's.")
  col_index_entrez <- which(colnames(df) == "entrezgene")
  species <- tolower(species[1])
  if(!(species %in% c("human", "mouse"))) 
    stop("species must be either 'human' or 'mouse'.")
  
  if(is.null(G)) {
    G <- get_reactome_matrix(species)
  }
  
  entrez <- drop(df[, col_index_entrez]) # Extract entrez id's as vector.
  pathway <- vector("character", length(entrez))
  
  index <- which(rownames(G) %in% entrez)
  G <- G[index, ]
  index <- which(entrez %in% rownames(G))

  names <- rownames(G)
  for(i in 1:nrow(G)) {
    pathway[entrez == names[i]] <- paste(colnames(G)[G[i, ] == 1], collapse = "; ")
  }
  
  df$pathways <- rep("", nrow(df))
  df$pathways <- pathway
  
  if(!is.numeric(df$entrezgene)) df$entrezgene <- as.numeric(df$entrezgene)
  
  return(df)
}
