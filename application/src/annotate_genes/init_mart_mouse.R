init_mart_mouse <- function() {
  # Obtain MGI gene symbols and other information for each Entrez ID. 
  # listMarts(); listAttributes(mart) # Useful functions to obtain more info.
  # listDatasets(useMart('ensembl'));
  
  # Add global variable.
  mart_mouse <<- useMart(biomart = "ensembl", 
                         dataset = "mmusculus_gene_ensembl")
}

