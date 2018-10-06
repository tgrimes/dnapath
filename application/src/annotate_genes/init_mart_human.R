init_mart_human <- function() {
  # Obtain MGI gene symbols and other information for each Entrez ID. 
  # listMarts(); listAttributes(mart) # Useful functions to obtain more info.
  # listDatasets(useMart('ensembl'));
  
  # Add global variable.
  mart_human <<- useMart(biomart = "ensembl",
                         dataset = "hsapiens_gene_ensembl")
}

