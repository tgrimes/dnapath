symbol_to_entrez <- function(df, species) {
  if(!any(names(df) %in% c("hgnc_symbol", "mgi_symbol")))
    stop("df must contain a column nammed 'hgnc_symbol' or 'mgi_symbol'.")
  species <- pmatch(tolower(species), c("human", "mouse"))
  if(is.na(species))
    stop("species must be either 'human' or 'mouse'.")
  
  if(species == 1) {
    # ! TODO: use instead -> as.list(org.Hs.egSYMBOL)
    # Human
    if(file.exists("data/annotations/hgnc_to_entrez.rds")) {
      load_file <- "data/annotations/hgnc_to_entrez.rds"
      cat("\t- loading gene info from", load_file, "\n")
      gene_info <- readRDS(load_file)
    } else {
      cat("Obtaining gene info from biomart.\n")
      
      if(!("mart_human" %in% ls())) init_mart_human()
      gene_info <- getBM(attributes = c("hgnc_symbol", "entrezgene"), 
                         mart = mart_human)
      gene_info %>%
        group_by(hgnc_symbol) %>%
        summarise(entrezgene = first(entrezgene)) ->
        gene_info
      
      gene_info$hgnc_symbol <- gsub("-", "_", gene_info$hgnc_symbol)
      
      save_file <- "data/annotations/hgnc_to_entrez.rds"
      cat("\t- saving gene info to", save_file, "\n")
      saveRDS(gene_info, save_file)
    }
    
  } else {
    # Mouse
    if(file.exists("data/annotations/mgi_to_entrez.rds")) {
      load_file <- "data/annotations/mgi_to_entrez.rds"
      cat("\t- loading gene info from", load_file, "\n")
      gene_info <- readRDS(load_file)
    } else {
      cat("Obtaining gene info from biomart.\n")
      
      if(!("mart_mouse" %in% ls())) init_mart_mouse()
      gene_info <- getBM(attributes = c("mgi_symbol", "entrezgene"), 
                         mart = mart_mouse)
      gene_info %>%
        group_by(mgi_symbol) %>%
        summarise(entrezgene = first(entrezgene)) ->
        gene_info
      
      gene_info$mgi_symbol <- gsub("-", "_", gene_info$mgi_symbol)
      
      save_file <- "data/annotations/mgi_to_entrez.rds"
      cat("\t- saving gene info to", save_file, "\n")
      saveRDS(gene_info, save_file)
    }
    
  }
  
  return(left_join(df, gene_info))
}