entrez_to_goid <- function(df, species,
                           dir_annotation = "data/annotations/") {
  if(!any(names(df) %in% c("entrezgene")))
    stop("df must contain a column nammed 'entrezgene'.")
  species <- pmatch(tolower(species), c("human", "mouse"))
  if(is.na(species))
    stop("species must be either 'human' or 'mouse'.")
  
  if(species == 1) {
    # Human
    if(!("mart_human" %in% ls())) init_mart_human()
    gene_info <- getBM(attributes = c("entrezgene", "go_id"), 
                       mart = mart_human)
  } else {
    # Mouse
    if(!("mart_mouse" %in% ls())) init_mart_mouse()
    gene_info <- getBM(attributes = c("entrezgene", "go_id"), 
                       mart = mart_mouse)
  }
  
  gene_info %>%
    filter(go_id != "", entrezgene != "") ->
    gene_info
  
  if(!is.numeric(df$entrezgene)) df$entrezgene <- as.numeric(df$entrezgene)
  
  return(left_join(df, gene_info, by = "entrezgene"))
}