goid_to_name <- function(df, species) {
  if(!any(names(df) %in% c("go_id")))
    stop("df must contain a column nammed 'go_id'.")
  species <- pmatch(tolower(species), c("human", "mouse"))
  if(is.na(species))
    stop("species must be either 'human' or 'mouse'.")
  
  if(species == 1) {
    # Human
    if(!("mart_human" %in% ls())) init_mart_human()
    gene_info <- getBM(attributes = c("go_id", "name_1006"), 
                       mart = mart_human)
  } else {
    # Mouse
    if(!("mart_mouse" %in% ls())) init_mart_mouse()
    gene_info <- getBM(attributes = c("go_id", "name_1006"), 
                       mart = mart_mouse)
  }
  
  gene_info %>%
    filter(go_id != "") ->
    gene_info
  
  return(left_join(df, gene_info))
}