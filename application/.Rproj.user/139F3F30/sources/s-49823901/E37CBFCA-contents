add_expression_counts <- function(df, counts_pair) {
  if(length(counts_pair) != 2) stop("counts_pair should be of length 2.")
  if(colnames(counts_pair[[1]])[1] != "entrezgene")
    stop("first column of counts_pair[[1]] is not entrezgene.")
  if(colnames(counts_pair[[2]])[1] != "entrezgene")
    stop("first column of counts_pair[[2]] is not entrezgene.")
  
  group_names <- names(counts_pair)
  if(length(unique(group_names)) != 2) {
    group_names <- c("First", "Second")
  }
  
  counts_joined <- full_join(counts_pair[[1]], counts_pair[[2]], "entrezgene")
  counts_joined[is.na(counts_joined)] <- 0
  n <- c(ncol(counts_pair[[1]]) - 1, ncol(counts_pair[[2]]) - 1)
  colnames(counts_joined) <- c("entrezgene", 
                               paste(group_names[1], 1:n[1], sep = "_"),
                               paste(group_names[2], 1:n[2], sep = "_"))
  counts_joined <- counts_joined %>%
    mutate(mean_expression_1 = dplyr::select(., starts_with(group_names[1])) %>% apply(1, mean),
           mean_expression_2 = dplyr::select(., starts_with(group_names[2])) %>% apply(1, mean),
           mean_expression = dplyr::select(., starts_with(group_names[1]), 
                                           starts_with(group_names[2])) %>% apply(1, mean),
           t_test = dplyr::select(., starts_with(group_names[1]), 
                                  starts_with(group_names[2])) %>% 
             apply(1, function(x) {
               y = x[(n[1] + 1):(n[1] + n[2])]
               x = x[1:n[1]]
               if(sd(x) == 0 | sd(y) == 0) return(NA)
               t.test(x, y)$p.val
             })) %>%
    dplyr::select(entrezgene, 
                  mean_expression_1, 
                  mean_expression_2, 
                  mean_expression, 
                  t_test)
  
  df <- left_join(df, counts_joined, by = "entrezgene")
  return(df)
}
