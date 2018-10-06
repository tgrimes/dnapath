summarize_results_from_dna <- function( dataset = c("facebase", "neuroblastoma"),
                                        measure = c("corr", "pcor"),
                                        lp = 2,
                                        pair_name = c("Lateral", "Medial", "Nasal", "Oral"),
                                        dir_output = c("output/facebase/zeroes0.34_normalized/",
                                                        "neuroblastoma/zeroes0.8/"),
                                        signif = 0.1,
                                        plotting = FALSE) {
  cat(pair_name, "--------------------------------------\n")

  # dataset = c("facebase", "neuroblastoma")
  # measure = c("corr", "pcor")
  # lp = 2
  # pair_name = c("Lateral", "Medial", "Nasal", "Oral")
  # dir_output = c("output/facebase/zeroes0.34_normalized/",
  #                "neuroblastoma/zeroes0.8/")
  # signif = 0.1
  # plotting = FALSE
  
  dataset <- tolower(dataset[1])
  measure <- tolower(measure[1])
  pair_name <- pair_name[1]
  dir_output <- dir_output[1]

  load_file <- paste0(dir_output, "dna_results/", measure, "/lp", lp, "_",
                      ifelse(dataset == "facebase", pair_name, "neuro"),
                      ".rds")
  results_list <- readRDS(load_file)
  
  
  
  #
  # Summarize pathways
  #
  cat("Pathways:\n")
  df_path <- summarize_pathways(results_list, signif = signif)
  
  cat("\t", sum(df_path$prop_genes_expressed < 0.8), "of", nrow(df_path),
      "pathways have <80% of their genes expressed. These are removed.\n")
  df_path <- df_path %>%
    filter(prop_genes_expressed >= 0.8) %>%
    mutate(mean_expr_perc_ant = rank(mean_expression_1) / n(),
           mean_expr_perc_post = rank(mean_expression_2) / n(),
           mean_expr_fold_ant = mean_expression_1 / mean(mean_expression_1),
           mean_expr_fold_post = mean_expression_2 / mean(mean_expression_2))
  
  cat("\tMean expression of pathways =", round(mean(df_path$mean_expression), 2), "\n")
  cat("\t", sum(df_path$p_values <= signif), "of", nrow(df_path), 
      "expressed pathways are significantly DC.\n")
  cat("\t", nrow(filter(df_path, 
                        p_values <= signif, 
                        prop_genes_significant == 0)),
      "pathways with p-value <=", signif, "but no significant genes.\n")
  
  if(plotting) {
    tryCatch({
      g <- df_path %>%
        ggplot(aes(x = n_genes_expressed, y = d_pathway)) + 
        geom_point(alpha = 0.2) + 
        facet_grid(~ (p_values <= signif)) +
        theme_bw() +
        labs(title = paste0(dataset, 
                            ifelse(dataset == "facebase", 
                                   paste0(" (", pair_name, ") "), 
                                   ""), 
                            " - ", measure, " lp", lp))
      print(g)
    }, error = function(e) print(e))
    tryCatch({
      g <- df_path %>%
        ggplot(aes(x = t_test, y = d_pathway)) + 
        geom_point(alpha = 0.2) + 
        facet_grid(~ (p_values <= signif)) +
        theme_bw() +
        labs(title = paste0(dataset, 
                            ifelse(dataset == "facebase", 
                                   paste0(" (", pair_name, ") "), 
                                   ""), 
                            " - ", measure, " lp", lp))
      print(g)
    }, error = function(e) print(e))
    tryCatch({
      labels <- paste0(c("F", "T")[(df_path$p_values <= signif) + 1], 
                    c("F", "T")[(df_path$t_test <= signif) + 1])
      col <- c("None", "DE", "DC", "Both")[sapply(labels, function(lab) 
        which(c("FF", "FT", "TF", "TT") %in% lab))]
      df_path$significance <- factor(col, levels = c("None", "DE", "DC", "Both"))
      g <- df_path %>%
        arrange(significance) %>%
        ggplot(aes(x = mean_expression, y = d_pathway, 
                   color = significance)) +
        scale_color_manual(values = c("black",
                                      "grey",
                                      "orange",
                                      "red")) +
        geom_hline(yintercept = median(df_path$d_pathway)) +
        geom_vline(xintercept = median(df_path$mean_expression)) +
        geom_point(alpha = 0.75) + 
        theme_bw() +
        labs(title = paste0(dataset, 
                            ifelse(dataset == "facebase", 
                                   paste0(" (", pair_name, ") "), 
                                   ""), 
                            " - ", measure, " lp", lp))
      print(g)
    }, error = function(e) print(e))
  }
  
  df_path_signif <- filter(df_path, p_values <= signif)
  
  
  
  #
  # Summarize genes
  #
  cat("Genes:\n")
  # results_list_filtered <- results_list
  results_list_filtered <- results_list[sapply(results_list, function(r) {
    r$n_genes_in_counts / r$n_genes_in_pathway >= 0.8 })]  
  df_genes <- summarize_genes_over_pathways(results_list_filtered, 
                                            signif = signif)
  
  # Add additional annotation to genes (from entrezgene ID).
  df_genes <- df_genes %>%
    entrez_to_symbol(ifelse(dataset == "facebase", "mouse", "human"))
  
  df_genes <- df_genes %>%
    filter(!is.na(mean_expression)) %>%
    mutate(mean_expr_perc_ant = rank(mean_expression_1) / n(),
           mean_expr_perc_post = rank(mean_expression_2) / n(),
           mean_expr_fold_ant = mean_expression_1 / mean(mean_expression_1),
           mean_expr_fold_post = mean_expression_2 / mean(mean_expression_2))
  cat("\tMean expression of genes =", round(mean(df_genes$mean_expression), 2), "\n")
  cat("\t", sum(df_genes$n_significant > 0), "of", nrow(df_genes), 
      "genes significantly DC.\n")
  cat("\t", sum(df_genes$n_significant_mono > 0), "of", nrow(df_genes), 
      "genes significantly DC with monotonization.\n")
  
  if(plotting) {
    tryCatch({
      g <- df_genes %>%
        ggplot(aes(x = mean_expression, y = d_genes)) + 
        geom_point(alpha = 0.2) + 
        facet_grid(~ (n_significant > 0)) +
        theme_bw() +
        labs(title = paste0(dataset, 
                            ifelse(dataset == "facebase", 
                                   paste0(" (", pair_name, ") "), 
                                   ""), 
                            " - ", measure, " lp", lp))
      print(g)
    }, error = function(e) print(e))
    tryCatch({
      g <- df_genes %>%
        ggplot(aes(x = t_test, y = d_genes)) + 
        geom_point(alpha = 0.2) + 
        facet_grid(~ (n_significant > 0)) +
        theme_bw() +
        labs(title = paste0(dataset, 
                            ifelse(dataset == "facebase", 
                                   paste0(" (", pair_name, ") "), 
                                   ""), 
                            " - ", measure, " lp", lp))
      print(g)
    }, error = function(e) print(e))
    tryCatch({
      labels <- paste0(c("F", "T")[(df_genes$n_significant > 0) + 1], 
                       c("F", "T")[(!is.na(df_genes$t_test) &
                                     (df_genes$t_test <= signif)) + 1])
      col <- c("None", "DE", "DC", "Both")[sapply(labels, function(lab) 
        which(c("FF", "FT", "TF", "TT") %in% lab))]
      df_genes$significance <- factor(col, levels = c("None", "DE", "DC", "Both"))
      g <- df_genes %>%
        arrange(significance) %>%
        ggplot(aes(x = mean_expression, y = d_genes, 
                   color = significance)) +
        scale_color_manual(values = c("black",
                                      "grey",
                                      "orange",
                                      "red")) +
        geom_hline(yintercept = median(df_genes$d_genes)) +
        geom_vline(xintercept = median(df_genes$mean_expression)) +
        geom_point(alpha = 0.75) + 
        theme_bw() +
        labs(title = paste0(dataset, 
                            ifelse(dataset == "facebase", 
                                   paste0(" (", pair_name, ") "), 
                                   ""), 
                            " - ", measure, " lp", lp))
      print(g)
    }, error = function(e) print(e))
  }
  
  symbol_name <- ifelse(dataset == "facebase", "mgi_symbol", "hgnc_symbol")
  df_genes_signif <- df_genes %>%
    filter(n_significant_mono > 0) %>%
    dplyr::select(entrezgene, !!symbol_name, d_genes, n_pathways,  
                  n_significant, n_significant_mono,
                  mean_expression_1, mean_expression_2, 
                  mean_expr_perc_ant, mean_expr_perc_post,
                  mean_expr_fold_ant, mean_expr_fold_post, t_test, pathways)
  
  
  
  #
  # Save results
  #
  dir_tables <- paste0(dir_output, "tables/", measure, "/")
  if(!dir.exists(dir_tables)) {
    dir.create(paste0(dir_tables, "pathways/"), recursive = TRUE)
    dir.create(paste0(dir_tables, "genes/"), recursive = TRUE)
  }
  
  save_file <- paste0(dir_tables, "pathways/lp", lp, "_", 
                      ifelse(dataset == "facebase", pair_name, "neuro.rds"), 
                      ".txt")
  cat("\t", sum(df_path_signif$t_test <= 0.1, na.rm = TRUE), "of", 
      nrow(df_path_signif), 
      "DC pathways are also significantly DE.\n")
  save_df <- df_path_signif %>% 
    filter((mean_expression_1 > quantile(df_path$mean_expression_1, 0.5)) | 
             (mean_expression_2 > quantile(df_path$mean_expression_2, 0.5)),
           d_pathway > quantile(df_path$d_pathway, 0.5)) %>%
    arrange(t_test > signif, -d_pathway) %>%
    mutate(n_genes = paste0(n_genes, "(", n_genes_significant, ")"),
           d_pathway = round(d_pathway, 2),
           expression_ant = round(mean_expr_fold_ant, 2),
           expression_post = round(mean_expr_fold_post, 2),
           t_test = round_signif(t_test, 3)) %>%
    dplyr::select(pathway, d_pathway, n_genes, expression_ant, expression_post, 
                  t_test)
  cat("... saving results for dc pathways in", pair_name, "to", save_file, "\n")
  write.table(save_df, 
              file = save_file, 
              quote = FALSE, sep = " \t& ", eol = " \\\\\n", row.names = FALSE)
  
  
  save_file <- paste0(dir_tables, "genes/lp", lp, "_", 
                      ifelse(dataset == "facebase", pair_name, "neuro.rds"), 
                      ".txt")
  cat("\t", sum(df_genes_signif$t_test <= 0.1, na.rm = TRUE), "of", 
      nrow(df_genes_signif), 
      "DC genes are also significantly DE.\n")
  save_df <- df_genes_signif %>% 
    filter((mean_expression_1 > quantile(df_genes$mean_expression_1, 0.5)) | 
             (mean_expression_2 > quantile(df_genes$mean_expression_2, 0.5)),
           d_genes > quantile(df_genes$d_genes, 0.5)) %>%
    arrange(t_test > signif, -d_genes) %>%
    mutate(n_pathways = paste0(n_pathways, "(", n_significant_mono, ")"),
           d_genes = round(d_genes, 2),
           expression_ant = round(mean_expr_fold_ant, 2),
           expression_post = round(mean_expr_fold_post, 2),
           # expression_ant = paste0(round(mean_expr_fold_ant, 2), 
           #                         "(", 100 * round(mean_expr_perc_ant, 2), ")"),
           # expression_post = paste0(round(mean_expr_fold_post, 2), 
           #                         "(", 100 * round(mean_expr_perc_post, 2), ")"),
           # mean_expression = round(mean_expression, 2),
           t_test = round_signif(t_test, 3)) %>%
    dplyr::select(!!symbol_name, d_genes, n_pathways, expression_ant, 
                  expression_post, t_test)
  cat("... saving results for dc genes to", save_file, "\n")
  write.table(save_df, 
              file = save_file, 
              quote = FALSE, sep = " \t& ", eol = " \\\\\n", row.names = FALSE)
  
  end_line <- paste0(paste(rep("=", nchar(pair_name) + 1), collapse = ""),
                     "======================================")
  cat(end_line, "\n")
}
