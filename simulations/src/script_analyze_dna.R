if(!exists("output_dir")) {
  stop("Variable `output_dir` does not exist in global environment.")
}
if(!exists("results_file")) {
  stop("Variable `results_file` does not exist in global environment.")
}
cat("... loading results from", results_file, "\n")
results <- as.tibble(read.table(results_file, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE))
results <- results %>%
  mutate(signal = factor(signal, levels = c("1", "0.9", "0.8", "NA"))) %>%
  dplyr::select(-time, -perm)

# sample_number
# network_number
# header <- c("time", "sample_number", "network_number", "network_seed", "component", 
#             "n", "p", "lp", "signal", "perm", "signif", "measure",
#             "sens", "spec", "tdr", "tndr", "f1")
perf_results <- results %>%
  dplyr::mutate(sens = TP / (TP + FN),
                spec = TN / (TN + FP),
                tdr = ifelse(TP + FP == 0, 0, TP / (TP + FP)),
                tndr = ifelse(TN + FN == 0, 0, TN / (TN + FN)),
                f1 = 2 / (1 / sens + 1 / tdr),
                MCC = (TP * TN - FP * FN) / exp(0.5 *
                                                  (log(TP + FP) + log(TP + FN) + log(TN + FP) + log(TN + FN))),
                MCC = ifelse(is.nan(MCC), 0, MCC)) %>%
  tidyr::gather(key = "performance", value = "value", 
                sens, spec, tdr, tndr, f1, MCC) %>%
  mutate(performance = factor(performance, 
                              levels = c("sens", "spec", "f1", 
                                         "tdr", "tndr", "MCC"))) %>%
  group_by(component, n, p, signal, signif, measure, performance, lp, monotonized) %>%
  summarise(lower = quantile(value, 0.25, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            upper = quantile(value, 0.75, na.rm = TRUE)) %>%
  ungroup()


#####################################
#
# Performance results
#
#####################################

p <- 500
signif <- 0.01
n <- 50
df <- filter(perf_results, p == !!p, signif == !!signif, n == !!n)
# Make figure.
m <- unique(perf_results$measure)[1]
comp <- unique(perf_results$component)[2]
mono <- unique(perf_results$monotonized)[2]
for(m in unique(perf_results$measure)) {
  for(comp in unique(perf_results$component)) {
    for(mono in unique(perf_results$monotonized)) {
      g <- df %>%
        filter(component == comp, measure == m, monotonized == mono) %>%
        ggplot(aes(x = lp, y = median, color = signal)) +
        facet_wrap(~ performance) + 
        geom_line(size = 1) +
        geom_point(size = 1, alpha = 0.5) +
        geom_errorbar(aes(ymin = lower, ymax = upper),
                      width = 0.05) +
        theme_bw() +
        labs(y = NULL, x = "Lp", color = "Pathway\nknowledge",
             title = paste0("Performance n", n, " (", m, ", ", comp, ") ",
                            "p = ", p, ", alpha = ", signif)) +
        lims(y = c(0, 1), x = c(0.5, 3.5)) +
        theme(legend.position = "bottom")
      
      print(g)
      
      save_file <- paste0(output_dir, "figures/performance_", m, "_", 
                          comp, "_p", p, "_lp", lp, "_alpha", signif, ".png")
      cat("... Saving performance figures to", save_file, "\n")
      ggsave(save_file, g, scale = 0.38, width = 540, height = 360, 
             units = "mm", device = "png")
    }
  }
}



p <- 500
lp <- 2
signif <- 0.01
df <- filter(perf_results, p == !!p, signif == !!signif, lp == !!lp)
# Make figure.
m <- unique(perf_results$measure)[1]
comp <- unique(perf_results$component)[2]
mono <- unique(perf_results$monotonized)[2]
for(m in unique(perf_results$measure)) {
  for(comp in unique(perf_results$component)) {
    for(mono in unique(perf_results$monotonized)) {
      
      if(mono == "TRUE" && comp == "paths") next
      
      g <- df %>%
        filter(component == comp, measure == m, monotonized == mono) %>%
        ggplot(aes(x = n, y = median, color = signal)) +
        facet_wrap(~ performance) + 
        geom_line(size = 1) +
        # geom_point(position = position_dodge(width = 1), size = 1, alpha = 0.5) +
        geom_errorbar(aes(ymin = lower, ymax = upper),
                      position = position_dodge(width = 1),
                      width = 75) +
        theme_bw() +
        labs(y = NULL, x = "Sample size", color = "Pathway\nknowledge",
             title = paste0("Performance Lp", lp, " (", m, ", ", comp, ") ",
                            "p = ", p, ", alpha = ", signif,
                            ifelse(mono == "TRUE", " - Monotonized", ""))) +
        lims(y = c(0, 1)) +
        theme(legend.position = "bottom")
      
      print(g)
      
      save_file <- paste0(output_dir, "figures/performance_", m, "_", 
                          comp, "_p", p, "_lp", lp, "_alpha", signif, 
                          ifelse(mono == "TRUE", "_mono", ""), ".png")
      cat("... Saving performance figures to", save_file, "\n")
      ggsave(save_file, g, scale = 0.38, width = 540, height = 360, 
             units = "mm", device = "png")
    }
  }
}

# # Make table.
# caption <- ""
# caption <- gsub("\n", "", caption)
# perf_results %>%
#   dplyr::select(-user, -system, -elapsed, -pathway_signal, -value, -lower, -upper) %>%
#   filter(sample_number == 1, n == 100) %>%
#   tidyr::spread(key = measure, value = median) %>%
#   dplyr::select(-sample_number, -network_number, -perm) %>%
#   save_df_as_latex_table(paste0(output_dir, "/tables/performance.txt"),
#                          caption, digits = 3)




p <- 500
lp <- 2
signif <- 0.01
df1 <- results %>%
  dplyr::mutate(sens = TP / (TP + FN),
                spec = TN / (TN + FP),
                tdr = ifelse(TP + FP == 0, 0, TP / (TP + FP)),
                tndr = ifelse(TN + FN == 0, 0, TN / (TN + FN)),
                f1 = 2 / (1 / sens + 1 / tdr),
                MCC = (TP * TN - FP * FN) / exp(0.5 *
                                                  (log(TP + FP) + log(TP + FN) + log(TN + FP) + log(TN + FN))),
                MCC = ifelse(is.nan(MCC), 0, MCC)) %>%
  tidyr::gather(key = "performance", value = "value", 
                sens, spec, tdr, tndr, f1, MCC) %>%
  mutate(performance = factor(performance, 
                              levels = c("sens", "spec", "f1", 
                                         "tdr", "tndr", "MCC")))
df <- filter(df1, p == !!p, signif == !!signif, lp == !!lp)
# Make figure.
m <- unique(perf_results$measure)[1]
comp <- unique(perf_results$component)[2]
mono <- unique(perf_results$monotonized)[1]
for(m in unique(perf_results$measure)) {
  for(comp in unique(perf_results$component)) {
    for(mono in unique(perf_results$monotonized)) {
      
      if(mono == "TRUE" && comp == "paths") next
      
      g <- df %>%
        filter(component == comp, measure == m, monotonized == mono) %>%
        ggplot(aes(x = as.factor(n), y = value, color = signal)) +
        facet_wrap(~ performance) + 
        geom_boxplot(size = 0.5, outlier.size = 0.2) +
        stat_summary(aes(group = signal), fun.y = median, geom = "line", 
                     size = 0.5) +
        theme_bw() +
        labs(y = NULL, x = "Sample size", color = "Pathway\nknowledge",
             title = paste0("Performance Lp", lp, " (", m, ", ", comp, ") ",
                            "p = ", p, ", alpha = ", signif,
                            ifelse(mono == "TRUE", " - Monotonized", ""))) +
        lims(y = c(0, 1)) +
        theme(legend.position = "right")
      
      print(g)
      
      save_file <- paste0(output_dir, "figures/performance_", m, "_", 
                          comp, "_p", p, "_lp", lp, "_alpha", 100*signif, 
                          ifelse(mono == "TRUE", "_mono", ""), ".png")
      cat("... Saving performance figures to", save_file, "\n")
      ggsave(save_file, g, scale = 0.38, width = 540, height = 320, 
             units = "mm", device = "png")
    }
  }
}


df <- df1 %>%
  filter(p == 500, 
         signif == 0.01, 
         lp == 2)


df %>%
  filter(n == 500,
         monotonized == "TRUE",
         measure == "pcor",
         component == "edges") %>%
  arrange(sample_number) %>%
  group_by(signal, performance) %>%
  mutate(m = cummean(value)) %>%
  ungroup() %>%
  ggplot(aes(x = sample_number, y = m, color = signal)) +
  facet_wrap(.~performance) +
  geom_line() +
  geom_line(aes(x = sample_number, y = value, color = signal), alpha = 0.2) +
  theme_bw()
