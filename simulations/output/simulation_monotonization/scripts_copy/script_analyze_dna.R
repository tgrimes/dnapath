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


#####################################
#
# Performance results
#
#####################################
df1 <- results %>%
  dplyr::mutate(Sensitivity = TP / (TP + FN),
                Specificity = TN / (TN + FP),
                TDR = ifelse(TP + FP == 0, 0, TP / (TP + FP)),
                TNDR = ifelse(TN + FN == 0, 0, TN / (TN + FN)),
                F1 = 2 / (1 / Sensitivity + 1 / TDR),
                MCC = (TP * TN - FP * FN) / 
                  exp(0.5 * (log(TP + FP) + log(TP + FN) + 
                               log(TN + FP) + log(TN + FN))),
                MCC = ifelse(is.nan(MCC), 0, MCC)) %>%
  tidyr::gather(key = "performance", value = "value", 
                Sensitivity, Specificity, TDR, TNDR, F1, MCC) %>%
  mutate(performance = factor(performance, 
                              levels = c("Sensitivity", "Specificity", 
                                         "F1", "TDR", "TNDR", "MCC")))


p <- 500
lp <- 2
signif <- 0.05
df <- filter(df1, p == !!p, signif == !!signif, lp == !!lp)
# Make figure.
m <- unique(results$measure)[2]
comp <- unique(results$component)[2]
mono <- unique(results$monotonized)[1]
for(m in unique(results$measure)) {
  for(comp in unique(results$component)) {
    for(mono in unique(results$monotonized)) {
      
      if(!(mono == "TRUE" && comp == "paths")) {
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
}

g <- df %>%
  filter(component == comp, measure == m, monotonized == mono) %>%
  filter(!(performance %in% c("Specificity", "TNDR"))) %>%
  mutate(performance = factor(performance, 
                              levels = c("Sensitivity", "TDR", "F1", "MCC"))) %>%
  ggplot(aes(x = as.factor(n), y = value, color = signal)) +
  facet_grid(. ~ performance) + 
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

save_file <- paste0(output_dir, "figures/sim_results_", m, "_", 
                    comp, "_p", p, "_lp", lp, "_alpha", 100*signif, 
                    ifelse(mono == "TRUE", "_mono", ""), ".png")
cat("... Saving performance figures to", save_file, "\n")
ggsave(save_file, g, scale = 0.38, width = 600, height = 160, 
       units = "mm", device = "png")



df <- df1 %>%
  filter(p == 500, 
         signif == 0.05, 
         lp == 2)
df %>%
  filter(n == 250,
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




p <- 500
signif <- 0.05
df <- filter(df1, p == !!p, signif == !!signif)
n <- 1000
m <- unique(results$measure)[2]
comp <- unique(results$component)[2]
mono <- unique(results$monotonized)[2]
# Make figure.
for(n in c(250, 1000)) {
  for(m in unique(results$measure)[2]) {
    for(comp in unique(results$component)[1]) {
      for(mono in unique(results$monotonized)[2]) {
        
        if(!(mono == "TRUE" && comp == "paths")) {
          g <- df %>%
            filter(component == comp, measure == m, monotonized == mono,
                   n == !!n, signal != "NA") %>%
            ggplot(aes(x = lp, y = value, color = signal)) +
            facet_wrap(~ performance) + 
            geom_line(aes(group = paste0(sample_number, signal), color = signal), 
                      alpha = 0.03) +
            stat_summary(aes(group = signal), fun.y = median, geom = "line", 
                         size = 0.75) +
            theme_bw() +
            labs(y = NULL, x = "Value of p in the test statistic", 
                 color = "Pathway\nknowledge",
                 title = paste0("Performance n = ", n, " (", m, ", ", comp, ") ",
                                "p = ", p, ", alpha = ", signif,
                                ifelse(mono == "TRUE", " - Monotonized", ""))) +
            lims(y = c(0, 1)) +
            scale_x_continuous(breaks = 1:3,
                               minor_breaks = setdiff(unique(df$lp), 1:3)) +
            theme(legend.position = "right")
          
          print(g)
          
          save_file <- paste0(output_dir, "figures/by_lp_performance_", m, "_", 
                              comp, "_p", p, "_n", n, "_alpha", 100*signif, 
                              ifelse(mono == "TRUE", "_mono", ""), ".png")
          cat("... Saving performance figures to", save_file, "\n")
          ggsave(save_file, g, scale = 0.38, width = 540, height = 320, 
                 units = "mm", device = "png")
        }
      }
    }
  }
}

