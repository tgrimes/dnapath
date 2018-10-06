#' Obtain counts from raw data
#' 
#' @param csv_file a string containing path and file name of csv file to load.
#' @return list containing counts for both networks.
create_counts <- function(csv_file) {
  
  cat("Loading .csv file from", csv_file, "\n")
  counts <- as_tibble(read.csv(csv_file, header = TRUE))
  
  attr(counts, "source") <- csv_file
  
  return(counts)
}

#' Filter genes using coefficient of variation
#' 
#' Genes with a coefficient of variation less than the threshold in any sample 
#' are removed from all samples.
#' @param threshold threshold used for low variability filter.
#' @return a list containing updated counts for both networks, and information
#' about this filter appended.
apply_filter_cv_at_least <- function(counts, threshold) {
  cat("\t - Filtering counts: cv at least", threshold, "\n")
  filters <- attr(counts, "filters")
  fn_cv <- function(x) ifelse(mean(x) == 0, 0, sd(x)/mean(x))
  
  p <- nrow(counts)
  counts <- counts %>%
    mutate(cv = dplyr::select(., starts_with("sample")) %>% apply(1, fn_cv)) %>%
    filter(cv >= threshold) %>%
    dplyr::select(-cv)
  cat("\t \t removed", p - nrow(counts), "genes with cv below", threshold, "\n")
  
  attr(counts, "filters") <- c(filters, list(cv = threshold))
  
  return(counts)
}

#' Filter genes with too many repeated counts
#' 
#' Genes with the same counts repeated in many observations are removed: If the 
#' maximum number of repeats is at or above the threshold in any one sample, 
#' remove that gene from all samples
#' @param threshold threshold used for low variability filter.
#' @return a list containing updated counts for both networks, and information
#' about this filter appended.
apply_filter_ties_at_most <- function(counts, threshold) {
  cat("\t - Filtering counts: ties at most", threshold, "\n")
  filters <- attr(counts, "filters")
  fn_ties <- function(x) max(table(x))
  
  p <- nrow(counts)
  counts <- counts %>%
    mutate(ties = dplyr::select(., starts_with("sample")) %>% apply(1, fn_ties)) %>%
    filter(ties <= threshold) %>%
    dplyr::select(-ties)
  cat("\t \t removed", p - nrow(counts), "genes with more than", threshold, "ties\n")
  
  attr(counts, "filters") <- c(filters, list(ties = threshold))
  
  return(counts)
}

#' Filter genes with too many repeated counts
#' 
#' Genes with too many zeroes are removed: If the percentage of zeros is
#' at or above the threshold, remove that gene.
#' @param threshold the maximum percentage of zero allowed.
#' @return a list containing updated counts for both networks, and information
#' about this filter appended.
apply_filter_percent_zero_at_most <- function(counts, threshold) {
  cat("\t - Filtering counts: percent zero at most", threshold, "\n")
  filters <- attr(counts, "filters")
  fn_zeroes <- function(x) mean(x == 0)
  
  p <- nrow(counts)
  counts <- counts %>%
    mutate(zeroes = dplyr::select(., -starts_with("sample")) %>% apply(1, fn_zeroes)) %>%
    filter(zeroes <= threshold) %>%
    dplyr::select(-zeroes)
  cat("\t \t removed", p - nrow(counts), "genes with more than", threshold, "zeroes\n")
  
  attr(counts, "filters") <- c(filters, list(zeroes = threshold))
  
  return(counts)
}

#' Filter genes with small range across samples
#' 
#' Remove genes if their counts do not vary much across all samples: If the 
#' range of counts across samples is less than the threshold, remove the gene.
#' @param threshold threshold used for low variability filter.
#' @return a list containing updated counts for both networks, and information
#' about this filter appended.
apply_filter_range_at_least <- function(counts, threshold) {
  cat("\t - Filtering counts: range at least", threshold, "\n")
  filters <- attr(counts, "filters")
  fn_range <- function(x) max(x) - min(x)
  
  p <- nrow(counts)
  counts <- counts %>%
    mutate(range = dplyr::select(., starts_with("sample")) %>% apply(1, fn_range)) %>%
    filter(range >= threshold) %>%
    dplyr::select(-range)
  cat("\t \t removed", p - nrow(counts), "genes with range below", threshold, "\n")
  
  attr(counts, "filters") <- c(filters, list(range = threshold))
  
  return(counts)
}


#' Filter genes with small counts
#' 
#' Sets a minimum threshold for gene counts. Any gene whose maximum count is 
#' less than this threshold are removed.
#' @param threshold gene must have at least one count above threshold
#' @return a list containing updated counts for both networks, and information
#' about this filter appended.
apply_filter_max_at_least <- function(counts, threshold) {
  cat("\t - Filtering counts: max at least", threshold, "\n")
  filters <- attr(counts, "filters")
  
  p <- nrow(counts)
  counts <- counts %>%
    mutate(max = dplyr::select(., starts_with("sample")) %>% apply(1, max)) %>%
    filter(max >= threshold) %>%
    dplyr::select(-max)
  cat("\t \t removed", p - nrow(counts), "genes with max below", threshold, "\n")
  
  attr(counts, "filters") <- c(filters, list(max = threshold))
  
  return(counts)
}


#' Filter genes with small counts
#' 
#' Sets a minimum threshold for gene counts. Any genes whose minimum count is 
#' less than this threshold are removed.
#' @param threshold gene must have all counts at or above threshold.
#' @return a list containing updated counts for both networks, and information
#' about this filter appended.
apply_filter_min_at_least <- function(counts, threshold) {
  cat("\t - Filtering counts: min at least", threshold, "\n")
  filters <- attr(counts, "filters")
  
  p <- nrow(counts)
  counts <- counts %>%
    mutate(min = dplyr::select(., starts_with("sample")) %>% apply(1, min)) %>%
    filter(min >= threshold) %>%
    dplyr::select(-min)
  cat("\t \t removed", p - nrow(counts), "genes with min below", threshold, "\n")
  
  attr(counts, "filters") <- c(filters, list(min = threshold))
  
  return(counts)
}


#' Save counts object 
#' 
#' @param counts The counts to be saved.
#' @param file_name The name of the file. If NULL, a name is generated based
#' on the filters used on the counts.
#' @param output_dir The directory where the counts are saved.
#' @param make_subdir Should a subdirecty in output_dir be made for the 
#' filtered data? If TRUE, the subdirecty name will be the filtered data name.
#' @return A string containing the name of the file saved.
save_counts <- function(counts,
                        dir_output = getwd(),
                        file_name = NULL,
                        make_subdir = TRUE) {
  if(is.null(file_name)) {
    if(is.data.frame(counts)) {
      filters <- attributes(counts)$filters
    } else {
      filters <- attributes(counts[[1]])$filters
    }
    n_filters <- length(filters)
    if(n_filters > 0) {
      filters <- sapply(filters, c)
      file_name <- paste(paste0(names(filters), filters), collapse = "_")
    } else {
      file_name <- "no_filters"
    }
  }
  
  if(make_subdir) {
    dir_output <- paste0(dir_output, file_name, "/")
  }
  
  if(!dir.exists(dir_output)) {
    dir.create(dir_output, recursive = TRUE)
  }
  
  save_file <- paste0(dir_output, file_name)
  cat("Saving filtered counts to", save_file, "\n")
  saveRDS(counts, save_file)
  
  return(file_name)
}


#' Load counts object 
#' 
#' @param dir The directory where the counts are stored.
#' @param file_name The name of the counts file to load. 
load_counts <- function(dir, file_name) {
  counts <- readRDS(file = paste0(dir, file_name))
  
  return(counts)
}
