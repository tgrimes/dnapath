timestamp <- function() {
  format(Sys.time(), "%Y_%m_%d_%Hh_%Mm_%Ss")
}

open_log <- function(log_file_name = NULL, 
                     log_file_dir = NULL,
                     add_timestamp = TRUE,
                     append = TRUE,
                     type = "output",
                     split = TRUE) {
  if(!is.null(log_file_name))
    log_file_name <- gsub(".txt", "", log_file_name)
  
  if(is.null(log_file_name)) {
    log_file <- paste0(log_file_dir, timestamp(), ".txt")
  } else if(add_timestamp) {
    log_file <- paste0(log_file_dir, log_file_name, "_", timestamp(), ".txt")
  } else {
    log_file <- paste0(log_file_dir, log_file_name, ".txt")
  }
  sink(log_file, append = append, type = type, split = split)
  cat("Log opened", format(Sys.time()), "\n\n")
}

close_log <- function() {
  cat("\n\nLog closed", format(Sys.time()))
  cat("\n\n***\nNotes:\n")
  sink()
}

source_dir <- function(dir) {
  for (name in list.files(dir)) {
    source(file.path(dir, name))
  }
}

#' Run a script
#'
#' Opens and closes log to record the details of the script.
#' @param script_name The name of the script to run.
run_script <- function(script_name) {
  success <- NULL
  
  open_log(paste0(output_dir, "log.txt"))
  cat("Running script:", script_name, "\n")
  success <- tryCatch({
    source(paste0(script_dir, script_name, ".R"), local = TRUE)
    TRUE
  }, error = function(e) {
    print(e)
    FALSE
  })
  close_log()
  
  gc()
  
  if(!success) print("Warning: script produced error.")
  rm(success)
}

save_df_as_latex_table <- function(df, save_file, caption = NULL, digits = 3) {
  my_round <- function(x, digits) {
    x <- sapply(x, function(x) {
      if(is.numeric(x)) return(round(x, digits))
      return(x)
    })
  }
  p <- ncol(df)
  n <- nrow(df)
  file.create(save_file)
  sink(save_file, append = TRUE)
  cat("\\begin{table}[h!]\n")
  cat("\\caption{")
  if(!is.null(caption)) cat(caption)
  cat("}\n")
  cat(paste0("\\begin{tabular}{l", 
             paste(rep("c", p - 1), collapse = ""), "}\n"))
  cat("\\hline\n")
  cat(paste0(paste(colnames(df), collapse = "\t&\t"), "\t\\\\\n"))
  cat("\\hline\n")
  for(i in 1:n) {
    cat(paste0(paste(my_round(df[i, ], digits), collapse = "\t&\t"), 
               "\t\\\\\n"))
  }
  cat("\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
  sink()
}

round_signif <- function(x, digits = 3) {
  if (length(x) > 1) {
    return(sapply(x, round_signif, simplify = "array"))
  }
  if (x == 0) return(0)
  sign_x <- sign(x)
  x <- abs(x)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits - 1)
  sign_x * base * 10^exponent
}