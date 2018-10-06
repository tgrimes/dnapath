init_cran_libraries <- function() {
  cat("Initializing CRAN libraries...\n")
  #Required files:
  packages <- c("corpcor",
                "dplyr",
                "tibble",
                "ggplot2",
                "Rcpp",
                "RcppArmadillo")

  #Install any missing packages.
  if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
  }
  print(sapply(packages, require, character.only = TRUE))
}

init_github_libraries <- function() {
  cat("Initializing GitHub libraries...\n")
  github_packages <- c("SeqNet")
  print(sapply(github_packages, require, character.only = TRUE))
}

init_local_R_files <- function(script_dir) {
  cat("Initializing local R files in", script_dir, "...\n")
  file_names <- c("dna.R", 
                  "miscellaneous.R")
  success <- rep(FALSE, length(file_names))
  names(success) <- c(file_names)
  
  #Load scripts from source files.
  for(file in file_names) {
    success[names(success) == file] <- 
    tryCatch({
        source(paste0(script_dir, file))
        TRUE
      }, error = function(e) {
        FALSE
      })
  }
  print(success)
}

init_local_cpp_files <- function(script_dir) {
  cat("Initializing local C++ files in", script_dir, "...\n")
  file_names <- c("dnaC.cpp")
  success <- rep(FALSE, length(file_names))
  names(success) <- c(file_names)
  
  #Load scripts from source files.
  for(file in file_names) {
    success[names(success) == file] <- 
      tryCatch({
        sourceCpp(paste0(script_dir, file))
        TRUE
      }, error = function(e) {
        FALSE
      })
  }
  print(success)
}


init_cran_libraries()
init_github_libraries()
init_local_R_files(script_dir)
init_local_cpp_files(script_dir)