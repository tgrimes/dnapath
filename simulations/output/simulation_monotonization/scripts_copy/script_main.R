# Directions:
# For each simulation run, do the following:
#   1. Specify `folder` name below. This is where the simulation will be stored.
#   2. Go to `script_create_networks.R` and modify networks, if desired.
#   3. Modify `generator_params` below; these are used in the SeqNet generator.
#   4. Set `n_set` and `n_iter_set` to desired sample size(s) and simulation size.
#   5. If re-running previous scripts, change script_dir (see below before
#      helper functions)
#   6. Run this script.

folder <- "simulation_monotonization" # Name of folder containing desired simulation.
output_dir <- "output/" # Output directory containing simulation folders.
script_dir <- "src/" # Directory containing simulation scripts.

# If the simulation folder does not exist, create it.
if(!(folder %in% list.files(output_dir)) |
   length(list.files(paste0(output_dir, folder))) == 0) {
  dir.create(paste0(output_dir, folder))
  file.create(paste0(output_dir, folder, "/log.txt"))
  dir.create(paste0(output_dir, folder, "/figures"))
  dir.create(paste0(output_dir, folder, "/generated_networks"))
  dir.create(paste0(output_dir, folder, "/generated_samples"))
  dir.create(paste0(output_dir, folder, "/scripts_copy"))
  dir.create(paste0(output_dir, folder, "/tables"))
  
  
  # Copy scripts used into the results folder.
  list_of_scripts <- list.files(script_dir, full.names = TRUE)
  file.copy(list_of_scripts, paste0(output_dir, folder, "/scripts_copy/"))
  
  cat("Created new directory:", paste0(output_dir, folder, "/:"), "\n",
      paste0(folder, "/"), "\n",
      "|- figures/\n",
      "|- generated_networks/\n",
      "|- generated_samples/\n",
      "|- log.txt\n",
      "|- scripts_copy/\n",
      "|- tables/\n")
}

# Move output director to chosen folder.
output_dir <- paste0(output_dir, folder, "/")

# # Note, if scripts are saved, may want to use the following for replication:
script_dir <- paste0(output_dir, "scripts_copy/")

# Source required R files and load packages.
source(paste0(script_dir, "init_required_packages_and_files.R"))

#######################################
#
# Helper functions
#
#######################################
get_networks <- function() {
  load_file <- paste0(output_dir, "generated_networks/", 
                      "network_seed_", seed_network, ".rds")
  cat("... loading list of networks from", load_file, "\n")
  network_list <- readRDS(load_file)
  
  return(list(`1` = network_list[[1]], `2` = network_list[[2]]))
}

get_generated_data <- function() {
  load_file <- paste0(output_dir, "generated_samples/", 
                      "gen_seed_", seed_network, ".rds")
  cat("... loading list of generated samples from", load_file, "\n")
  gen_list <- readRDS(load_file)
  
  return(list(`1` = gen_list[[1]], `2` = gen_list[[2]]))
}

#######################################
#
# Simulation parameters
#
#######################################
iter_max <- 100
n_max <- 1000
n_iter_set <- 51:74 # Index for repeated simulations with fixed settings.
n_set <- c(50, 100, 250, 500, 750, 1000) # max(n_set) * max(n_iter_set) cannot exceed generated n.
network_iter_set <- 1 # < 1000; Number of networks to simulate per network setup. 
network_params_list <- list(
  small = list(p = 100, n_pathways = 4, mu = 20, sigma = 10, name = "small"),
  medium = list(p = 250, n_pathways = 10, mu = 20, sigma = 10, name = "medium"),
  large = list(p = 500, n_pathways = 20, mu = 20, sigma = 10, name = "large")
)
lp_set <- seq(0.5, 3.5, 0.25)
measure_set <- list(
  corr = function(x) {
    scores <- corC(x)
    colnames(scores) <- colnames(x)
    return(scores)
  },
  pcor = function(x) {
    scores <- pcor_shrinkC(x)
    colnames(scores) <- colnames(x)
    return(scores)
  }) # Association measures.
signal_set <- c(NA, 1, 0.9, 0.8) # Proportion of correct nodes in pathway in dna.
n_perm <- 100 # Number of permutations to perform in dna.
signif_set <- c(0.01, 0.05) 
generator_params <- list(
  n = n_max * iter_max # Create enough samples for entire simulation.
)
results_file <- paste0(output_dir, "sim_results_dna.txt")

#######################################
#
# Run simulation
#
#######################################
n_simulations <- max(network_iter_set) * max(n_iter_set) * length(n_set) * 1.4 # Adjust for signal reps.
times = list(small = c(1, 2),
             medium = c(4, 11),
             large = c(20, 50))
est <- lapply(times, function(time) sum(time * n_simulations) / 60 / 60) # in hours
marks <- c(clisymbols::symbol$cross,
           clisymbols::symbol$tick)[names(times) %in% names(network_params_list) + 1]
cat("Total simulations to perform:", n_simulations, 
    "\nEstimated run time:",
    "\n\t", marks[1], names(times)[1], ":", round(est[[1]], 1), "hours",
    "\n\t", marks[2], names(times)[2], ":", round(est[[2]], 1), "hours",
    "\n\t", marks[3], names(times)[3], ":", round(est[[3]], 1), "hours\n")

for(i in 3) {
  network_params <- network_params_list[[i]]
  for(network_number in network_iter_set) {
    seed_network <- 10000000 + 1000 * i + network_number
    
    # Create network pair.
    run_script("script_create_networks")
    run_script("script_plot_networks")
    
    # Generate samples for both networks.
    run_script("script_generate_samples")
    
    # Perform differential network analysis on samples.
    run_script("script_run_dna")
  }
}


#######################################
#
# Analyze the results
#
#######################################
# run_script("script_analyze_dna")
