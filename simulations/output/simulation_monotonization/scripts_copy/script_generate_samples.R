if(!exists("generator_params")) {
  stop("Variable `generator_params` needs to be specified.")
}
if(!exists("seed_network")) {
  stop("Variable `seed_network` needs to be specified.")
}
if(!exists("output_dir")) {
  stop("Variable `output_dir` needs to be specified.")
}

save_file <- paste0(output_dir, "generated_samples/", 
                    "gen_seed_", seed_network, ".rds")
if(file.exists(save_file)) {
  stop(paste("generated samples already exist for this network seed at", save_file))
}

network_list <- get_networks()

# Generate samples.
cat("Generating sample data from network 1 and 2 with\n",
    "\t n =", generator_params$n, "\n")
gen_list <- gen_gaussian(generator_params$n, network_list = network_list)
cat("... saving pair of generated samples to", save_file, "\n")
saveRDS(gen_list, save_file)

rm(network_list, gen_list)
gc()
