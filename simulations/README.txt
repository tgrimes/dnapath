Integrating gene regulatory pathways into differential co-expression analysis of gene expression data

Simulation study
------------------------------------------------------
Project orginization:
  output/
  	figures/ : contains .png files illustrating the generated differential network and its pathways, and the performance results for various simulation settings.
  	generated_network/ : contains the generated network(s) as an .rds file(s).
  	generated_samples/ : contains the generated gene expression samples as an .rds file.
  	log.txt : log entries from each step in the simulation process.
  	scripts_copy/ : this folder contains a copy of all the scripts used in the simulation.
  	sim_results_dna.txt : main output from the simulation.

  src/
     scripts_main.R - run this script to perform the entire simulation.
     script_analyze_dna.R - used to create various plots of the simulation output.
     
  README :
  differential.Rproj : opens up project into RStudio.
  submit.sh : used to submit jobs to HiPerGator
  
