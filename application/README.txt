Integrating gene regulatory pathways into differential co-expression analysis of gene expression data

Application to craniofacial and neuroblastoma datasets
------------------------------------------------------

Project orginization:
  data/ 
  	annotations/
  		- Stores the various annotations to entrezgenes that were obtained from BioMart
  	facebase/
  		counts/ : .csv files containing the gene expression counts
  		gene_length/ : the length of each entrezgene, used in the normalization of counts.
  		normalized/ : .csv files containing normalized expression counts.
  		raw/ : the raw BAM files are available on the FaceBase website.
  	neuroblastoma/
  		counts/ : .csv files containing normalized counts for the two tumor groups.
  		raw/ : the raw gene expression files are available in the GEO database.
  
  logs/ : contains logs from running each step in the analysis for both datasets.
  
  output/ : contains .txt files with the summary tables for each analysis.
  
  src/
    - main_facebase.R : running this file will rerun the entire analysis for the craniofacial dataset.
    - main_neuroblastoma.R : running this file performs the analysis for the neuroblastoma dataset.
    - dna.R : contains the main functions for performing the differential network analysis.
    - dnaC.cpp : this Rcpp/Armadillo code contains an implementation of the sparse covariance estimation, along with some other helper functions for the differential network analysis. This code was compiled on a macOS 10.13.6 machine. 
  
