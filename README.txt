This repository contain the scripts used in performing the analysis of the scientific reports manuscript entitled:
"Dynamics of SARS-CoV-2 mutations reveals regional-specificity and similar trends of N501 and High-Frequency mutation N501Y in different levels of control measures"
Justo et al. 2021

Description of files:
Preparation/ = This folder contain the python scripts and files necessary to perform the filtrations of low quality genomes from the GISAID database and the alignment against the reference.
Frequency_calculation/ = This folder contains the python script to calculate the relative frequency of each nucleotide in each position of the genome
Retrieve_genomes_mutated/ = This folder contains the python script to take sequences from an alignment with an specific mutation
master_script.R = This R script contains all the analysis performed in the manuscript, files generated with the python scripts of the other folders are necessary to perform the analysis