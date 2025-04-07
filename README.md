# Gene_vs_sites_model

Here are the scripts used to perform analyses in the manuscripts entitled "A gene-based model of fitness and its implications for genetic variation: Linkage disequilibrium" (https://www.biorxiv.org/content/10.1101/2024.09.12.612686v2.abstract) and "A gene-based model of fitness and its implications for genetic variation: Genetics and inbreeding loads" (https://www.biorxiv.org/content/10.1101/2025.02.19.639162v1.abstract) by Parul Johri and Brian Charlesworth.\ 

Please direct any questions to either author. Contact information: pjohri@unc.edu; Brian.Charlesworth@ed.ac.uk.

Description of scripts:
>>1_two_site_simulations:\
>>BC_ComputerCode_1.f90\
Program for Properties of Deterministic Equilibrium Populations:
    This program uses iteration to calculate the equilibrium haplotype frequencies and load statistics for the two-locus additive gene and sites models with no recombination, allowing for epistasis.\
>>
>>BC_ComputerCode_2.f90\
>>Program for Properties of Equilibrium Multi-site Gene Model with no epistasis:
This program uses the approximations in the last section of the Appendix of the LD paper to calculate the properties of the gene model with multiple sites.\
>>
>>BC_ComputerCode_3.f90\
>>Program for Properties of Multiplicative Fitness Model:
This program uses iteration to calculate the equilibrium haplotype frequencies and load 
statistics for the two-locus gene and sites multiplicative models with no recombination, allowing for 
epistasis.\
>>
>>BC_ComputerCode_4.f90\
>>Program for Simulating Two Segregating Loci:
This program simulates a haploid population with two segregating sites with no   
recombination, allowing weights to be applied towards low allele frequencies by the 
methods of Garcia & Lohmueller (2021) or Good (2022).\
>>
>>BC_ComputerCode_5.f90\
>>Program for Simulating Samples from a Population with Two Segregating Loci:
This program simulates repeatied sampling from a population with two loci with the same allele frequency at each locus 
and a specified value of  D. The full set of D values for a fixed number of replicates can be stored.\

>>2_multiple_site_simulations:\
SLiM scripts that were used to run the simulations are provided. Bash scripts that contain the command lines used to run the SLiM scripts are also provided.
Note that load statistics were calculated within SLiM and stored in output files. 

>>3_calculation_of_statistics_multiple_site_simulations:\
Python scripts are provided that were used to calculate mean frequencies of selected alleles, LD between selected and minor alleles, and load statatistics.

>>4_make_tables:\
R scripts were used to obtain the final mean and SE values of LD and load statistics across replicates. These are provided here. 

>>5_make_figures:\
R scripts used to make the final figures are provided.
