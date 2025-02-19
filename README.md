# Gene_vs_sites_model
Description of scripts:
>>1_two_site_simulations:\
>>BC_ComputerCode_1.f90\
Program for Properties of Deterministic Equilibrium Populations 
    This program uses exact recursions equations to calculate the equilibrium haplotype frequencies and load 
     statistics for the two-locus gene and sites fitness models with no recombination, allowing for 
     epistasis. It also calculates  approximate measures of LD.\
>>
>>BC_ComputerCode_2.f90\
>>Program for Simulating Two Segregating Loci
   This program simulates a haploid population with two segregating sites with no   
   recombination, allowing weights to be applied towards low allele frequencies by the 
   methods of Garcia & Lohmueller (2021) or Good (2022)\
>>
>>BC_ComputerCode_3.f90\
>>Program for Simulating Samples from a Population with Two Segregating Loci
   This program simulates repeatied sampling from a population with two loci with the same allele frequency at each locus 
    and a specified value of  D. The full set of D values for a fixed number of replicates can be stored.\

>>2_multiple_site_simulations:\
SLiM scripts that were used to run the simulations are provided. Bash scripts that contain the command lines used to run the SLiM scripts are also provided.

>>3_calculation_of_statistics_multiple_site_simulations:\
//get LD across all SNPs and mean allele frequencies for all sites\
-run_calculate_LD_af_stats.sh\
-run_calculate_LD_af_stats_epistasis.sh\
-run_calculate_LD_af_stats_humans.sh\
//get load stats:\
-summarize_load_stats_epistasis.py\
-summarize_load_stats_humans.py\
-summarize_load_stats.py\

>>4_make_tables:\
>>make_table_LD_stats_multiple_sites_epistasis_v2.r\
make_table_LD_stats_multiple_sites_v2_humans.r\
make_table_LD_stats_multiple_sites_v2_multiplicative.r\
make_table_LD_stats_multiple_sites_v2_N5000.r\
make_table_LD_stats_multiple_sites_v2.r\
>>allele frequency and load stats:\
>>make_table_load_stats_multiple_sites_epistasis.r\
make_table_load_stats_multiple_sites_humans.r\
make_table_load_stats_multiple_sites_N5000.r\
make_table_load_stats_multiple_sites.r\

>>5_make_figures:
