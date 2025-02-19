Here is hte description of each script:

Scripts to perform forward simulations using SLiM:
-two_site_dintd_fitness_v4.slim
This is a script to model fitness using only two sites. This script was NOT used to perform any analyses in the paper, but is useful to understand the fitness function for two sites.

-multiple_site_dintd_fitness_LD.slim
Simulates 1000 linked selected sites with an option to calcualte fitness in 3 different fitness models- mutliplicative, sites, or gene model. Assumes a constant and single value of "s".

-multiple_site_dintd_fitness_LD_DFE.slim
Simulates 1000 linked selected sites with an option to calcualte fitness in 3 different fitness models- mutliplicative, sites, or gene model. Assumes that selection coefficients follow a DFE (modelled by a gamma distribution).

-multiple_site_dintd_fitness_LD_DFE_N5000.slim
Simulates 1000 linked selected sites with an option to calcualte fitness in 3 different fitness models- mutliplicative, sites, or gene model. Assumes that selection coefficients follow a DFE (modelled by a gamma distribution). Here we use 5000 diploid individuals to simulate the population instead of 1000 diploid individuals (as assumed in other scripts).

-multiple_site_dintd_fitness_LD_DFE_epistasis_v3.slim
Simulates 1000 linked selected sites with an option to calcualte fitness in 3 different fitness models- mutliplicative, sites, or gene model. Assumes that selection coefficients follow a DFE (modelled by a gamma distribution). Includes the effect of epistasis, as described in the paper (see Methods).

-multiple_site_dintd_fitness_LD_DFE_epistasis_v3_N5000.slim
Simulates 1000 linked selected sites with an option to calcualte fitness in 3 different fitness models- mutliplicative, sites, or gene model. Assumes that selection coefficients follow a DFE (modelled by a gamma distribution). Includes the effect of epistasis, as described in the paper (see Methods). Here we use 5000 diploid individuals to simulate the population instead of 1000 diploid individuals (as assumed in other scripts).

-run_multiple_site_dintd_fitness_constant_s.sh
Command lines used to run simulations with fixed selective effect of deleterious mutations at all sites.

-run_multiple_site_dintd_fitness_with_DFE.sh
Command lines used to run simulations that incorporated a DFE. There is no epistasis here.

-run_multiple_site_dintd_fitness_with_DFE_N5000.sh
Command lines used to run simulations that incorporated a DFE and a lower rescaling factor.

-run_multiple_site_dintd_fitness_with_DFE_gamma100_epistasis.sh
Command lines used to run simulations that incorporated a DFE and epistasis.

-run_multiple_site_dintd_fitness_with_DFE_gamma100_epistasis_lowrec.sh
Command lines used to run simulations that incorporated a DFE and epistasis for very low rates of recombination (0.1x mean).

-run_multiple_site_dintd_fitness_with_DFE_humans.sh
Command lines used to run simulations of human-like parameters. There is a DFE but no epistasis.

