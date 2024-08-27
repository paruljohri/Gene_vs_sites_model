#!/bin/bash
#SBATCH --mail-user=pjohri@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -p general
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-2:00
#SBATCH --mem=20m
#SBATCH -o /nas/longleaf/home/pjohri/LOGFILES/LDstats%j.out
#SBATCH -e /nas/longleaf/home/pjohri/LOGFILES/LDstats%j.err

#module load python/3.9.6

num_sites=1000 #1000/100/2000
sites="thousand_site" #thousand_site/hundred_site/humans
theta="theta0_005" #theta0_005/theta0_00025
fitness_model="additive_site" #additive_site/additive_gene
epistatic="noEpistasis"
infolder="/work/users/p/j/pjohri/FitnessNote/simulations/${sites}/${fitness_model}/${theta}/${epistatic}"
outfolder="/work/users/p/j/pjohri/FitnessNote/statistics/${sites}/${theta}/${epistatic}"
ListOfGamma=( "gamma100") #gamma2/gamma20/gamma2_mean/gamma20_mean/gamma100_mean/gamma1000_mean
ListOfDom=( "h_0_0" "h_0_2" "h_0_5")

cd /nas/longleaf/home/pjohri/SlimStats
for gamma in "${ListOfGamma[@]}"
do
	for dom in "${ListOfDom[@]}"
	do
		echo "calculating stats for " ${gamma} " and " ${dom}
		#using my own script:
		python statistics_LD_general_reps.py -fileExt ".ms" -sampleSize 100 -regionLen ${num_sites} -input_folder ${infolder}/${gamma}/${dom} -output_folder ${outfolder} -output_prefix ${fitness_model}_${gamma}_${dom}
	done
done


echo "Finished"


