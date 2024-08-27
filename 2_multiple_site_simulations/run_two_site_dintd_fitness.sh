#!/bin/bash
#SBATCH --mail-user=pjohri@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -p general
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-0:20
#SBATCH --mem=20m
#SBATCH -a 1-1000%100
#SBATCH -o /nas/longleaf/home/pjohri/LOGFILES/fitness_%A_rep%a.out
#SBATCH -e /nas/longleaf/home/pjohri/LOGFILES/fitness_%A_rep%a.err

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#To be run on LONGLEAF!
#10 simulations per submitted job!
#######################################


declare -i repID=0+$SLURM_ARRAY_TASK_ID

module load slim

cd /nas/longleaf/home/pjohri/FitnessNote/programs

echo "starting simulation " $repID

fitness_model="multiplicative" #multiplicative/additive_site/additive_gene
folder="/work/users/p/j/pjohri/FitnessNote/simulations/${fitness_model}"

##gamma=-2; h=0.5
slim -d d_seed=${repID} -d "d_fitness_model='${fitness_model}'" -d d_gamma=-2.0 -d d_dom=0.5 -d "d_folder='${folder}/gamma2/h_0_5'" -d "d_repID='${repID}'" two_site_dintd_fitness_v4.slim

##gamma=-2; h=0.2
slim -d d_seed=${repID} -d "d_fitness_model='${fitness_model}'" -d d_gamma=-2.0 -d d_dom=0.2 -d "d_folder='${folder}/gamma2/h_0_2'" -d "d_repID='${repID}'" two_site_dintd_fitness_v4.slim

##gamma=-2; h=0.0
slim -d d_seed=${repID} -d "d_fitness_model='${fitness_model}'" -d d_gamma=-2.0 -d d_dom=0.0 -d "d_folder='${folder}/gamma2/h_0_0'" -d "d_repID='${repID}'" two_site_dintd_fitness_v4.slim

##gamma=-20; h=0.5
slim -d d_seed=${repID} -d "d_fitness_model='${fitness_model}'" -d d_gamma=-20.0 -d d_dom=0.5 -d "d_folder='${folder}/gamma20/h_0_5'" -d "d_repID='${repID}'" two_site_dintd_fitness_v4.slim

##gamma=-20; h=0.2
slim -d d_seed=${repID} -d "d_fitness_model='${fitness_model}'" -d d_gamma=-20.0 -d d_dom=0.2 -d "d_folder='${folder}/gamma20/h_0_2'" -d "d_repID='${repID}'" two_site_dintd_fitness_v4.slim

##gamma=-20; h=0.0
slim -d d_seed=${repID} -d "d_fitness_model='${fitness_model}'" -d d_gamma=-20.0 -d d_dom=0.0 -d "d_folder='${folder}/gamma20/h_0_0'" -d "d_repID='${repID}'" two_site_dintd_fitness_v4.slim

echo "Finished simulation " $repID





